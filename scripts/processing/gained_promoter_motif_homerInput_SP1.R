# ##############################################################################
# filename:    gained_promoter_motif_homerInput_SP1.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Test whether gained loop anchor promoters are enriched for
#              SP/KLF motifs relative to expression-matched background
#              promoters. Focus = promoters overlapping gained loop anchors;
#              background = expression-matched promoters via nullranges::matchRanges()
#              using rowMeans of DESeq2 normalized counts across all timecourse
#              samples as the matching covariate.
# ##############################################################################

# Libraries ----
library(mariner)
library(InteractionSet)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(nullranges)
library(plyranges)
library(DESeq2)
library(ggplot2)

# Parameters ----
proj_dir       <- "/work/users/j/p/jpflores/projects/STRS"
difflps_rds    <- file.path(proj_dir, "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
rnaseq_rds     <- file.path(proj_dir, "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
out_dir        <- file.path(proj_dir, "data/processed/cutntag/homer_input")
plot_dir       <- file.path(proj_dir, "plots")
tss_upstream   <- 2000L
tss_downstream <-  500L
match_method   <- "stratified"  # preferred over "rejection" when pool/focal ratio is small
match_replace  <- FALSE

dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,  showWarnings = FALSE, recursive = TRUE)

# Data import ----
cat("Loading differential loops...\n")
diffLoops <- readRDS(difflps_rds)

gainedLoops <- diffLoops[
  which(rowData(diffLoops)$padj < 0.1 &
          rowData(diffLoops)$log2FoldChange > 0)
]
cat("  Gained loops:", length(gainedLoops), "\n")

# Analysis ----

## Extract gained loop anchors ----

## Flatten both anchors; reduce() merges overlapping anchors so no region
## is double-counted in the focal set
gained_anchors <- c(
  anchors(gainedLoops, "first"),
  anchors(gainedLoops, "second")
) |>
  sort() |>
  GenomicRanges::reduce()

cat("  Unique gained anchors (post-reduce):", length(gained_anchors), "\n")

## Build promoter GRanges from TxDb ----
cat("Building promoter GRanges from TxDb (hg38)...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## promoters() is strand-aware
all_tx <- transcripts(txdb, columns = c("tx_id", "gene_id"))

all_promoters <- promoters(
  all_tx,
  upstream   = tss_upstream,
  downstream = tss_downstream
)

## gene_id returns as CharacterList from TxDb; flatten to plain character
all_promoters$entrez_id <- as.character(all_promoters$gene_id)

## Map Entrez -> gene symbol
entrez_to_symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(all_promoters$entrez_id),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENTREZID"
)

all_promoters$gene_symbol <- entrez_to_symbol$SYMBOL[
  match(all_promoters$entrez_id, entrez_to_symbol$ENTREZID)
]

## Drop alt/random/unplaced contigs
all_promoters <- keepStandardChromosomes(all_promoters, pruning.mode = "coarse")

cat("  Total promoter regions:", length(all_promoters), "\n")

## Derive per-gene mean expression across the timecourse ----
##
## The LRT DESeqDataSet has no single meaningful baseMean from results().
## Instead: normalized counts (counts(dds, normalized=TRUE)) rowMeans across
## all samples — one number per gene reflecting overall expression level
## timepoint-agnostically.
cat("Loading RNA-seq LRT timecourse DESeqDataSet...\n")
dds <- readRDS(rnaseq_rds)

stopifnot(is(dds, "DESeqDataSet"))
if (!all(sizeFactors(dds) > 0)) {
  stop("Size factors not estimated. Run estimateSizeFactors(dds) first.")
}

norm_counts      <- counts(dds, normalized = TRUE)
timecourse_mean  <- rowMeans(norm_counts)

cat("  Timecourse samples:", ncol(norm_counts), "\n")
cat("  Genes in DESeqDataSet:", nrow(norm_counts), "\n")

rna_df <- data.frame(
  gene_id         = names(timecourse_mean),
  timecourse_mean = timecourse_mean,
  stringsAsFactors = FALSE
)

## Map Ensembl -> Entrez for joining with TxDb-derived promoters
rna_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rna_df$gene_id,
  columns = c("ENSEMBL", "ENTREZID"),
  keytype = "ENSEMBL"
)

rna_df$entrez_id <- rna_map$ENTREZID[match(rna_df$gene_id, rna_map$ENSEMBL)]
rna_df <- rna_df[!is.na(rna_df$entrez_id) & !is.na(rna_df$timecourse_mean), ]

cat("  Genes with Entrez mapping:", nrow(rna_df), "\n")

## Attach timecourse mean to promoter GRanges ----
all_promoters$timecourse_mean <- rna_df$timecourse_mean[
  match(all_promoters$entrez_id, rna_df$entrez_id)
]
all_promoters <- all_promoters[!is.na(all_promoters$timecourse_mean)]

cat("  Promoters with expression data:", length(all_promoters), "\n")

## Split focal vs pool ----
cat("Overlapping promoters with gained loop anchors...\n")

ol        <- findOverlaps(all_promoters, gained_anchors)
focus_idx <- unique(queryHits(ol))
pool_idx  <- setdiff(seq_along(all_promoters), focus_idx)

focus_promoters <- all_promoters[focus_idx]
pool_promoters  <- all_promoters[pool_idx]

cat("  Focal promoters (at gained anchors):", length(focus_promoters), "\n")
cat("  Pool  promoters (not at gained anchors):", length(pool_promoters), "\n")
cat("  Pool/focal ratio:", round(length(pool_promoters) / length(focus_promoters), 1), "\n\n")

if (length(pool_promoters) < 5 * length(focus_promoters)) {
  warning("Pool is less than 5x focal size. Consider replace = TRUE.")
}

## matchRanges: expression-matched background ----
##
## Fits logistic regression (focal vs pool ~ timecourse_mean) for propensity
## scores, then uses stratified sampling to mirror the focal expression
## distribution in the selected background set.
cat("Running matchRanges (method =", match_method,
    "| replace =", match_replace, ")...\n")

set.seed(42)
mgr <- matchRanges(
  focal   = focus_promoters,
  pool    = pool_promoters,
  covar   = ~ timecourse_mean,
  method  = match_method,
  replace = match_replace
)

matched_bg <- matched(mgr)
cat("  Matched background regions:", length(matched_bg), "\n\n")

## Quality assessment ----
cat("Matching quality (overview):\n")
print(overview(mgr))

# Visualization ----
p_propensity <- plotPropensity(mgr)
ggsave(
  file.path(plot_dir, "matchRanges_qc_propensity_SP1.pdf"),
  p_propensity, width = 6, height = 4
)
cat("  Propensity score QC plot saved.\n")

p_covariate <- plotCovariate(mgr)
ggsave(
  file.path(plot_dir, "matchRanges_qc_covariate_SP1.pdf"),
  p_covariate, width = 6, height = 4
)
cat("  Covariate distribution QC plot saved.\n\n")

## Format BED files for HOMER ----
## HOMER BED: 0-based, 6 columns (chr, start, end, name, score, strand)
ranges_to_bed <- function(gr, label_prefix) {
  data.frame(
    chr    = as.character(seqnames(gr)),
    start  = start(gr) - 1L,
    end    = end(gr),
    name   = paste0(label_prefix, "_", seq_along(gr)),
    score  = 0,
    strand = as.character(strand(gr))
  )
}

focus_bed <- ranges_to_bed(focus_promoters, "focus")
bg_bed    <- ranges_to_bed(matched_bg,      "background")

# Save outputs ----
focus_out <- file.path(out_dir, "gained_anchor_promoters_focus_SP1.bed")
bg_out    <- file.path(out_dir, "gained_anchor_promoters_background_SP1.bed")

write.table(focus_bed, focus_out,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(bg_bed, bg_out,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("BED files written:\n")
cat("  Focus:      ", focus_out, "\n")
cat("  Background: ", bg_out, "\n")
cat("\nFocus n =", nrow(focus_bed), "| Background n =", nrow(bg_bed), "\n")

sessionInfo()

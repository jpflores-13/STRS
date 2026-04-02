## gained_promoter_motif_homerInput_SP1.R
##
## Goal: Test whether gained loop anchor promoters are enriched for SP/KLF motifs
##       relative to expression-matched background promoters.
##
## Strategy:
##   Focus      -> promoters (TSS ± upstream/downstream) overlapping gained loop anchors
##   Background -> promoters of genes NOT at gained anchors, matched on mean
##                 expression across the timecourse using nullranges::matchRanges()
##
## Expression covariate:
##   The RNA-seq object is a timecourse DESeq2 dataset fit with LRT.
##   There is no single "baseMean" from results() that fairly represents
##   expression across all timepoints. Instead, we compute rowMeans() on
##   DESeq2 normalized counts (size-factor corrected) across all samples,
##   giving one per-gene mean that reflects overall expression level
##   regardless of timepoint or treatment.
##
## Output:
##   homer_input/gained_anchor_promoters_focus_SP1.bed
##   homer_input/gained_anchor_promoters_background_SP1.bed
##   plots/matchRanges_qc_propensity_SP1.pdf
##   plots/matchRanges_qc_covariate_SP1.pdf

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

# ── Parameters ────────────────────────────────────────────────────────────────

PROJ_DIR    <- "/work/users/j/p/jpflores/projects/STRS"
DIFFLPS_RDS <- file.path(PROJ_DIR, "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
RNASEQ_RDS  <- file.path(PROJ_DIR, "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
OUT_DIR     <- file.path(PROJ_DIR, "data/processed/cutntag/homer_input")
PLOT_DIR    <- file.path(PROJ_DIR, "plots")

## Promoter window: how far upstream/downstream of TSS to extend
TSS_UPSTREAM   <- 2000L
TSS_DOWNSTREAM <-  500L

## matchRanges settings
## "stratified" without replacement is preferred over "rejection" when the
## pool is not orders of magnitude larger than focal (Davis et al. 2023)
MATCH_METHOD  <- "stratified"
MATCH_REPLACE <- FALSE

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load differential loops ────────────────────────────────────────────────

cat("Loading differential loops...\n")
diffLoops <- readRDS(DIFFLPS_RDS)

## Gained loops: padj < 0.1 and positive log2FC (sorbitol > control)
gainedLoops <- diffLoops[
  which(rowData(diffLoops)$padj < 0.1 &
          rowData(diffLoops)$log2FoldChange > 0)
]

cat("  Gained loops:", length(gainedLoops), "\n")

# ── 2. Extract gained loop anchors ────────────────────────────────────────────

## Flatten both anchors into one GRanges; reduce() merges overlapping anchors
## so the same genomic region isn't double-counted in the focal set
gained_anchors <- c(
  anchors(gainedLoops, "first"),
  anchors(gainedLoops, "second")
) |>
  sort() |>
  GenomicRanges::reduce()

cat("  Unique gained anchors (post-reduce):", length(gained_anchors), "\n")

# ── 3. Build promoter GRanges from TxDb ──────────────────────────────────────

cat("Building promoter GRanges from TxDb (hg38)...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## promoters() is strand-aware: upstream/downstream are interpreted relative
## to TSS direction, so (+) and (-) strand genes are handled correctly
all_tx <- transcripts(txdb, columns = c("tx_id", "gene_id"))

all_promoters <- promoters(
  all_tx,
  upstream   = TSS_UPSTREAM,
  downstream = TSS_DOWNSTREAM
)

## gene_id returns as a CharacterList from TxDb; flatten to plain character
all_promoters$entrez_id <- as.character(all_promoters$gene_id)

## Map Entrez -> gene symbol for readability
entrez_to_symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(all_promoters$entrez_id),
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENTREZID"
)

all_promoters$gene_symbol <- entrez_to_symbol$SYMBOL[
  match(all_promoters$entrez_id, entrez_to_symbol$ENTREZID)
]

## Drop alt/random/unplaced contigs — HOMER doesn't handle them well
all_promoters <- keepStandardChromosomes(all_promoters, pruning.mode = "coarse")

cat("  Total promoter regions:", length(all_promoters), "\n")

# ── 4. Derive per-gene mean expression across the timecourse ─────────────────
##
## The LRT timecourse DESeqDataSet has no single meaningful baseMean from
## results() — that would reflect only one contrast. Instead:
##
##   1. Extract DESeq2 normalized counts (counts(dds, normalized = TRUE)):
##      these are raw counts divided by sample-specific size factors,
##      making samples comparable despite differences in sequencing depth.
##   2. Take rowMeans() across ALL samples (all timepoints, all replicates).
##      This gives one number per gene: its average normalized expression
##      across the entire experiment, timepoint-agnostic.
##
## This is the fairest expression summary for matching because:
##   - It captures steady-state expression level, not a timepoint-specific value
##   - It mirrors how DESeq2 itself computes baseMean internally (same formula)
##   - Genes with similar rowMeans had similar overall activity, making them
##     valid comparators regardless of whether they changed over the timecourse

cat("Loading RNA-seq LRT timecourse DESeqDataSet...\n")
dds <- readRDS(RNASEQ_RDS)

## Sanity check: confirm this is a DESeqDataSet with size factors estimated
stopifnot(is(dds, "DESeqDataSet"))
if (!all(sizeFactors(dds) > 0)) {
  stop("Size factors not estimated. Run estimateSizeFactors(dds) first.")
}

## Normalized counts matrix: genes x samples
norm_counts <- counts(dds, normalized = TRUE)

cat("  Timecourse samples:", ncol(norm_counts), "\n")
cat("  Genes in DESeqDataSet:", nrow(norm_counts), "\n")

## One mean per gene across all samples
timecourse_mean <- rowMeans(norm_counts)

## Build a tidy data frame for joining; rownames are Ensembl IDs
rna_df <- data.frame(
  gene_id          = names(timecourse_mean),
  timecourse_mean  = timecourse_mean,
  stringsAsFactors = FALSE
)

## Map Ensembl -> Entrez to join with TxDb-derived promoters
rna_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rna_df$gene_id,
  columns = c("ENSEMBL", "ENTREZID"),
  keytype = "ENSEMBL"
)

rna_df$entrez_id <- rna_map$ENTREZID[match(rna_df$gene_id, rna_map$ENSEMBL)]
rna_df <- rna_df[!is.na(rna_df$entrez_id) & !is.na(rna_df$timecourse_mean), ]

cat("  Genes with Entrez mapping:", nrow(rna_df), "\n")

# ── 5. Attach timecourse mean to promoter GRanges ────────────────────────────

all_promoters$timecourse_mean <- rna_df$timecourse_mean[
  match(all_promoters$entrez_id, rna_df$entrez_id)
]

all_promoters <- all_promoters[!is.na(all_promoters$timecourse_mean)]

cat("  Promoters with expression data:", length(all_promoters), "\n")

# ── 6. Split into focal vs pool ───────────────────────────────────────────────

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

# ── 7. matchRanges: expression-matched background ────────────────────────────
##
## matchRanges() fits a logistic regression (focal vs pool ~ timecourse_mean)
## to get propensity scores, then uses stratified sampling to select pool
## ranges whose timecourse_mean distribution mirrors the focal set.

cat("Running matchRanges (method =", MATCH_METHOD,
    "| replace =", MATCH_REPLACE, ")...\n")

set.seed(42)
mgr <- matchRanges(
  focal   = focus_promoters,
  pool    = pool_promoters,
  covar   = ~ timecourse_mean,
  method  = MATCH_METHOD,
  replace = MATCH_REPLACE
)

## matched() extracts the matched GRanges from the MatchedGRanges object
## equivalent to: pool(mgr)[indices(mgr, set = "matched")]
matched_bg <- matched(mgr)

cat("  Matched background regions:", length(matched_bg), "\n\n")

# ── 8. Quality assessment ─────────────────────────────────────────────────────
##
## overview() reports mean/SD of covariates and propensity scores for focal,
## matched, pool, and unmatched sets, plus the focal-matched difference.
## The focal-matched difference for timecourse_mean should be close to zero.

cat("Matching quality (overview):\n")
print(overview(mgr))

p_propensity <- plotPropensity(mgr)
ggsave(
  file.path(PLOT_DIR, "matchRanges_qc_propensity_SP1.pdf"),
  p_propensity,
  width = 6, height = 4
)
cat("\n  Propensity score QC plot saved.\n")

p_covariate <- plotCovariate(mgr)
ggsave(
  file.path(PLOT_DIR, "matchRanges_qc_covariate_SP1.pdf"),
  p_covariate,
  width = 6, height = 4
)
cat("  Covariate distribution QC plot saved.\n\n")

# ── 9. Format BED files for HOMER ─────────────────────────────────────────────

## HOMER BED: 0-based, 6 columns (chr, start, end, name, score, strand)
## GRanges are 1-based, so subtract 1 from start

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

# ── 10. Write output ──────────────────────────────────────────────────────────

focus_out <- file.path(OUT_DIR, "gained_anchor_promoters_focus_SP1.bed")
bg_out    <- file.path(OUT_DIR, "gained_anchor_promoters_background_SP1.bed")

write.table(focus_bed, focus_out,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

write.table(bg_bed, bg_out,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

cat("BED files written:\n")
cat("  Focus:      ", focus_out, "\n")
cat("  Background: ", bg_out, "\n")
cat("\nFocus n =", nrow(focus_bed), "| Background n =", nrow(bg_bed), "\n")
cat("\nDone! Run homer_gained_promoter_motifs_SP1.sh next.\n")
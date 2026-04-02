# ##############################################################################
# filename:    extract_gained_anchor_promoters_homerInput.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Identify gene promoters (TSS ± 2000/500 bp) overlapping gained
#              loop anchors; export BED files and gene lists for HOMER analysis
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Parameters ----
diff_loops_rds       <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
homer_output_dir     <- "data/processed/cutntag/homer_input"
gained_promoters_bed <- file.path(homer_output_dir, "gained_anchor_promoters.bed")
all_promoters_bed    <- file.path(homer_output_dir, "all_promoters.bed")
gained_genes_txt     <- file.path(homer_output_dir, "gained_anchor_genes.txt")
tss_upstream         <- 2000
tss_downstream       <- 500

# Data import ----
cat("Loading differential loops...\n")
diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions()

# Analysis ----

## Extract gained loops
gained_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 &
                                        diff_loopCounts$log2FoldChange > 0)]
cat("Number of gained loops:", length(gained_loops), "\n")

## Get all loop anchors from gained loops
gained_anchors <- c(anchors(gained_loops, type = "first"),
                    anchors(gained_loops, type = "second")) |>
  unique()

seqlevelsStyle(gained_anchors) <- "UCSC"
cat("Number of unique gained anchors:", length(gained_anchors), "\n")

## Define promoters by genes
cat("Defining gene promoters...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr    <- genes(txdb)
promoters_gr <- promoters(genes_gr, upstream = tss_upstream, downstream = tss_downstream)

seqlevelsStyle(promoters_gr) <- "UCSC"
promoters_gr$gene_id <- names(promoters_gr)

## Keep only standard chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
promoters_gr   <- promoters_gr[seqnames(promoters_gr) %in% standard_chrs]
gained_anchors <- gained_anchors[seqnames(gained_anchors) %in% standard_chrs]

cat("Promoters on standard chromosomes:", length(promoters_gr), "\n")
cat("Gained anchors on standard chromosomes:", length(gained_anchors), "\n")

## Find promoters overlapping gained anchors
cat("Finding promoters at gained anchors...\n")
overlaps              <- findOverlaps(promoters_gr, gained_anchors)
gained_promoter_indices <- unique(queryHits(overlaps))

gained_promoters <- promoters_gr[gained_promoter_indices]
all_promoters    <- promoters_gr
cat("Gained anchor promoters found:", length(gained_promoters), "\n")

## Convert to data frames with unique IDs
cat("\nPreparing BED files...\n")
gained_df <- as.data.frame(gained_promoters) |>
  dplyr::select(seqnames, start, end, gene_id, strand) |>
  mutate(
    start    = start - 1,  # 0-based for BED
    score    = 0,
    seqnames = as.character(seqnames),
    strand   = as.character(strand),
    name     = paste0(seqnames, ":", start, "-", end)
  )

all_df <- as.data.frame(all_promoters) |>
  dplyr::select(seqnames, start, end, gene_id, strand) |>
  mutate(
    start    = start - 1,
    score    = 0,
    seqnames = as.character(seqnames),
    strand   = as.character(strand),
    name     = paste0(seqnames, ":", start, "-", end)
  )

## Remove duplicates
if (any(duplicated(gained_df$name))) {
  warning("Duplicates found in gained promoters! Removing...")
  gained_df <- gained_df[!duplicated(gained_df$name), ]
}

if (any(duplicated(all_df$name))) {
  warning("Duplicates found in all promoters! Removing...")
  all_df <- all_df[!duplicated(all_df$name), ]
}

## Add gene symbols
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = gained_promoters$gene_id,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")
gained_promoters$symbol <- gene_symbols$SYMBOL[match(gained_promoters$gene_id,
                                                     gene_symbols$ENTREZID)]

# Save outputs ----
dir.create(homer_output_dir, showWarnings = FALSE, recursive = TRUE)

## BED6 format: chr, start, end, name, score, strand
write.table(gained_df[, c("seqnames", "start", "end", "name", "score", "strand")],
            gained_promoters_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(all_df[, c("seqnames", "start", "end", "name", "score", "strand")],
            all_promoters_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(data.frame(gene_id = gained_promoters$gene_id,
                       symbol  = gained_promoters$symbol),
            gained_genes_txt,
            quote = FALSE, sep = "\t", row.names = FALSE)

cat("\nSummary\n")
cat("Total promoters:", length(all_promoters), "\n")
cat("Gained anchor promoters:", length(gained_promoters), "\n")
cat("Percentage:", round(100 * length(gained_promoters) / length(all_promoters), 2), "%\n")
cat("\nOutput files written to", homer_output_dir, "\n")

sessionInfo()

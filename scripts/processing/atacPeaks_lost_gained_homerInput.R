# ##############################################################################
# filename:    atacPeaks_lost_gained_homerInput.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Identify ATAC peaks overlapping gained and lost loop anchors
#              and export BED files for HOMER motif analysis
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(InteractionSet)

# Parameters ----
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
atac_peaks_file   <- "data/processed/atac/output/peaks/YAPP_HEK_cont_0h_peaks.narrowPeak"
homer_output_dir  <- "data/processed/atac/homer_input"
gained_atac_bed   <- file.path(homer_output_dir, "gained_anchor_atac.bed")
lost_atac_bed     <- file.path(homer_output_dir, "lost_anchor_atac.bed")

# Data import ----

## Load differential loops
diff_loopCounts <- readRDS(diff_loops_rds) |>
  interactions()

# Analysis ----

## Extract GAINED loops
gained_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 &
                                        diff_loopCounts$log2FoldChange > 0)]

## Extract LOST loops
lost_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 &
                                      diff_loopCounts$log2FoldChange < 0)]

## Get all loop anchors from gained loops
gained_anchors <- c(anchors(gained_loops, type = "first"),
                    anchors(gained_loops, type = "second")) |>
  unique()

# CRITICAL: Ensure chromosome names have 'chr' prefix
seqlevelsStyle(gained_anchors) <- "UCSC"

## Get all loop anchors from lost loops
lost_anchors <- c(anchors(lost_loops, type = "first"),
                  anchors(lost_loops, type = "second")) |>
  unique()

# CRITICAL: Ensure chromosome names have 'chr' prefix
seqlevelsStyle(lost_anchors) <- "UCSC"

## Load ATAC peaks
atac_peaks <- plyranges::read_narrowpeaks(atac_peaks_file)

# CRITICAL: Ensure chromosome names match
seqlevelsStyle(atac_peaks) <- "UCSC"

## Keep only standard chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
atac_peaks    <- atac_peaks[seqnames(atac_peaks) %in% standard_chrs]
gained_anchors <- gained_anchors[seqnames(gained_anchors) %in% standard_chrs]
lost_anchors   <- lost_anchors[seqnames(lost_anchors) %in% standard_chrs]

## Find ATAC peaks overlapping gained anchors
gained_overlaps     <- findOverlaps(atac_peaks, gained_anchors)
gained_atac_indices <- unique(queryHits(gained_overlaps))
gained_atac         <- atac_peaks[gained_atac_indices]

## Find ATAC peaks overlapping lost anchors
lost_overlaps     <- findOverlaps(atac_peaks, lost_anchors)
lost_atac_indices <- unique(queryHits(lost_overlaps))
lost_atac         <- atac_peaks[lost_atac_indices]

## Convert to data frames for BED format
## BED format: chr, start (0-based), end, name, score, strand
gained_df <- as.data.frame(gained_atac) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(lost_atac) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

## Check for duplicates
cat("\nChecking for duplicate IDs...\n")
cat("Gained anchor ATAC peaks - unique IDs:", length(unique(gained_df$name)),
    "out of", nrow(gained_df), "\n")
cat("Lost anchor ATAC peaks - unique IDs:", length(unique(lost_df$name)),
    "out of", nrow(lost_df), "\n")

if (any(duplicated(gained_df$name))) {
  warning("Duplicates found in gained anchor ATAC peaks! Removing...")
  gained_df <- gained_df[!duplicated(gained_df$name), ]
}

if (any(duplicated(lost_df$name))) {
  warning("Duplicates found in lost anchor ATAC peaks! Removing...")
  lost_df <- lost_df[!duplicated(lost_df$name), ]
}

# Save outputs ----
dir.create(homer_output_dir, showWarnings = FALSE, recursive = TRUE)

write.table(gained_df, gained_atac_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, lost_atac_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## Verification
gained_check <- read.table(gained_atac_bed, header = FALSE, stringsAsFactors = FALSE)
lost_check   <- read.table(lost_atac_bed,   header = FALSE, stringsAsFactors = FALSE)

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("Total ATAC peaks:", length(atac_peaks), "\n")
cat("ATAC peaks at gained anchors:", nrow(gained_check), "\n")
cat("ATAC peaks at lost anchors:", nrow(lost_check), "\n")
cat("Gained percentage:", round(100 * nrow(gained_check) / length(atac_peaks), 2), "%\n")
cat("Lost percentage:", round(100 * nrow(lost_check) / length(atac_peaks), 2), "%\n")
cat("\nOutput files written to", homer_output_dir, "\n")
cat("  - gained_anchor_atac.bed\n")
cat("  - lost_anchor_atac.bed\n")

sessionInfo()

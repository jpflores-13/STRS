# ##############################################################################
# filename:    ctcfPeaks_lost_gained_homerInput.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Extract retained vs lost CTCF peaks and export BED files
#              for HOMER motif analysis
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(dplyr)

# Parameters ----
ctcf_control_narrowpeak  <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
ctcf_sorbitol_narrowpeak <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak"
homer_output_dir         <- "data/processed/cutntag/homer_input"
ctcf_retained_bed        <- file.path(homer_output_dir, "ctcf_retained_peaks.bed")
ctcf_lost_bed            <- file.path(homer_output_dir, "ctcf_lost_peaks.bed")

# Data import ----
ctcf_control  <- plyranges::read_narrowpeaks(ctcf_control_narrowpeak)
ctcf_sorbitol <- plyranges::read_narrowpeaks(ctcf_sorbitol_narrowpeak)

# Analysis ----

## Retained peaks: control peaks that overlap with sorbitol peaks
overlaps        <- findOverlaps(ctcf_control, ctcf_sorbitol)
retained_indices <- unique(queryHits(overlaps))

ctcf_retained <- ctcf_control[retained_indices]
ctcf_lost     <- ctcf_control[-retained_indices]

## Convert to data frames for BED format
## BED format: chr, start (0-based), end, name, score, strand
retained_df <- as.data.frame(ctcf_retained) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(ctcf_lost) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

# Save outputs ----
dir.create(homer_output_dir, showWarnings = FALSE, recursive = TRUE)

write.table(retained_df, ctcf_retained_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, ctcf_lost_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("Output files written to", homer_output_dir, "\n")
cat("  - ctcf_retained_peaks.bed\n")
cat("  - ctcf_lost_peaks.bed\n")

sessionInfo()

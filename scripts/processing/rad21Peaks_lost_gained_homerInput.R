# ##############################################################################
# filename:    rad21Peaks_lost_gained_homerInput.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Extract retained vs lost RAD21 peaks and export BED files
#              for HOMER motif analysis
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(dplyr)

# Parameters ----
rad21_control_narrowpeak  <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
rad21_sorbitol_narrowpeak <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_peaks.narrowPeak"
homer_output_dir          <- "data/processed/cutntag/homer_input"
rad21_retained_bed        <- file.path(homer_output_dir, "rad21_retained_peaks.bed")
rad21_lost_bed            <- file.path(homer_output_dir, "rad21_lost_peaks.bed")

# Data import ----
rad21_control  <- plyranges::read_narrowpeaks(rad21_control_narrowpeak)
rad21_sorbitol <- plyranges::read_narrowpeaks(rad21_sorbitol_narrowpeak)

# Analysis ----

## Retained peaks: control peaks that overlap with sorbitol peaks
overlaps        <- findOverlaps(rad21_control, rad21_sorbitol)
retained_indices <- unique(queryHits(overlaps))

rad21_retained <- rad21_control[retained_indices]
rad21_lost     <- rad21_control[-retained_indices]

## Convert to data frames for BED format
## BED format: chr, start (0-based), end, name, score, strand
retained_df <- as.data.frame(rad21_retained) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(rad21_lost) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

# Save outputs ----
dir.create(homer_output_dir, showWarnings = FALSE, recursive = TRUE)

write.table(retained_df, rad21_retained_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, rad21_lost_bed,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat("Output files written to", homer_output_dir, "\n")
cat("  - rad21_retained_peaks.bed\n")
cat("  - rad21_lost_peaks.bed\n")

sessionInfo()

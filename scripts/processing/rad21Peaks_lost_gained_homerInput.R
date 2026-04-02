## Extract retained vs lost RAD21 peaks for HOMER motif analysis

library(GenomicRanges)
library(plyranges)
library(dplyr)

# Load RAD21 peaks ---------------------------------------------------------

## Control RAD21 peaks 
rad21_control <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
)

## RAD21 peaks after sorbitol treatment
rad21_sorbitol <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_peaks.narrowPeak"
)

# Identify retained and lost peaks ----------------------------------------

## Retained peaks: control peaks that overlap with sorbitol peaks
## Using findOverlaps to identify which control peaks have overlaps
overlaps <- findOverlaps(rad21_control, rad21_sorbitol)
retained_indices <- unique(queryHits(overlaps))

rad21_retained <- rad21_control[retained_indices]
rad21_lost <- rad21_control[-retained_indices]

# Export as BED files for HOMER -------------------------------------------

dir.create("data/processed/cutntag/homer_input", 
           showWarnings = FALSE, 
           recursive = TRUE)

## Convert to data frames for BED format
## BED format: chr, start (0-based), end, name, score, strand
retained_df <- as.data.frame(rad21_retained) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(rad21_lost) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

## Write BED files
write.table(retained_df, 
            "data/processed/cutntag/homer_input/rad21_retained_peaks.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, 
            "data/processed/cutntag/homer_input/rad21_lost_peaks.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

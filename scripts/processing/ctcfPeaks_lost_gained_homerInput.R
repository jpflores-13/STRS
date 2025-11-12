## Extract retained vs lost CTCF peaks for HOMER motif analysis

library(GenomicRanges)
library(plyranges)
library(dplyr)

# Load CTCF peaks ---------------------------------------------------------

## Control CTCF peaks (0h, pre-treatment)
ctcf_control <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)

## Sorbitol-treated CTCF peaks (you'll need to specify the correct file)
ctcf_sorbitol <- plyranges::read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak"
)

# Identify retained and lost peaks ----------------------------------------

## Retained peaks: control peaks that overlap with sorbitol peaks
## Using findOverlaps to identify which control peaks have overlaps
overlaps <- findOverlaps(ctcf_control, ctcf_sorbitol)
retained_indices <- unique(queryHits(overlaps))

ctcf_retained <- ctcf_control[retained_indices]
ctcf_lost <- ctcf_control[-retained_indices]

# Export as BED files for HOMER -------------------------------------------

dir.create("data/processed/cutntag/homer_input", 
           showWarnings = FALSE, 
           recursive = TRUE)

## Convert to data frames for BED format
## BED format: chr, start (0-based), end, name, score, strand
retained_df <- as.data.frame(ctcf_retained) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(ctcf_lost) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

## Write BED files
write.table(retained_df, 
            "data/processed/cutntag/homer_input/ctcf_retained_peaks.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, 
            "data/processed/cutntag/homer_input/ctcf_lost_peaks.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

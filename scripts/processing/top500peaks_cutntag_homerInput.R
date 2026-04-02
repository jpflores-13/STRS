# ##############################################################################
# filename:    top500peaks_cutntag_homerInput.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Select top 500 CUT&Tag peaks by signal value for each protein
#              and condition; export BED files for HOMER motif analysis
# ##############################################################################

# Libraries ----
library(plyranges)
library(dplyr)

# Parameters ----
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
homer_output_dir  <- "data/processed/cutntag/homer_input"
target            <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition         <- c("control", "sorbitol")

# Data import ----
cutntag_files <- list.files(cutntag_peaks_dir,
                            full.names = TRUE,
                            pattern = ".narrowPeak")

# Analysis ----
dir.create(homer_output_dir, showWarnings = FALSE, recursive = TRUE)

for (protein in target) {
  for (cond in condition) {

    sample_name <- paste0(protein, "_", cond)
    cat("\nProcessing:", sample_name, "\n")

    # Find the matching file
    pattern   <- paste0(protein, "_", gsub("control", "cont", cond))
    peak_file <- grep(pattern, cutntag_files, value = TRUE, ignore.case = TRUE)

    if (length(peak_file) == 0) {
      cat("  WARNING: No file found for", sample_name, "\n")
      next
    }

    # Read peaks
    peaks <- plyranges::read_narrowpeaks(peak_file[1])

    # Select top 500 by signal value
    top_peaks <- peaks |>
      as.data.frame() |>
      arrange(desc(signalValue)) |>
      head(500) |>
      dplyr::select(seqnames, start, end, name, score, strand) |>
      mutate(start = start - 1)  # Convert to 0-based

    cat("  Selected", nrow(top_peaks), "peaks\n")

    # Write BED file
    output_file <- file.path(homer_output_dir,
                             paste0("top500_", tolower(protein), "_", cond, ".bed"))

    write.table(top_peaks, output_file,
                quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    cat("  Saved to:", output_file, "\n")
  }
}

# Save outputs ----
cat("\nAll BED files created in:", homer_output_dir, "\n")

sessionInfo()

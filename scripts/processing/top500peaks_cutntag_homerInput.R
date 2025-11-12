## Export top 500 peaks for each protein for HOMER motif analysis

library(plyranges)
library(dplyr)

# Load peaks --------------------------------------------------------------

cutntag_files <- list.files("data/processed/cutntag/output/peaks/",
                            full.names = TRUE,
                            pattern = ".narrowPeak")

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")

dir.create("data/processed/cutntag/homer_input", showWarnings = FALSE, recursive = TRUE)

# Process each sample -----------------------------------------------------

for (protein in target) {
  for (cond in condition) {
    
    sample_name <- paste0(protein, "_", cond)
    cat("\nProcessing:", sample_name, "\n")
    
    # Find the matching file
    pattern <- paste0(protein, "_", gsub("control", "cont", cond))
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
    output_file <- paste0("data/processed/cutntag/homer_input/top500_", 
                          tolower(protein), "_", cond, ".bed")
    
    write.table(top_peaks, output_file,
                quote = FALSE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
    
    cat("  Saved to:", output_file, "\n")
  }
}

cat("\n\nAll BED files created in: data/processed/cutntag/homer_input/\n")

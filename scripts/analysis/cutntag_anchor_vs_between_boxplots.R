## Compare protein binding at loop anchors vs between-anchor regions

library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(bamsignals)
library(GenomicRanges)

# Data --------------------------------------------------------------------

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

## Gained loops
gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]

cat("Number of gained loops:", length(gainedLoops), "\n")

## Load CUT&Tag peaks - CORRECTED VERSION
target <- c("CTCF", "RAD21")
condition <- c("control", "sorbitol")

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    # Match by protein name in filename
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(
      "data/processed/cutntag/output/peaks/",
      full.names = TRUE,
      pattern = pattern
    )
    
    # Filter to only .narrowPeak files
    file_path <- file_path[grepl("\\.narrowPeak$", file_path)]
    
    if (length(file_path) == 1) {
      peak_name <- paste0(t, "_", cond)
      peak_list[[peak_name]] <- read_narrowpeaks(file_path)
      cat("Loaded peaks:", peak_name, "from", basename(file_path), "\n")
    } else if (length(file_path) == 0) {
      warning("No peak file found for ", t, "_", cond)
    } else {
      warning("Multiple peak files found for ", t, "_", cond, ": ", paste(basename(file_path), collapse = ", "))
    }
  }
}

## Load BAM files - CORRECTED VERSION
bam_files <- character()
for (t in target) {
  for (cond in condition) {
    # Match by protein name in filename
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(
      "data/processed/cutntag/output/mergeAlign/",
      full.names = TRUE,
      pattern = pattern
    )
    
    # Filter to only .bam files (not .bai)
    file_path <- file_path[grepl("\\.bam$", file_path)]
    
    if (length(file_path) == 1) {
      bam_name <- paste0(t, "_", cond)
      bam_files[bam_name] <- file_path
      cat("Loaded BAM:", bam_name, "from", basename(file_path), "\n")
    } else if (length(file_path) == 0) {
      warning("No BAM file found for ", t, "_", cond)
    } else {
      warning("Multiple BAM files found for ", t, "_", cond, ": ", paste(basename(file_path), collapse = ", "))
    }
  }
}

# Helper Functions --------------------------------------------------------

## Extract loop anchors
get_loop_anchors <- function(loops) {
  first_anchors <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  
  all_anchors <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_id <- c(
    paste0(seq_along(first_anchors), "_anchor1"),
    paste0(seq_along(second_anchors), "_anchor2")
  )
  
  return(all_anchors)
}

## Extract regions BETWEEN loop anchors
get_between_regions <- function(loops) {
  first_anchors <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  
  # Create regions between anchors (from end of first to start of second)
  between_regions <- GRanges(
    seqnames = seqnames(first_anchors),
    ranges = IRanges(
      start = end(first_anchors) + 1,
      end = start(second_anchors) - 1
    ),
    loop_id = paste0(seq_along(first_anchors), "_between")
  )
  
  # Remove any invalid ranges (where anchors overlap or are adjacent)
  between_regions <- between_regions[width(between_regions) > 0]
  
  return(between_regions)
}

## Find overlapping peaks and calculate log2FC
analyze_regions <- function(regions, peaks_control, peaks_treat, 
                            bam_control, bam_treat, region_type) {
  
  # Find peaks overlapping regions
  overlaps <- findOverlaps(peaks_control, regions)
  
  if (length(overlaps) == 0) {
    return(data.frame())
  }
  
  overlapping_peaks <- peaks_control[queryHits(overlaps)]
  mcols(overlapping_peaks)$region_id <- mcols(regions)$loop_id[subjectHits(overlaps)]
  
  # Create unique peak IDs
  peak_ids <- paste0(
    as.character(seqnames(overlapping_peaks)), ":",
    start(overlapping_peaks), "-",
    end(overlapping_peaks)
  )
  
  # Calculate signal from BAM files
  control_signal <- bamCount(bam_control, overlapping_peaks, paired.end = "midpoint")
  treat_signal <- bamCount(bam_treat, overlapping_peaks, paired.end = "midpoint")
  
  # Calculate log2FC
  log2FC <- case_when(
    control_signal == 0 & treat_signal == 0 ~ 0,
    control_signal == 0 ~ log2(treat_signal + 1),
    treat_signal == 0 ~ -log2(control_signal + 1),
    TRUE ~ log2((treat_signal + 1) / (control_signal + 1))
  )
  
  # Create data frame with all information
  result_df <- data.frame(
    region_type = region_type,
    region_id = mcols(overlapping_peaks)$region_id,
    peak_id = peak_ids,
    log2FC = log2FC,
    control_signal = control_signal,
    treat_signal = treat_signal,
    peak_score = mcols(overlapping_peaks)$signalValue,
    stringsAsFactors = FALSE
  )
  
  # Remove duplicate peaks, keeping the one with highest peak score
  result_df <- result_df |>
    group_by(peak_id) |>
    slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    ungroup()
  
  return(result_df)
}

# Analysis ----------------------------------------------------------------

## Get anchors and between regions for gained loops
gained_anchors <- get_loop_anchors(gainedLoops)
gained_between <- get_between_regions(gainedLoops)

cat("\nGained loop regions:\n")
cat("  Anchors:", length(gained_anchors), "\n")
cat("  Between regions:", length(gained_between), "\n")
cat("  Mean between region size:", mean(width(gained_between)), "bp\n")

## Analyze CTCF and RAD21
results_list <- list()

for (protein in target) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing:", protein, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  control_peaks <- peak_list[[paste0(protein, "_control")]]
  sorbitol_peaks <- peak_list[[paste0(protein, "_sorbitol")]]
  control_bam <- bam_files[paste0(protein, "_control")]
  sorbitol_bam <- bam_files[paste0(protein, "_sorbitol")]
  
  # Analyze anchors
  cat("Analyzing loop anchors...\n")
  anchor_results <- analyze_regions(
    gained_anchors, control_peaks, sorbitol_peaks,
    control_bam, sorbitol_bam, "Loop Anchor"
  )
  
  # Analyze between regions
  cat("Analyzing between regions...\n")
  between_results <- analyze_regions(
    gained_between, control_peaks, sorbitol_peaks,
    control_bam, sorbitol_bam, "Between Anchors"
  )
  
  # Combine results
  combined_results <- rbind(anchor_results, between_results) |>
    mutate(protein = protein)
  
  results_list[[protein]] <- combined_results
  
  cat("  Anchor peaks (unique):", nrow(anchor_results), "\n")
  cat("  Between peaks (unique):", nrow(between_results), "\n")
}

## Combine all results
all_results <- bind_rows(results_list)

## Save results
saveRDS(all_results, "data/processed/cutntag/anchor_vs_between_analysis.rds")
write.csv(all_results, "data/processed/cutntag/anchor_vs_between_analysis.csv", 
          row.names = FALSE)

# Visualization -----------------------------------------------------------

## Create box plots
pdf("plots/cutntag_anchor_vs_between_boxplots.pdf", width = 8, height = 6)

for (protein in target) {
  protein_data <- all_results |> filter(protein == !!protein)
  
  # Calculate statistics
  stats <- protein_data |>
    group_by(region_type) |>
    summarize(
      median = median(log2FC, na.rm = TRUE),
      mean = mean(log2FC, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  cat("\n", protein, "statistics:\n")
  print(stats)
  
  # Wilcoxon test
  anchor_vals <- protein_data |> filter(region_type == "Loop Anchor") |> pull(log2FC)
  between_vals <- protein_data |> filter(region_type == "Between Anchors") |> pull(log2FC)
  
  if (length(anchor_vals) > 0 && length(between_vals) > 0) {
    test_result <- wilcox.test(anchor_vals, between_vals)
    cat("Wilcoxon test p-value:", test_result$p.value, "\n")
  }
  
  # Plot
  p <- ggplot(protein_data, aes(x = region_type, y = log2FC, fill = region_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = c("Loop Anchor" = "#E64B35", 
                                 "Between Anchors" = "#4DBBD5")) +
    labs(
      title = paste(protein, "Binding: Loop Anchors vs Between Regions"),
      subtitle = paste("Gained loops (n =", length(gainedLoops), ")"),
      x = "",
      y = "log2FoldChange(Treated/Control)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 11)
    ) +
    coord_cartesian(ylim = c(-3, 3))
  
  print(p)
}

dev.off()

cat("\n\nAnalysis complete!\n")
cat("Results saved to: data/processed/cutntag/anchor_vs_between_analysis.csv\n")
cat("Plots saved to: plots/cutntag_anchor_vs_between_boxplots.pdf\n")

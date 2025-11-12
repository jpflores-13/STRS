## MA plots for Cut&Tag data with peak filtering

# Load libraries ----------------------------------------------------------

library(bamsignals)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(ggplot2)
library(dplyr)

# Data --------------------------------------------------------------------

## Read narrowPeak files and keep only standard chromosomes
cutntag_raw <- list.files("data/processed/cutntag/output/peaks/",
                          full.names = TRUE,
                          pattern = ".narrowPeak") |> 
  lapply(read_narrowpeaks) |> 
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

## Define our proteins of interest and conditions
target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag_raw) <- paste0(rep(target, each = 2), "_", condition)

## Filter peaks by score >= 200 (removes low-confidence peaks)
cutntag <- lapply(cutntag_raw, function(x) x[x$score >= 200])

cat("Peak counts after filtering (score >= 200):\n")
lapply(cutntag, length) |> print()

## Set up paths to BAM files
bam_files <- list.files("data/processed/cutntag/output/mergeAlign/",
                        full.names = TRUE,
                        pattern = ".bam$")
names(bam_files) <- paste0(rep(target, each = 2), "_", condition)

# Helper Functions --------------------------------------------------------

## Function to merge peaks for a single protein
merge_protein_peaks <- function(protein, peak_list) {
  # Find indices for control and treatment
  control_idx <- which(names(peak_list) == paste0(protein, "_control"))
  treatment_idx <- which(names(peak_list) == paste0(protein, "_sorbitol"))
  
  # Combine & merge peaks
  combined_peaks <- GenomicRanges::reduce(
    c(peak_list[[control_idx]], peak_list[[treatment_idx]])
  )
  
  return(combined_peaks)
}

## Function to quantify peaks using bamCount
get_peak_counts <- function(peaks, bam_file) {
  counts <- bamCount(bam_file, 
                     peaks,
                     paired.end = "midpoint")
  
  return(counts)
}

## Function to analyze count distributions
analyze_count_distribution <- function(peak_list, bam_files, target) {
  # Create an empty list to store count data for each protein
  count_data_list <- list()
  
  # Process each protein
  for (protein in target) {
    # Merge peaks from control and treatment conditions
    merged_peaks <- merge_protein_peaks(protein, peak_list)
    
    # Get the correct BAM file paths
    control_bam <- bam_files[names(bam_files) == paste0(protein, "_control")]
    treatment_bam <- bam_files[names(bam_files) == paste0(protein, "_sorbitol")]
    
    # Get counts
    control_counts <- get_peak_counts(merged_peaks, control_bam)
    treatment_counts <- get_peak_counts(merged_peaks, treatment_bam)
    
    # Calculate mean counts
    mean_counts <- (control_counts + treatment_counts) / 2
    
    # Store in data frame
    count_data_list[[protein]] <- data.frame(
      mean_count = mean_counts,
      protein = protein
    )
  }
  
  # Combine all data
  all_count_data <- do.call(rbind, count_data_list)
  
  # Create distribution plot
  count_dist_plot <- ggplot(all_count_data, aes(x = mean_count + 1)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    scale_x_log10(
      breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000),
      labels = c("1", "5", "10", "50", "100", "500", "1000", "5000")
    ) +
    geom_vline(xintercept = c(5, 10, 20), 
               linetype = "dashed", 
               color = "red", 
               alpha = 0.5) +
    facet_wrap(~protein, scales = "free_y") +
    labs(
      title = "Distribution of Mean Peak Counts",
      subtitle = "Filtered peaks (score ≥ 200). Red lines: potential thresholds (5, 10, 20)",
      x = "Mean Count (log scale)",
      y = "Number of Peaks"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )
  
  # Create summary statistics
  summary_stats <- all_count_data |> 
    dplyr::group_by(protein) |> 
    dplyr::summarise(
      total_peaks = dplyr::n(),
      peaks_above_5 = sum(mean_count >= 5),
      peaks_above_10 = sum(mean_count >= 10),
      peaks_above_20 = sum(mean_count >= 20),
      pct_above_5 = round(peaks_above_5 / total_peaks * 100, 1),
      pct_above_10 = round(peaks_above_10 / total_peaks * 100, 1),
      pct_above_20 = round(peaks_above_20 / total_peaks * 100, 1),
      .groups = "drop"
    )
  
  return(list(plot = count_dist_plot, stats = summary_stats))
}

## Function to create MA plot data for a single protein (updated with filtering)
create_ma_data <- function(protein, peak_list, bam_files, min_mean_count = 5) {
  # Merge peaks from control and treatment conditions
  merged_peaks <- merge_protein_peaks(protein, peak_list)
  
  # Get the correct BAM file paths using exact name matching
  control_bam <- bam_files[names(bam_files) == paste0(protein, "_control")]
  treatment_bam <- bam_files[names(bam_files) == paste0(protein, "_sorbitol")]
  
  # Get counts for control and treatment
  control_counts <- get_peak_counts(merged_peaks, control_bam)
  treatment_counts <- get_peak_counts(merged_peaks, treatment_bam)
  
  # Calculate mean counts for filtering
  mean_counts <- (control_counts + treatment_counts) / 2
  
  # Filter out low count peaks
  keep_peaks <- mean_counts >= min_mean_count
  
  # Apply filtering to all relevant vectors
  control_counts <- control_counts[keep_peaks]
  treatment_counts <- treatment_counts[keep_peaks]
  mean_counts <- mean_counts[keep_peaks]
  
  # Calculate log2 fold changes for remaining peaks
  log2fc <- log2((treatment_counts + 1) / (control_counts + 1))
  
  # Classify changes based on log2 fold change thresholds
  change_type <- case_when(
    log2fc >= 1 ~ "Increased",
    log2fc <= -1 ~ "Decreased",
    TRUE ~ "Unchanged"
  )
  
  # Create data frame for plotting
  plot_data <- data.frame(
    mean_signal = mean_counts,
    log2FC = log2fc,
    change_type = change_type,
    protein = protein
  )
  
  return(plot_data)
}

## Function to calculate summary statistics
get_summary_stats <- function(data) {
  # Convert to df
  data <- as.data.frame(data)
  
  # Calculate counts and percentages for each change type
  stats <- data |>
    dplyr::group_by(change_type) |>
    dplyr::summarise(
      n = dplyr::n(),
      percent = (dplyr::n() / nrow(data)) * 100,
      .groups = "drop"
    )
  
  # Ensure we have all change types, even if count is 0
  all_types <- data.frame(
    change_type = c("Increased", "Decreased", "Unchanged"),
    n = 0,
    percent = 0
  )
  
  # Combine with actual results, keeping all types
  stats <- dplyr::full_join(stats, all_types, by = "change_type") |>
    dplyr::mutate(
      n = coalesce(n.x, n.y),
      percent = coalesce(percent.x, percent.y)
    ) |>
    dplyr::select(change_type, n, percent)
  
  return(stats)
}

## Create individual MA plots (updated to include filtering info in subtitle)
create_ma_plot <- function(data, protein_name, min_mean_count) {
  # Calculate stats
  summary_stats <- get_summary_stats(data)
  
  # Subtitle text + stats
  up_stats <- summary_stats |> 
    dplyr::filter(change_type == "Increased")
  down_stats <- summary_stats |> 
    dplyr::filter(change_type == "Decreased")
  
  subtitle_text <- sprintf(
    "Filtered peaks (score ≥ 200, mean count ≥ %d)\n%d up (%.1f%%), %d down (%.1f%%)",
    min_mean_count,
    up_stats$n,
    up_stats$percent,
    down_stats$n,
    down_stats$percent
  )
  
  # Plot
  ggplot(data, aes(x = mean_signal, y = log2FC, color = change_type)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0,
               linetype = 2,
               color = "grey40") +
    geom_hline(yintercept = c(-1, 1),
               linetype = 3,
               color = "red",
               alpha = 0.5) +
    scale_color_manual(values = c(
      "Increased" = "#E64B35",
      "Decreased" = "#4DBBD5",
      "Unchanged" = "grey80"
    )) +
    scale_x_log10() +
    labs(
      title = paste(protein_name),
      subtitle = subtitle_text,
      x = "mean of normalized signal",
      y = "log2FoldChange(sorbitol/control)",
      color = "Change Type"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# Analysis ----------------------------------------------------------------

## First, analyze count distributions to help choose filtering threshold
count_analysis <- analyze_count_distribution(cutntag, bam_files, target)

## Save count distribution plot
pdf("plots/count_distribution_filtered.pdf", width = 10, height = 8)
print(count_analysis$plot)
dev.off()

## Print summary statistics to help choose threshold
cat("\nCount distribution summary (after score filtering):\n")
print(count_analysis$stats)

## Set filtering threshold based on count distribution analysis
min_mean_count <- 5  # Adjust this value based on your count distribution analysis

## Process all proteins and create MA plots
ma_plot_list <- list()

## Process each protein
for (prot in target) {
  ma_plot_list[[prot]] <- create_ma_data(prot, cutntag, bam_files, min_mean_count)
}

## Combine data frames
ma_plot_data <- do.call(rbind, ma_plot_list)

## Print summary of final filtered data
cat("\nSummary of final filtered data:\n")
final_summary <- ma_plot_data |>
  group_by(protein) |>
  summarize(
    total_peaks = n(),
    increased = sum(change_type == "Increased"),
    decreased = sum(change_type == "Decreased"),
    unchanged = sum(change_type == "Unchanged"),
    .groups = "drop"
  )
print(final_summary)

## Create individual plot objects for each protein
plot_ctcf <- create_ma_plot(subset(ma_plot_data, protein == "CTCF"), "CTCF", min_mean_count)
plot_h3k27ac <- create_ma_plot(subset(ma_plot_data, protein == "H3K27ac"), "H3K27ac", min_mean_count)
plot_rad21 <- create_ma_plot(subset(ma_plot_data, protein == "RAD21"), "RAD21", min_mean_count)
plot_yap1 <- create_ma_plot(subset(ma_plot_data, protein == "YAP1"), "YAP1", min_mean_count)

## Store in list
ma_plots <- list(
  CTCF = plot_ctcf,
  H3K27ac = plot_h3k27ac,
  RAD21 = plot_rad21,
  YAP1 = plot_yap1
)

## Save as PDF
pdf("plots/cutntag_MAplots_filtered.pdf", width = 8, height = 6)
for (plot in ma_plots) {
  print(plot)
}
dev.off()

## Save processed data
saveRDS(ma_plot_data, "data/processed/cutntag/ma_plot_data_filtered.rds")

cat("\nAnalysis complete!\n")
cat("Count distribution plot: plots/count_distribution_filtered.pdf\n")
cat("MA plots: plots/cutntag_MAplots_filtered.pdf\n")
cat("Processed data: data/processed/cutntag/ma_plot_data_filtered.rds\n")

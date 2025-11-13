## MA plots for Cut&Tag data with peak filtering
## Updated to match differential loop MA plot style

# Load libraries ----------------------------------------------------------

library(bamsignals)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(ggplot2)
library(dplyr)

# Data --------------------------------------------------------------------

## Read narrowPeak files and keep only standard chromosomes
cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE,
                      pattern = ".narrowPeak") |> 
  lapply(read_narrowpeaks) |> 
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

## Define our proteins of interest and conditions
target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

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

## Function to create MA plot data for a single protein (updated with filtering)
create_ma_data <- function(protein, peak_list, bam_files, min_mean_count = 5, log2fc_threshold = 1) {
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
  # Using terminology matching the differential loops analysis
  change_type <- case_when(
    log2fc > log2fc_threshold ~ "TRUE - upreg",
    log2fc < -log2fc_threshold ~ "TRUE - downreg",
    TRUE ~ "FALSE"
  )
  
  # Create data frame for plotting
  plot_data <- data.frame(
    mean_signal = mean_counts,
    log2FC = log2fc,
    change_type = change_type,
    protein = protein
  )
  
  # Arrange by change_type to ensure plotting order (FALSE plotted first)
  plot_data <- plot_data |> 
    arrange(change_type)
  
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
    change_type = c("TRUE - upreg", "TRUE - downreg", "FALSE"),
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

## Create individual MA plots (updated to match differential loop style)
create_ma_plot <- function(data, protein_name, min_mean_count, log2fc_threshold = 1, ylim_range = c(-4, 4)) {
  # Calculate stats
  summary_stats <- get_summary_stats(data)
  
  # Get counts for annotation
  up_count <- summary_stats |> 
    dplyr::filter(change_type == "TRUE - upreg") |> 
    pull(n)
  down_count <- summary_stats |> 
    dplyr::filter(change_type == "TRUE - downreg") |> 
    pull(n)
  
  # Determine annotation position based on data range
  max_x <- max(data$mean_signal, na.rm = TRUE)
  annotation_x <- max_x * 0.7  # Position at 70% of max x value
  
  # Create plot matching differential loop style
  plot_obj <- ggplot(data, aes(x = mean_signal, y = log2FC, color = change_type)) +
    geom_point(alpha = 1) +
    geom_hline(yintercept = 0,
               linetype = 2,
               color = "grey40") +
    scale_color_manual(values = c(
      "TRUE - upreg" = "#F8766D",
      "TRUE - downreg" = "#619CFF",
      "FALSE" = "grey80"
    )) +
    ylim(ylim_range) +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    labs(
      title = protein_name,
      y = paste0(protein_name, " log2FoldChange (sorbitol/control)"),
      x = "mean of normalized counts"
    ) +
    theme_classic() +
    theme(
      legend.position = "NONE",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    annotate(geom = "text",
             label = "Gained",
             x = annotation_x,
             y = ylim_range[2] * 0.875,  # 87.5% of upper limit
             color = "#F8766D") +
    annotate(geom = "text",
             label = paste0("n = ", up_count),
             x = annotation_x,
             y = ylim_range[2] * 0.775,  # 77.5% of upper limit
             color = "#F8766D") +
    annotate(geom = "text",
             label = "Lost",
             x = annotation_x,
             y = ylim_range[1] * 0.875,  # -87.5% of lower limit
             color = "#619CFF") +
    annotate(geom = "text",
             label = paste0("n = ", down_count),
             x = annotation_x,
             y = ylim_range[1] * 0.975,  # -97.5% of lower limit
             color = "#619CFF")
  
  return(plot_obj)
}

## Create density plot companion (matching differential loop style)
create_density_plot <- function(data, ylim_range = c(-4, 4)) {
  ggplot(data, aes(y = log2FC)) +
    geom_density(
      color = "#619CFF",
      fill = 4,
      alpha = 0.25
    ) +
    geom_hline(yintercept = 0,
               linetype = 2,
               color = "grey40") +
    ylim(ylim_range) +
    xlim(c(0, 1.5)) +
    theme_classic() +
    theme(
      legend.position = "NONE",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line.x = element_blank()
    )
}

# Analysis ----------------------------------------------------------------

## Set filtering threshold based on your needs
min_mean_count <- 5  # Adjust this value based on your count distribution needs
log2fc_threshold <- 1  # log2FC threshold for classification (same as differential loops)

## Process all proteins and create MA plots
ma_plot_list <- list()

## Process each protein
for (prot in target) {
  ma_plot_list[[prot]] <- create_ma_data(prot, cutntag, bam_files, 
                                         min_mean_count, log2fc_threshold)
}

## Combine data frames
ma_plot_data <- do.call(rbind, ma_plot_list)

## Print summary of final filtered data
cat("\nSummary of final filtered data:\n")
final_summary <- ma_plot_data |>
  group_by(protein) |>
  summarize(
    total_peaks = n(),
    gained = sum(change_type == "TRUE - upreg"),
    lost = sum(change_type == "TRUE - downreg"),
    unchanged = sum(change_type == "FALSE"),
    .groups = "drop"
  )
print(final_summary)

## Determine appropriate y-axis limits for each protein
## Calculate 95th percentile of absolute log2FC to set reasonable limits
ylim_values <- ma_plot_data |>
  group_by(protein) |>
  summarize(
    max_abs_log2fc = quantile(abs(log2FC), 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(ylim = ceiling(max_abs_log2fc))

## Create individual plot objects for each protein with appropriate ylims
ma_plots <- list()
density_plots <- list()

for (prot in target) {
  protein_data <- subset(ma_plot_data, protein == prot)
  protein_ylim <- ylim_values |> 
    filter(protein == prot) |> 
    pull(ylim)
  
  # Use at least 4 for ylim
  protein_ylim <- max(protein_ylim, 4)
  ylim_range <- c(-protein_ylim, protein_ylim)
  
  ma_plots[[prot]] <- create_ma_plot(protein_data, prot, 
                                     min_mean_count, log2fc_threshold, ylim_range)
  density_plots[[prot]] <- create_density_plot(protein_data, ylim_range)
}

## Save MA plots as PDF
pdf("plots/cutntag_MAplots.pdf", width = 8, height = 5)
for (plot in ma_plots) {
  print(plot)
}
dev.off()

## Save density plots as PDF
pdf("plots/cutntag_MAplots_density.pdf", width = 3, height = 5)
for (plot in density_plots) {
  print(plot)
}
dev.off()

cat("\nAnalysis complete!\n")
cat("MA plots: plots/cutntag_MAplots.pdf\n")
cat("Density plots: plots/cutntag_MAplots_density.pdf\n")
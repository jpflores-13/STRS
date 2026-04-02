# ##############################################################################
# filename:    cutntag_MAplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: MA plots for CUT&Tag data with peak filtering; matches
#              differential loop MA plot style; produces per-protein MA and
#              density plots
# ##############################################################################

# Libraries ----
library(bamsignals)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(ggplot2)
library(dplyr)

# Parameters ----
cutntag_peaks_dir   <- "data/processed/cutntag/output/peaks/"
bam_files_dir       <- "data/processed/cutntag/output/mergeAlign/"
output_ma_pdf       <- "plots/cutntag_MAplots.pdf"
output_density_pdf  <- "plots/cutntag_MAplots_density.pdf"
min_mean_count      <- 5
log2fc_threshold    <- 1
target              <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition           <- c("control", "sorbitol")

# Data import ----
cutntag <- list.files(cutntag_peaks_dir, full.names = TRUE,
                      pattern = ".narrowPeak") |>
  lapply(plyranges::read_narrowpeaks) |>
  lapply(keepStandardChromosomes, pruning.mode = "coarse")

names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

bam_files <- list.files(bam_files_dir, full.names = TRUE, pattern = ".bam$")
names(bam_files) <- paste0(rep(target, each = 2), "_", condition)

# Helper functions ----

merge_protein_peaks <- function(protein, peak_list) {
  control_idx   <- which(names(peak_list) == paste0(protein, "_control"))
  treatment_idx <- which(names(peak_list) == paste0(protein, "_sorbitol"))

  GenomicRanges::reduce(
    c(peak_list[[control_idx]], peak_list[[treatment_idx]])
  )
}

get_peak_counts <- function(peaks, bam_file) {
  bamCount(bam_file, peaks, paired.end = "midpoint")
}

create_ma_data <- function(protein, peak_list, bam_files,
                           min_mean_count = 5, log2fc_threshold = 1) {
  merged_peaks  <- merge_protein_peaks(protein, peak_list)

  control_bam   <- bam_files[names(bam_files) == paste0(protein, "_control")]
  treatment_bam <- bam_files[names(bam_files) == paste0(protein, "_sorbitol")]

  control_counts   <- get_peak_counts(merged_peaks, control_bam)
  treatment_counts <- get_peak_counts(merged_peaks, treatment_bam)
  mean_counts      <- (control_counts + treatment_counts) / 2

  keep_peaks       <- mean_counts >= min_mean_count
  control_counts   <- control_counts[keep_peaks]
  treatment_counts <- treatment_counts[keep_peaks]
  mean_counts      <- mean_counts[keep_peaks]

  log2fc      <- log2((treatment_counts + 1) / (control_counts + 1))
  change_type <- case_when(
    log2fc >  log2fc_threshold ~ "TRUE - upreg",
    log2fc < -log2fc_threshold ~ "TRUE - downreg",
    TRUE ~ "FALSE"
  )

  data.frame(mean_signal = mean_counts, log2FC = log2fc,
             change_type = change_type, protein = protein) |>
    arrange(change_type)
}

get_summary_stats <- function(data) {
  data <- as.data.frame(data)

  stats <- data |>
    dplyr::group_by(change_type) |>
    dplyr::summarise(n = dplyr::n(),
                     percent = (dplyr::n() / nrow(data)) * 100,
                     .groups = "drop")

  all_types <- data.frame(
    change_type = c("TRUE - upreg", "TRUE - downreg", "FALSE"),
    n = 0, percent = 0
  )

  dplyr::full_join(stats, all_types, by = "change_type") |>
    dplyr::mutate(n       = coalesce(n.x,       n.y),
                  percent = coalesce(percent.x, percent.y)) |>
    dplyr::select(change_type, n, percent)
}

create_ma_plot <- function(data, protein_name, min_mean_count,
                           log2fc_threshold = 1, ylim_range = c(-4, 4)) {
  summary_stats <- get_summary_stats(data)

  up_count   <- summary_stats |> dplyr::filter(change_type == "TRUE - upreg")   |> pull(n)
  down_count <- summary_stats |> dplyr::filter(change_type == "TRUE - downreg") |> pull(n)

  max_x        <- max(data$mean_signal, na.rm = TRUE)
  annotation_x <- max_x * 0.7

  ggplot(data, aes(x = mean_signal, y = log2FC, color = change_type)) +
    geom_point(alpha = 1) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    scale_color_manual(values = c("TRUE - upreg"   = "#F8766D",
                                  "TRUE - downreg" = "#619CFF",
                                  "FALSE"          = "grey80")) +
    ylim(ylim_range) +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    labs(title = protein_name,
         y = paste0(protein_name, " log2FoldChange (sorbitol/control)"),
         x = "mean of normalized counts") +
    theme_classic() +
    theme(legend.position = "NONE",
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    annotate(geom = "text", label = "Gained",
             x = annotation_x, y = ylim_range[2] * 0.875, color = "#F8766D") +
    annotate(geom = "text", label = paste0("n = ", up_count),
             x = annotation_x, y = ylim_range[2] * 0.775, color = "#F8766D") +
    annotate(geom = "text", label = "Lost",
             x = annotation_x, y = ylim_range[1] * 0.875, color = "#619CFF") +
    annotate(geom = "text", label = paste0("n = ", down_count),
             x = annotation_x, y = ylim_range[1] * 0.975, color = "#619CFF")
}

create_density_plot <- function(data, ylim_range = c(-4, 4)) {
  ggplot(data, aes(y = log2FC)) +
    geom_density(color = "#619CFF", fill = 4, alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    ylim(ylim_range) +
    xlim(c(0, 1.5)) +
    theme_classic() +
    theme(legend.position = "NONE",
          axis.text   = element_blank(),
          axis.title  = element_blank(),
          axis.ticks  = element_blank(),
          axis.line.x = element_blank())
}

# Analysis ----
ma_plot_list <- list()

for (prot in target) {
  ma_plot_list[[prot]] <- create_ma_data(prot, cutntag, bam_files,
                                         min_mean_count, log2fc_threshold)
}

ma_plot_data <- do.call(rbind, ma_plot_list)

cat("\nSummary of final filtered data:\n")
final_summary <- ma_plot_data |>
  group_by(protein) |>
  summarize(total_peaks = n(),
            gained      = sum(change_type == "TRUE - upreg"),
            lost        = sum(change_type == "TRUE - downreg"),
            unchanged   = sum(change_type == "FALSE"),
            .groups = "drop")
print(final_summary)

ylim_values <- ma_plot_data |>
  group_by(protein) |>
  summarize(max_abs_log2fc = quantile(abs(log2FC), 0.95, na.rm = TRUE),
            .groups = "drop") |>
  mutate(ylim = ceiling(max_abs_log2fc))

ma_plots      <- list()
density_plots <- list()

for (prot in target) {
  protein_data  <- subset(ma_plot_data, protein == prot)
  protein_ylim  <- ylim_values |> filter(protein == prot) |> pull(ylim)
  protein_ylim  <- max(protein_ylim, 4)
  ylim_range    <- c(-protein_ylim, protein_ylim)

  ma_plots[[prot]]      <- create_ma_plot(protein_data, prot,
                                          min_mean_count, log2fc_threshold, ylim_range)
  density_plots[[prot]] <- create_density_plot(protein_data, ylim_range)
}

# Save outputs ----
pdf(output_ma_pdf, width = 8, height = 5)
for (plot in ma_plots) print(plot)
dev.off()

pdf(output_density_pdf, width = 3, height = 5)
for (plot in density_plots) print(plot)
dev.off()

cat("\nAnalysis complete!\n")
cat("MA plots:", output_ma_pdf, "\n")
cat("Density plots:", output_density_pdf, "\n")

sessionInfo()

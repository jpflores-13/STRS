# ##############################################################################
# filename:    log2FC_vs_CPMbins_boxplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: log2FC (sorbitol/control) vs CPM-quantile bins for CTCF and
#              RAD21; median trend lines for both proteins on the same plot
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)
library(GenomeInfoDb)
library(Rsamtools)

# Parameters ----
ctcf_peaks_file  <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
rad21_peaks_file <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
ctcf_bam_ctl     <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_bam_trt     <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"
rad21_bam_ctl    <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_nodups_sorted.bam"
rad21_bam_trt    <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_nodups_sorted.bam"
output_pdf       <- "plots/log2FC_vs_CPMbins_combined_100bins.pdf"
n_bins           <- 100

# Helper functions ----
load_peaks <- function(file) {
  gr <- plyranges::read_narrowpeaks(file)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

bam_libsize <- function(bam) {
  sum(idxstatsBam(bam)$mapped, na.rm = TRUE)
}

analyze_by_cpm_bins <- function(peaks_gr,
                                bam_control,
                                bam_treat,
                                protein_name    = "CTCF",
                                n_bins          = 20,
                                pseudocount_cpm = 1e-6,
                                min_peaks_per_bin = 20,
                                paired_end_mode = "midpoint") {
  message("Peaks input: ", length(peaks_gr))

  ctl_counts <- bamCount(bam_control, peaks_gr, paired.end = paired_end_mode)
  trt_counts <- bamCount(bam_treat,   peaks_gr, paired.end = paired_end_mode)

  ctl_cpm <- ctl_counts / (bam_libsize(bam_control) / 1e6)
  trt_cpm <- trt_counts / (bam_libsize(bam_treat)   / 1e6)
  log2fc  <- log2((trt_cpm + pseudocount_cpm) / (ctl_cpm + pseudocount_cpm))

  ctl_cpm_pos <- pmax(ctl_cpm, pseudocount_cpm)
  qbreaks     <- unique(quantile(ctl_cpm,
                                 probs = seq(0, 1, length.out = n_bins + 1),
                                 na.rm = TRUE, type = 7))
  bin_factor  <- cut(ctl_cpm, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)

  df <- tibble(
    protein        = protein_name,
    bin            = bin_factor,
    ctl_cpm        = ctl_cpm,
    trt_cpm        = trt_cpm,
    log2fc         = log2fc,
    log10_ctl_cpm  = log10(ctl_cpm_pos)
  ) |> filter(!is.na(bin))

  bin_summ <- df |>
    group_by(bin) |>
    summarise(n                    = n(),
              median_log10_ctl_cpm = median(log10_ctl_cpm),
              median_log2fc        = median(log2fc, na.rm = TRUE),
              .groups              = "drop") |>
    arrange(as.integer(bin))

  keep_bins <- bin_summ$bin[bin_summ$n >= min_peaks_per_bin]
  df        <- df       |> filter(bin %in% keep_bins)
  bin_summ  <- bin_summ |> filter(bin %in% keep_bins)
  bin_summ$protein <- protein_name

  list(data = df, bins = bin_summ)
}

# Data import ----
ctcf_peaks  <- load_peaks(ctcf_peaks_file)
rad21_peaks <- load_peaks(rad21_peaks_file)

# Analysis ----
ctcf_res  <- analyze_by_cpm_bins(ctcf_peaks,  ctcf_bam_ctl,  ctcf_bam_trt,
                                 protein_name = "CTCF",  n_bins = n_bins)
rad21_res <- analyze_by_cpm_bins(rad21_peaks, rad21_bam_ctl, rad21_bam_trt,
                                 protein_name = "RAD21", n_bins = n_bins)

combined_bins <- bind_rows(ctcf_res$bins, rad21_res$bins) |>
  group_by(protein) |>
  mutate(percentile = (as.integer(bin) - 0.5) / max(as.integer(bin)) * 100) |>
  ungroup()

# Visualization ----
protein_colors <- c("CTCF" = "#008B8B", "RAD21" = "#FF6B6B")

p_combined <- ggplot(combined_bins,
                     aes(x = percentile, y = median_log2fc,
                         color = protein, group = protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.6) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  annotate("text", x = 98,  y = -1.2, label = "CTCF",
           color = protein_colors["CTCF"],  fontface = "bold", size = 4.5,
           hjust = 0, vjust = 0.5) +
  annotate("text", x = 100, y =  0.25, label = "RAD21",
           color = protein_colors["RAD21"], fontface = "bold", size = 4.5,
           hjust = 0.5, vjust = 0.5) +
  scale_color_manual(values = protein_colors) +
  scale_x_continuous(limits = c(58, 100), expand = c(0, 0)) +
  coord_cartesian(xlim = c(58, 100), ylim = c(-2.3, 0.5), clip = "off") +
  labs(x = "CPM percentile", y = "log2FC (sorbitol / control)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line       = element_line(color = "black"),
    plot.margin     = margin(5.5, 40, 5.5, 5.5, "pt")
  )

# Save outputs ----
pdf(output_pdf, width = 9, height = 6, useDingbats = FALSE)
print(p_combined)
dev.off()

sessionInfo()

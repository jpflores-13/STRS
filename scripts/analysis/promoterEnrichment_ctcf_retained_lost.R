# ##############################################################################
# filename:    promoterEnrichment_ctcf_retained_lost.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Density, boxplot, and violin plots comparing log2FC
#              (sorbitol/control) of CTCF peaks at promoters vs non-promoters;
#              Wilcoxon, KS test, and Cohen's d effect size
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)
library(ggsignif)
library(patchwork)

# Parameters ----
ctcf_peaks_file   <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
ctcf_control_bam  <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_sorbitol_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"
output_density    <- "plots/ctcf_promoter_density.pdf"
output_boxplot    <- "plots/ctcf_promoter_boxplot.pdf"
output_violin     <- "plots/ctcf_promoter_violin.pdf"
output_combined   <- "plots/ctcf_promoter_combined.pdf"
min_counts        <- 10
pseudocount       <- 1

# Data import ----
ctcf_control <- plyranges::read_narrowpeaks(ctcf_peaks_file) |>
  keepStandardChromosomes(pruning.mode = "coarse")

seqlevelsStyle(ctcf_control) <- "UCSC"

cat("CTCF control peaks:", length(ctcf_control), "\n")

txdb             <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter_regions <- promoters(genes(txdb), upstream = 2000, downstream = 2000)

cat("Promoter regions:", length(promoter_regions), "\n")

# Analysis ----
promoter_overlaps <- countOverlaps(ctcf_control, promoter_regions) > 0

mcols(ctcf_control)$overlaps_promoter  <- promoter_overlaps
mcols(ctcf_control)$promoter_category  <- ifelse(promoter_overlaps, "Promoter", "Non-promoter")

cat(sprintf("Peaks at promoters:     %d (%.1f%%)\n",
            sum(promoter_overlaps),  100 * mean(promoter_overlaps)))
cat(sprintf("Peaks at non-promoters: %d (%.1f%%)\n",
            sum(!promoter_overlaps), 100 * mean(!promoter_overlaps)))

control_counts  <- bamCount(ctcf_control_bam,  ctcf_control, paired.end = "midpoint")
sorbitol_counts <- bamCount(ctcf_sorbitol_bam, ctcf_control, paired.end = "midpoint")

mcols(ctcf_control)$control_counts  <- control_counts
mcols(ctcf_control)$sorbitol_counts <- sorbitol_counts
mcols(ctcf_control)$log2fc          <- log2((sorbitol_counts + pseudocount) /
                                            (control_counts  + pseudocount))

ctcf_filtered <- ctcf_control[control_counts >= min_counts]

cat(sprintf("Peaks after filtering (>= %d counts): %d of %d (%.1f%%)\n",
            min_counts, length(ctcf_filtered), length(ctcf_control),
            100 * length(ctcf_filtered) / length(ctcf_control)))

peak_data <- as.data.frame(mcols(ctcf_filtered))

summary_stats <- peak_data |>
  group_by(promoter_category) |>
  summarise(n_peaks       = n(),
            median_log2fc = median(log2fc),
            mean_log2fc   = mean(log2fc),
            sd_log2fc     = sd(log2fc),
            q25           = quantile(log2fc, 0.25),
            q75           = quantile(log2fc, 0.75),
            .groups       = "drop")

print(summary_stats)

## Statistical tests ----
promoter_log2fc    <- peak_data$log2fc[peak_data$promoter_category == "Promoter"]
nonpromoter_log2fc <- peak_data$log2fc[peak_data$promoter_category == "Non-promoter"]

wilcox_result <- wilcox.test(log2fc ~ promoter_category, data = peak_data,
                              alternative = "two.sided")
ks_result     <- ks.test(promoter_log2fc, nonpromoter_log2fc)

n_prom    <- length(promoter_log2fc)
n_nonprom <- length(nonpromoter_log2fc)
pooled_sd <- sqrt(((n_prom - 1)    * sd(promoter_log2fc)^2 +
                   (n_nonprom - 1) * sd(nonpromoter_log2fc)^2) /
                  (n_prom + n_nonprom - 2))
cohens_d  <- (mean(promoter_log2fc) - mean(nonpromoter_log2fc)) / pooled_sd

cat(sprintf("\nWilcoxon: W = %.0f, p = %.2e\n", wilcox_result$statistic, wilcox_result$p.value))
cat(sprintf("KS test:  D = %.4f, p = %.2e\n",   ks_result$statistic,    ks_result$p.value))
cat(sprintf("Cohen's d = %.4f\n", cohens_d))

get_significance_label <- function(p) {
  if (p < 0.0001) "****" else if (p < 0.001) "***" else if (p < 0.01) "**" else
    if (p < 0.05) "*" else "ns"
}

sig_label <- get_significance_label(wilcox_result$p.value)

# Visualization ----
fill_colors <- c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF")

p_density <- ggplot(peak_data, aes(x = log2fc, fill = promoter_category,
                                   color = promoter_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values  = fill_colors, name = "") +
  scale_color_manual(values = fill_colors, name = "") +
  labs(title    = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
       subtitle = sprintf("Wilcoxon p = %.2e (%s)", wilcox_result$p.value, sig_label),
       x = "log2 Fold Change (sorbitol/control)", y = "Density") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        panel.grid.minor = element_blank())

p_boxplot <- ggplot(peak_data, aes(x = promoter_category, y = log2fc,
                                   fill = promoter_category)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_signif(comparisons = list(c("Promoter", "Non-promoter")),
              map_signif_level = TRUE, test = "wilcox.test",
              textsize = 5, vjust = 0.5, tip_length = 0.02) +
  scale_fill_manual(values = fill_colors) +
  labs(title    = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
       subtitle = sprintf("Cohen's d = %.3f", cohens_d),
       x = "", y = "log2 Fold Change (sorbitol/control)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank())

p_violin <- ggplot(peak_data, aes(x = promoter_category, y = log2fc,
                                  fill = promoter_category)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0.5, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_signif(comparisons = list(c("Promoter", "Non-promoter")),
              map_signif_level = TRUE, test = "wilcox.test",
              textsize = 5, vjust = 0.5, tip_length = 0.02) +
  scale_fill_manual(values = fill_colors) +
  labs(title    = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
       subtitle = sprintf("KS test p = %.2e", ks_result$p.value),
       x = "", y = "log2 Fold Change (sorbitol/control)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none",
        panel.grid.minor = element_blank())

p_combined <- p_density / (p_boxplot | p_violin) +
  plot_annotation(title = "CTCF Peak Response to Sorbitol Treatment",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

# Save outputs ----
ggsave(output_density,  p_density,  width = 8,  height = 6)
ggsave(output_boxplot,  p_boxplot,  width = 7,  height = 6)
ggsave(output_violin,   p_violin,   width = 7,  height = 6)
ggsave(output_combined, p_combined, width = 12, height = 10)

sessionInfo()

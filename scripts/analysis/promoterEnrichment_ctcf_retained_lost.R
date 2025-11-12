## CTCF Peaks: Promoter Overlap vs Log2FC

# Load libraries ----------------------------------------------------------

library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)
library(ggsignif)

# Helper function to read narrowPeak files --------------------------------

read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name", 
                                    "score", "strand", "signalValue", 
                                    "pValue", "qValue", "peak"))
  
  GRanges(seqnames = peaks$seqnames,
          ranges = IRanges(start = peaks$start + 1, end = peaks$end),
          strand = "*",
          score = peaks$score,
          signalValue = peaks$signalValue,
          pValue = peaks$pValue,
          qValue = peaks$qValue,
          peak = peaks$peak)
}

# Load CTCF control peaks -------------------------------------------------

cat("Loading CTCF control peaks...\n")

ctcf_control <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")

## Harmonize chromosome naming
seqlevelsStyle(ctcf_control) <- "UCSC"

cat("CTCF control peaks:", length(ctcf_control), "\n")

# Load promoter annotations -----------------------------------------------

cat("\nLoading promoter annotations...\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)

## Define promoters as TSS ± 2kb (default: -2000/+200)
## Using symmetric window for clarity
promoter_regions <- promoters(genes, upstream = 2000, downstream = 2000)

cat("Number of promoter regions:", length(promoter_regions), "\n")

# Load BAM files for quantification ---------------------------------------

cat("\nBAM files for signal quantification:\n")

## Control BAM (0h)
ctcf_control_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"

## Sorbitol BAM (adjust filename to match your data)
ctcf_sorbitol_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"

cat("  Control:", ctcf_control_bam, "\n")
cat("  Sorbitol:", ctcf_sorbitol_bam, "\n")

# Classify peaks by promoter overlap --------------------------------------

cat("\nClassifying peaks by promoter overlap...\n")

## Find peaks that overlap promoters
promoter_overlaps <- countOverlaps(ctcf_control, promoter_regions) > 0

## Add annotation to peaks
mcols(ctcf_control)$overlaps_promoter <- promoter_overlaps
mcols(ctcf_control)$promoter_category <- ifelse(
  promoter_overlaps,
  "Promoter",
  "Non-promoter"
)

cat(sprintf("  Peaks at promoters: %d (%.1f%%)\n",
            sum(promoter_overlaps),
            100 * mean(promoter_overlaps)))
cat(sprintf("  Peaks at non-promoters: %d (%.1f%%)\n",
            sum(!promoter_overlaps),
            100 * mean(!promoter_overlaps)))

# Calculate log2FC for all peaks -----------------------------------------

cat("\nCalculating log2 fold changes...\n")

## Count reads in both conditions
control_counts <- bamCount(ctcf_control_bam, ctcf_control, 
                           paired.end = "midpoint")
sorbitol_counts <- bamCount(ctcf_sorbitol_bam, ctcf_control, 
                            paired.end = "midpoint")

## Add to peak object
mcols(ctcf_control)$control_counts <- control_counts
mcols(ctcf_control)$sorbitol_counts <- sorbitol_counts

## Calculate log2 fold change with pseudocount
pseudocount <- 1
mcols(ctcf_control)$log2fc <- log2(
  (sorbitol_counts + pseudocount) / (control_counts + pseudocount)
)

# Filter peaks by minimum signal -----------------------------------------

cat("\nFiltering peaks...\n")

## Require minimum signal in control to avoid noisy estimates
min_counts <- 10

ctcf_filtered <- ctcf_control[control_counts >= min_counts]

cat(sprintf("Peaks after filtering (>=%d counts): %d of %d (%.1f%%)\n",
            min_counts,
            length(ctcf_filtered),
            length(ctcf_control),
            100 * length(ctcf_filtered) / length(ctcf_control)))

# Prepare data for plotting -----------------------------------------------

## Convert to data frame
peak_data <- as.data.frame(mcols(ctcf_filtered))

## Summary statistics
summary_stats <- peak_data |>
  group_by(promoter_category) |>
  summarise(
    n_peaks = n(),
    median_log2fc = median(log2fc),
    mean_log2fc = mean(log2fc),
    sd_log2fc = sd(log2fc),
    q25 = quantile(log2fc, 0.25),
    q75 = quantile(log2fc, 0.75),
    .groups = "drop"
  )

cat("\nSummary statistics by promoter category:\n")
print(summary_stats)

# Statistical tests -------------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("STATISTICAL TESTS\n")
cat(rep("=", 60), "\n", sep = "")

## 1. Wilcoxon rank-sum test (Mann-Whitney U test)
## Non-parametric test for difference in distributions
wilcox_result <- wilcox.test(
  log2fc ~ promoter_category, 
  data = peak_data,
  alternative = "two.sided"
)

cat("\n1. Wilcoxon Rank-Sum Test (Mann-Whitney U):\n")
cat(sprintf("   - Test statistic (W): %.0f\n", wilcox_result$statistic))
cat(sprintf("   - P-value: %.2e\n", wilcox_result$p.value))
cat(sprintf("   - Interpretation: %s\n",
            ifelse(wilcox_result$p.value < 0.05,
                   "*** SIGNIFICANT difference in log2FC distributions ***",
                   "No significant difference in log2FC distributions")))

## 2. Kolmogorov-Smirnov test
## Tests if two samples come from the same distribution
promoter_log2fc <- peak_data$log2fc[peak_data$promoter_category == "Promoter"]
nonpromoter_log2fc <- peak_data$log2fc[peak_data$promoter_category == "Non-promoter"]

ks_result <- ks.test(promoter_log2fc, nonpromoter_log2fc)

cat("\n2. Kolmogorov-Smirnov Test:\n")
cat(sprintf("   - Test statistic (D): %.4f\n", ks_result$statistic))
cat(sprintf("   - P-value: %.2e\n", ks_result$p.value))
cat(sprintf("   - Interpretation: %s\n",
            ifelse(ks_result$p.value < 0.05,
                   "*** SIGNIFICANT difference in distributions ***",
                   "No significant difference in distributions")))

## 3. Effect size: Cohen's d
## Measures the magnitude of difference between groups
mean_promoter <- mean(promoter_log2fc)
mean_nonpromoter <- mean(nonpromoter_log2fc)
sd_promoter <- sd(promoter_log2fc)
sd_nonpromoter <- sd(nonpromoter_log2fc)
n_promoter <- length(promoter_log2fc)
n_nonpromoter <- length(nonpromoter_log2fc)

## Pooled standard deviation
pooled_sd <- sqrt(((n_promoter - 1) * sd_promoter^2 + 
                     (n_nonpromoter - 1) * sd_nonpromoter^2) / 
                    (n_promoter + n_nonpromoter - 2))

cohens_d <- (mean_promoter - mean_nonpromoter) / pooled_sd

cat("\n3. Effect Size (Cohen's d):\n")
cat(sprintf("   - Cohen's d: %.4f\n", cohens_d))
cat(sprintf("   - Interpretation: %s\n",
            ifelse(abs(cohens_d) < 0.2, "Negligible effect",
                   ifelse(abs(cohens_d) < 0.5, "Small effect",
                          ifelse(abs(cohens_d) < 0.8, "Medium effect", "Large effect")))))

cat("\n", rep("=", 60), "\n", sep = "")

# Helper function to convert p-value to significance stars ----------------

get_significance_label <- function(p_value) {
  if (p_value < 0.0001) {
    return("****")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

sig_label <- get_significance_label(wilcox_result$p.value)

# Create plots ------------------------------------------------------------

## 1. Density plot
p_density <- ggplot(peak_data, 
                    aes(x = log2fc, 
                        fill = promoter_category,
                        color = promoter_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
    name = ""
  ) +
  scale_color_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
    name = ""
  ) +
  labs(
    title = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
    subtitle = sprintf("Wilcoxon p = %.2e (%s)", wilcox_result$p.value, sig_label),
    x = "log2 Fold Change (sorbitol/control)",
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

## 2. Boxplot with significance brackets
p_boxplot <- ggplot(peak_data, 
                    aes(x = promoter_category, 
                        y = log2fc,
                        fill = promoter_category)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_signif(
    comparisons = list(c("Promoter", "Non-promoter")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 5,
    vjust = 0.5,
    tip_length = 0.02
  ) +
  scale_fill_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF")
  ) +
  labs(
    title = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
    subtitle = sprintf("Cohen's d = %.3f", cohens_d),
    x = "",
    y = "log2 Fold Change (sorbitol/control)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

## 3. Violin plot with significance brackets
p_violin <- ggplot(peak_data, 
                   aes(x = promoter_category, 
                       y = log2fc,
                       fill = promoter_category)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0.5,
               fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_signif(
    comparisons = list(c("Promoter", "Non-promoter")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 5,
    vjust = 0.5,
    tip_length = 0.02
  ) +
  scale_fill_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF")
  ) +
  labs(
    title = "CTCF Response to Sorbitol: Promoter vs Non-promoter Peaks",
    subtitle = sprintf("KS test p = %.2e", ks_result$p.value),
    x = "",
    y = "log2 Fold Change (sorbitol/control)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Save plots --------------------------------------------------------------

dir.create("plots", showWarnings = FALSE)

## Save individual plots
ggsave(
  filename = "plots/ctcf_promoter_density.pdf",
  plot = p_density,
  width = 8,
  height = 6
)

ggsave(
  filename = "plots/ctcf_promoter_boxplot.pdf",
  plot = p_boxplot,
  width = 7,
  height = 6
)

ggsave(
  filename = "plots/ctcf_promoter_violin.pdf",
  plot = p_violin,
  width = 7,
  height = 6
)

## Save combined plot
library(patchwork)
p_combined <- p_density / (p_boxplot | p_violin) +
  plot_annotation(
    title = "CTCF Peak Response to Sorbitol Treatment",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  filename = "plots/ctcf_promoter_combined.pdf",
  plot = p_combined,
  width = 12,
  height = 10
)

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE!\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nPlots saved to:\n")
cat("  - plots/ctcf_promoter_density.pdf\n")
cat("  - plots/ctcf_promoter_boxplot.pdf\n")
cat("  - plots/ctcf_promoter_violin.pdf\n")
cat("  - plots/ctcf_promoter_combined.pdf\n")
cat("\nThis analysis:\n")
cat("  1. Classifies ALL control peaks by promoter overlap\n")
cat("  2. Calculates log2FC for each peak (no DESeq2 classification)\n")
cat("  3. Compares log2FC distributions between promoter and non-promoter peaks\n")
cat("  4. Tests statistical significance with multiple approaches\n")
cat("  5. Adds significance stars to visualizations\n")
cat("  6. Completely avoids power bias from DESeq2\n")

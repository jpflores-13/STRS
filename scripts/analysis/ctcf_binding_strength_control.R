# ##############################################################################
# filename:    ctcf_binding_strength_control.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Density plots of CTCF and RAD21 occupancy (control BAM signal)
#              at retained vs lost control peaks; Wilcoxon tests for each
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)

# Parameters ----
ctcf_narrowpeak      <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
rad21_narrowpeak     <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
ctcf_deseq_rds       <- "data/processed/cutntag/deseq2/diff_CTCF_counts.rds"
rad21_deseq_rds      <- "data/processed/cutntag/deseq2/diff_RAD21_counts.rds"
ctcf_control_bam     <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
rad21_control_bam    <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_nodups_sorted.bam"
output_pdf           <- "plots/control_peaks_density_plots.pdf"
count_threshold      <- 10

# Data import ----

## Helper to read narrowPeak files
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

cat("Loading CONTROL peak calls...\n")

ctcf_control  <- read_narrowpeaks(ctcf_narrowpeak)
ctcf_control  <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

rad21_control <- read_narrowpeaks(rad21_narrowpeak)
rad21_control <- keepStandardChromosomes(rad21_control, pruning.mode = "coarse")
cat("RAD21 control peaks:", length(rad21_control), "\n")

cat("\nLoading DESeq2 differential analysis results...\n")
ctcf_deseq  <- readRDS(ctcf_deseq_rds)
rad21_deseq <- readRDS(rad21_deseq_rds)

# Analysis ----

## Classify and quantify peaks ----
analyze_peak_retention <- function(control_peaks, deseq_results, control_bam,
                                   protein_name) {

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing", protein_name, "\n")
  cat(rep("=", 60), "\n", sep = "")

  overlaps <- findOverlaps(control_peaks, deseq_results)
  cat("Control peaks matched to DESeq2 results:", length(overlaps), "\n")

  matched_control_peaks <- control_peaks[queryHits(overlaps)]

  mcols(matched_control_peaks)$log2FoldChange <-
    mcols(deseq_results)$log2FoldChange[subjectHits(overlaps)]
  mcols(matched_control_peaks)$padj <-
    mcols(deseq_results)$padj[subjectHits(overlaps)]

  cat("Calculating signal at control peaks...\n")
  control_signal <- bamCount(control_bam, matched_control_peaks,
                             paired.end = "midpoint")
  mcols(matched_control_peaks)$control_signal <- control_signal

  mcols(matched_control_peaks)$peak_status <- case_when(
    !is.na(matched_control_peaks$padj) &
      matched_control_peaks$padj < 0.1 &
      matched_control_peaks$log2FoldChange > 0 &
      control_signal > count_threshold ~ "Gained",
    !is.na(matched_control_peaks$padj) &
      matched_control_peaks$padj < 0.1 &
      matched_control_peaks$log2FoldChange < 0 &
      control_signal > count_threshold ~ "Lost",
    control_signal > count_threshold ~ "Static",
    TRUE ~ "Filtered"
  )

  mcols(matched_control_peaks)$retention_category <- ifelse(
    matched_control_peaks$peak_status %in% c("Gained", "Static"),
    "Retained",
    ifelse(matched_control_peaks$peak_status == "Lost", "Lost", "Filtered")
  )

  cat("\nPeak classification:\n")
  cat("  Gained:", sum(matched_control_peaks$peak_status == "Gained"), "\n")
  cat("  Lost:", sum(matched_control_peaks$peak_status == "Lost"), "\n")
  cat("  Static:", sum(matched_control_peaks$peak_status == "Static"), "\n")
  cat("  Retained (Gained+Static):",
      sum(matched_control_peaks$retention_category == "Retained"), "\n")
  cat("  Filtered (low signal):",
      sum(matched_control_peaks$peak_status == "Filtered"), "\n")

  return(matched_control_peaks)
}

ctcf_classified  <- analyze_peak_retention(ctcf_control,  ctcf_deseq,  ctcf_control_bam,  "CTCF")
rad21_classified <- analyze_peak_retention(rad21_control, rad21_deseq, rad21_control_bam, "RAD21")

## Prepare data for plotting ----
ctcf_filtered  <- ctcf_classified[ctcf_classified$retention_category  %in% c("Retained", "Lost")]
rad21_filtered <- rad21_classified[rad21_classified$retention_category %in% c("Retained", "Lost")]

ctcf_df  <- as.data.frame(mcols(ctcf_filtered))  |> mutate(log10_signal = log10(control_signal + 1))
rad21_df <- as.data.frame(mcols(rad21_filtered)) |> mutate(log10_signal = log10(control_signal + 1))

## Statistics ----
ctcf_wilcox  <- wilcox.test(log10_signal ~ retention_category, data = ctcf_df)
rad21_wilcox <- wilcox.test(log10_signal ~ retention_category, data = rad21_df)

cat("\nCTCF (control peaks only):\n")
cat(sprintf("  Retained: n=%d, median signal=%.2f\n",
            sum(ctcf_df$retention_category == "Retained"),
            median(ctcf_df$log10_signal[ctcf_df$retention_category == "Retained"])))
cat(sprintf("  Lost: n=%d, median signal=%.2f\n",
            sum(ctcf_df$retention_category == "Lost"),
            median(ctcf_df$log10_signal[ctcf_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e\n", ctcf_wilcox$p.value))

cat("\nRAD21 (control peaks only):\n")
cat(sprintf("  Retained: n=%d, median signal=%.2f\n",
            sum(rad21_df$retention_category == "Retained"),
            median(rad21_df$log10_signal[rad21_df$retention_category == "Retained"])))
cat(sprintf("  Lost: n=%d, median signal=%.2f\n",
            sum(rad21_df$retention_category == "Lost"),
            median(rad21_df$log10_signal[rad21_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e\n", rad21_wilcox$p.value))

# Visualization ----

p_ctcf <- ggplot(ctcf_df,
                 aes(x = log10_signal, fill = retention_category,
                     color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values  = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "Peak Status") +
  scale_color_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                     name = "Peak Status") +
  labs(title = "CTCF Occupancy at Control Peaks: Retained vs Lost",
       subtitle = sprintf("Wilcoxon p = %.2e", ctcf_wilcox$p.value),
       x = "log10(Control Signal + 1)", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

p_rad21 <- ggplot(rad21_df,
                  aes(x = log10_signal, fill = retention_category,
                      color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values  = c("Retained" = "#41B6C4", "Lost" = "#E34A33"),
                    name = "Peak Status") +
  scale_color_manual(values = c("Retained" = "#41B6C4", "Lost" = "#E34A33"),
                     name = "Peak Status") +
  labs(title = "RAD21 Occupancy at Control Peaks: Retained vs Lost",
       subtitle = sprintf("Wilcoxon p = %.2e", rad21_wilcox$p.value),
       x = "log10(Control Signal + 1)", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

# Save outputs ----
dir.create("plots", showWarnings = FALSE)

pdf(output_pdf, width = 8, height = 5)
print(p_ctcf)
print(p_rad21)
dev.off()

cat("\nDensity plots saved to:", output_pdf, "\n")

sessionInfo()

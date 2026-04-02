# ##############################################################################
# filename:    ctcf_retention_by_rad21_fc.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Density, violin, and boxplot comparisons of RAD21 log2FC
#              (sorbitol/control) at retained vs lost CTCF peaks; Wilcoxon
#              test for significance
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(ggsignif)

# Parameters ----
ctcf_narrowpeak <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
ctcf_deseq_rds  <- "data/processed/cutntag/deseq2/diff_CTCF_counts.rds"
rad21_deseq_rds <- "data/processed/cutntag/deseq2/diff_RAD21_counts.rds"
plots_dir       <- "plots"

# Helper functions ----

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

# Data import ----
cat("Loading CONTROL peak calls...\n")

ctcf_control <- read_narrowpeaks(ctcf_narrowpeak)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

cat("\nLoading DESeq2 differential analysis results...\n")
ctcf_deseq  <- readRDS(ctcf_deseq_rds)
rad21_deseq <- readRDS(rad21_deseq_rds)

# Analysis ----

## Classify CTCF peaks as Retained vs Lost ----
ctcf_overlaps <- findOverlaps(ctcf_control, ctcf_deseq)
cat("CTCF control peaks matched to DESeq2 results:", length(ctcf_overlaps), "\n")

ctcf_matched <- ctcf_control[queryHits(ctcf_overlaps)]

mcols(ctcf_matched)$log2FoldChange_CTCF <-
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(ctcf_overlaps)]
mcols(ctcf_matched)$padj_CTCF <-
  mcols(ctcf_deseq)$padj[subjectHits(ctcf_overlaps)]

mcols(ctcf_matched)$peak_status <- case_when(
  !is.na(ctcf_matched$padj_CTCF) &
    ctcf_matched$padj_CTCF < 0.1 &
    ctcf_matched$log2FoldChange_CTCF > 0 ~ "Gained",
  !is.na(ctcf_matched$padj_CTCF) &
    ctcf_matched$padj_CTCF < 0.1 &
    ctcf_matched$log2FoldChange_CTCF < 0 ~ "Lost",
  !is.na(ctcf_matched$padj_CTCF) ~ "Static",
  TRUE ~ "Filtered"
)

mcols(ctcf_matched)$retention_category <- ifelse(
  ctcf_matched$peak_status %in% c("Gained", "Static"),
  "Retained",
  ifelse(ctcf_matched$peak_status == "Lost", "Lost", "Filtered")
)

cat("\nCTCF peak classification:\n")
cat("  Gained:", sum(ctcf_matched$peak_status == "Gained"), "\n")
cat("  Lost:", sum(ctcf_matched$peak_status == "Lost"), "\n")
cat("  Static:", sum(ctcf_matched$peak_status == "Static"), "\n")
cat("  Retained (Gained+Static):",
    sum(ctcf_matched$retention_category == "Retained"), "\n")
cat("  Filtered (no DESeq2 data):",
    sum(ctcf_matched$peak_status == "Filtered"), "\n")

## Match CTCF peaks to RAD21 log2FC ----
ctcf_rad21_overlaps <- findOverlaps(ctcf_matched, rad21_deseq)
cat("CTCF peaks overlapping with RAD21 DESeq2 results:",
    length(unique(queryHits(ctcf_rad21_overlaps))), "\n")

rad21_fc_by_ctcf <- mcols(rad21_deseq)$log2FoldChange[subjectHits(ctcf_rad21_overlaps)] |>
  split(queryHits(ctcf_rad21_overlaps)) |>
  sapply(mean, na.rm = TRUE)

mcols(ctcf_matched)$log2FC_RAD21 <- NA_real_
mcols(ctcf_matched)$log2FC_RAD21[as.integer(names(rad21_fc_by_ctcf))] <-
  rad21_fc_by_ctcf

cat("CTCF peaks with RAD21 log2FC data:",
    sum(!is.na(ctcf_matched$log2FC_RAD21)), "\n")

## Prepare plot data ----
ctcf_filtered <- ctcf_matched[
  ctcf_matched$retention_category %in% c("Retained", "Lost") &
    !is.na(ctcf_matched$log2FC_RAD21)
]

cat("\nFinal dataset:\n")
cat("  Retained CTCF peaks with RAD21 data:",
    sum(ctcf_filtered$retention_category == "Retained"), "\n")
cat("  Lost CTCF peaks with RAD21 data:",
    sum(ctcf_filtered$retention_category == "Lost"), "\n")

plot_df <- as.data.frame(mcols(ctcf_filtered)) |>
  select(retention_category, log2FC_RAD21)

## Statistics ----
wilcox_result <- wilcox.test(log2FC_RAD21 ~ retention_category, data = plot_df)

cat("\nRAD21 log2FC at CTCF peaks:\n")
cat(sprintf("  Retained: n=%d, median=%.3f\n",
            sum(plot_df$retention_category == "Retained"),
            median(plot_df$log2FC_RAD21[plot_df$retention_category == "Retained"])))
cat(sprintf("  Lost: n=%d, median=%.3f\n",
            sum(plot_df$retention_category == "Lost"),
            median(plot_df$log2FC_RAD21[plot_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e\n", wilcox_result$p.value))

# Visualization ----

## Density plot ----
p_density <- ggplot(plot_df,
                    aes(x = log2FC_RAD21,
                        fill  = retention_category,
                        color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values  = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  scale_color_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                     name = "CTCF Peak Status") +
  labs(title    = "RAD21 Changes at Retained vs Lost CTCF Peaks",
       subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
       x = "RAD21 log2FC (sorbitol/control)", y = "Density") +
  theme_bw() +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10),
        legend.position = "top")

## Violin plot ----
p_violin <- ggplot(plot_df,
                   aes(x = retention_category,
                       y = log2FC_RAD21,
                       fill = retention_category)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_signif(comparisons = list(c("Lost", "Retained")),
              map_signif_level = TRUE, test = "wilcox.test",
              textsize = 3.5, vjust = 0.5) +
  scale_fill_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  labs(title    = "RAD21 Changes at Retained vs Lost CTCF Peaks",
       subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
       x = "CTCF Peak Status", y = "RAD21 log2FC (sorbitol/control)") +
  theme_bw() +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10),
        legend.position = "none")

## Boxplot ----
p_box <- ggplot(plot_df,
                aes(x = retention_category,
                    y = log2FC_RAD21,
                    fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_signif(comparisons = list(c("Lost", "Retained")),
              map_signif_level = TRUE, test = "wilcox.test",
              textsize = 3.5, vjust = 0.5) +
  scale_fill_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  labs(title    = "RAD21 Changes at Retained vs Lost CTCF Peaks",
       subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
       x = "CTCF Peak Status", y = "RAD21 log2FC (sorbitol/control)") +
  theme_bw() +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10),
        legend.position = "none") +
  scale_y_continuous(limits = c(-2.5, 1))

# Save outputs ----
dir.create(plots_dir, showWarnings = FALSE)

pdf(file.path(plots_dir, "ctcf_retention_by_rad21_fc_density.pdf"), width = 8, height = 5)
print(p_density)
dev.off()

pdf(file.path(plots_dir, "ctcf_retention_by_rad21_fc_violin.pdf"), width = 7, height = 6)
print(p_violin)
dev.off()

pdf(file.path(plots_dir, "ctcf_retention_by_rad21_fc_box.pdf"), width = 7, height = 6)
print(p_box)
dev.off()

sessionInfo()

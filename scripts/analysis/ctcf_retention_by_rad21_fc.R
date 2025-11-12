## Density Plots: RAD21 log2FC at Retained vs Lost CTCF Peaks
## Shows how RAD21 changes relate to CTCF peak retention

# Load libraries ----------------------------------------------------------

library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
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

# Load control peak calls -------------------------------------------------

cat("Loading CONTROL peak calls...\n")

## CTCF control peaks
ctcf_control <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

# Load DESeq2 results -----------------------------------------------------

cat("\nLoading DESeq2 differential analysis results...\n")

ctcf_deseq <- readRDS("data/processed/cutntag/deseq2/diff_CTCF_counts.rds")
rad21_deseq <- readRDS("data/processed/cutntag/deseq2/diff_RAD21_counts.rds")

# Classify CTCF peaks as Retained vs Lost --------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Classifying CTCF peaks\n")
cat(rep("=", 60), "\n", sep = "")

## Find overlaps between CTCF control peaks and CTCF DESeq2 results
ctcf_overlaps <- findOverlaps(ctcf_control, ctcf_deseq)

cat("CTCF control peaks matched to DESeq2 results:", length(ctcf_overlaps), "\n")

## Extract matched control peaks
ctcf_matched <- ctcf_control[queryHits(ctcf_overlaps)]

## Add DESeq2 classification
mcols(ctcf_matched)$log2FoldChange_CTCF <- 
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(ctcf_overlaps)]
mcols(ctcf_matched)$padj_CTCF <- 
  mcols(ctcf_deseq)$padj[subjectHits(ctcf_overlaps)]

## Classify peaks as Retained vs Lost using DESeq2 results
## Retained = Gained (padj<0.1, log2FC>0) OR Static (padj>=0.1)
## Lost = Lost (padj<0.1, log2FC<0)

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

## Retention category
mcols(ctcf_matched)$retention_category <- ifelse(
  ctcf_matched$peak_status %in% c("Gained", "Static"),
  "Retained",
  ifelse(ctcf_matched$peak_status == "Lost", "Lost", "Filtered")
)

## Print summary
cat("\nCTCF peak classification:\n")
cat("  Gained:", sum(ctcf_matched$peak_status == "Gained"), "\n")
cat("  Lost:", sum(ctcf_matched$peak_status == "Lost"), "\n")
cat("  Static:", sum(ctcf_matched$peak_status == "Static"), "\n")
cat("  Retained (Gained+Static):", 
    sum(ctcf_matched$retention_category == "Retained"), "\n")
cat("  Filtered (no DESeq2 data):", 
    sum(ctcf_matched$peak_status == "Filtered"), "\n")

# Match CTCF peaks to RAD21 peaks -----------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Matching CTCF peaks to RAD21 log2FC\n")
cat(rep("=", 60), "\n", sep = "")

## Find overlaps between CTCF peaks and RAD21 DESeq2 results
## This will transfer RAD21 log2FC values to CTCF peaks
ctcf_rad21_overlaps <- findOverlaps(ctcf_matched, rad21_deseq)

cat("CTCF peaks overlapping with RAD21 DESeq2 results:", 
    length(unique(queryHits(ctcf_rad21_overlaps))), "\n")

## For CTCF peaks that overlap multiple RAD21 regions, take the mean log2FC
## This handles cases where a CTCF peak spans multiple RAD21 peaks
rad21_fc_by_ctcf <- mcols(rad21_deseq)$log2FoldChange[subjectHits(ctcf_rad21_overlaps)] |>
  split(queryHits(ctcf_rad21_overlaps)) |>
  sapply(mean, na.rm = TRUE)

## Initialize RAD21 log2FC column with NA
mcols(ctcf_matched)$log2FC_RAD21 <- NA_real_

## Assign mean RAD21 log2FC to corresponding CTCF peaks
mcols(ctcf_matched)$log2FC_RAD21[as.integer(names(rad21_fc_by_ctcf))] <- 
  rad21_fc_by_ctcf

cat("CTCF peaks with RAD21 log2FC data:", 
    sum(!is.na(ctcf_matched$log2FC_RAD21)), "\n")

# Prepare data for plotting -----------------------------------------------

## Filter to only Retained vs Lost CTCF peaks that have RAD21 data
ctcf_filtered <- ctcf_matched[
  ctcf_matched$retention_category %in% c("Retained", "Lost") &
    !is.na(ctcf_matched$log2FC_RAD21)
]

cat("\nFinal dataset for plotting:\n")
cat("  Retained CTCF peaks with RAD21 data:", 
    sum(ctcf_filtered$retention_category == "Retained"), "\n")
cat("  Lost CTCF peaks with RAD21 data:", 
    sum(ctcf_filtered$retention_category == "Lost"), "\n")

## Create data frame
plot_df <- as.data.frame(mcols(ctcf_filtered)) |>
  select(retention_category, log2FC_RAD21)

# Statistics --------------------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("STATISTICS\n")
cat(rep("=", 60), "\n", sep = "")

## Statistical test
wilcox_result <- wilcox.test(log2FC_RAD21 ~ retention_category, data = plot_df)

cat("\nRAD21 log2FC at CTCF peaks:\n")
cat(sprintf("  Retained CTCF: n=%d, median RAD21 log2FC=%.3f\n",
            sum(plot_df$retention_category == "Retained"),
            median(plot_df$log2FC_RAD21[plot_df$retention_category == "Retained"])))
cat(sprintf("  Lost CTCF: n=%d, median RAD21 log2FC=%.3f\n",
            sum(plot_df$retention_category == "Lost"),
            median(plot_df$log2FC_RAD21[plot_df$retention_category == "Lost"])))
cat(sprintf("  Median difference: %.3f\n",
            median(plot_df$log2FC_RAD21[plot_df$retention_category == "Retained"]) -
              median(plot_df$log2FC_RAD21[plot_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e\n", wilcox_result$p.value))

# Create density plot -----------------------------------------------------

p_density <- ggplot(plot_df, 
                    aes(x = log2FC_RAD21, 
                        fill = retention_category,
                        color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  scale_color_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "RAD21 Changes at Retained vs Lost CTCF Peaks",
    subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
    x = "RAD21 log2FC (sorbitol/control)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "top"
  )

# Create violin plot -----------------------------------------------------

p_violin <- ggplot(plot_df, 
                   aes(x = retention_category, 
                       y = log2FC_RAD21,
                       fill = retention_category)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_signif(
    comparisons = list(c("Lost", "Retained")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 3.5,
    vjust = 0.5
  ) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "RAD21 Changes at Retained vs Lost CTCF Peaks",
    subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
    x = "CTCF Peak Status",
    y = "RAD21 log2FC (sorbitol/control)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none"
  )

# Create boxplot ----------------------------------------------------------

p_box <- ggplot(plot_df, 
                aes(x = retention_category, 
                    y = log2FC_RAD21,
                    fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5,
               outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_signif(
    comparisons = list(c("Lost", "Retained")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 3.5,
    vjust = 0.5
  ) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "RAD21 Changes at Retained vs Lost CTCF Peaks",
    subtitle = sprintf("Wilcoxon p = %.2e", wilcox_result$p.value),
    x = "CTCF Peak Status",
    y = "RAD21 log2FC (sorbitol/control)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-2.5, 1))

# Save plots --------------------------------------------------------------

dir.create("plots", showWarnings = FALSE)

pdf("plots/ctcf_retention_by_rad21_fc_density.pdf", width = 8, height = 5)
print(p_density)
dev.off()

pdf("plots/ctcf_retention_by_rad21_fc_violin.pdf", width = 7, height = 6)
print(p_violin)
dev.off()

pdf("plots/ctcf_retention_by_rad21_fc_box.pdf", width = 7, height = 6)
print(p_box)
dev.off()

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE!\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nPlots saved to:\n")
cat("  - plots/ctcf_retention_by_rad21_fc_density.pdf\n")
cat("  - plots/ctcf_retention_by_rad21_fc_violin.pdf\n")
cat("  - plots/ctcf_retention_by_rad21_fc_box.pdf\n")
cat("\nThis analysis shows:\n")
cat("  How RAD21 changes (log2FC) differ between retained vs lost CTCF peaks\n")
cat("  Density plot: Distribution comparison with p-value\n")
cat("  Violin plot: Distribution shape with significance bracket\n")
cat("  Box plot: Median, quartiles, and significance bracket\n")

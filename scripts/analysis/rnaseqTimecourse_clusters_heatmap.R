# ##############################################################################
# filename:    rnaseqTimecourse_clusters_heatmap.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: ComplexHeatmap of significantly time-varying genes from the
#              RNA-seq LRT timecourse; variance-stabilized LFCs, replicate
#              averaging, hierarchical row clustering
# ##############################################################################

# Libraries ----
library(DESeq2)
library(DEGreport)
library(tibble)
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)
library(circlize)

# Parameters ----
dds_rds    <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
output_pdf <- "plots/rnaseqTimecourse-heatmap.pdf"
lfc_cutoff <- 2
padj_cutoff <- 0.05

# Data import ----
source("scripts/utils/make_norm_matrix.R")
dds <- readRDS(dds_rds)

# Analysis ----
coef <- resultsNames(dds) |> stringr::str_subset("Time")

message("Performing LFC shrinkage...")
resShrunk <- lapply(coef, function(c) {
  df <- lfcShrink(dds, coef = c, type = "apeglm")
  df$timepoint <- c
  df
})
shrunkenLFCs <- do.call(rbind, resShrunk)

message("Identifying significantly changing genes...")
sigLRT_genes <- shrunkenLFCs |>
  data.frame() |>
  mutate(timepoint = case_when(
    timepoint == "Time_1h_vs_0h"  ~ "1h",
    timepoint == "Time_3h_vs_0h"  ~ "3h",
    timepoint == "Time_6h_vs_0h"  ~ "6h",
    timepoint == "Time_9h_vs_0h"  ~ "9h",
    timepoint == "Time_12h_vs_0h" ~ "12h",
    timepoint == "Time_24h_vs_0h" ~ "24h"
  ),
  regulation = case_when(
    padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "upregulated",
    padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "downregulated"
  )) |>
  filter(regulation %in% c("upregulated", "downregulated")) |>
  rownames_to_column(var = "gene") |>
  pull(gene)

message("Found ", length(sigLRT_genes), " significantly changing genes")

count_matrix_norm <- make_norm_matrix(dds)
count_matrix_norm <- count_matrix_norm[rownames(count_matrix_norm) %in% sigLRT_genes, ]

message("Averaging biological replicates...")
n_tp       <- ncol(count_matrix_norm) / 2
combinedMat <- matrix(NA, nrow = nrow(count_matrix_norm), ncol = n_tp)
rownames(combinedMat) <- rownames(count_matrix_norm)

for (i in seq_len(n_tp)) {
  combinedMat[, i] <- (count_matrix_norm[, 2 * i - 1] + count_matrix_norm[, 2 * i]) / 2
}

colnames(combinedMat) <- unique(as.vector(colData(dds)$Time))

# Visualization ----
hm <- Heatmap(
  combinedMat,
  col = colorRamp2(c(-2, 0, 2), c("deepskyblue2", "black", "gold")),
  column_order = colnames(combinedMat),
  heatmap_legend_param = list(title = "Expression", at = c(-2, 0, 2)),
  show_column_dend = FALSE,
  show_row_dend    = TRUE,
  show_row_names   = FALSE,
  column_names_rot = 0
)

# Save outputs ----
pdf(output_pdf)
draw(hm)
dev.off()

message("Heatmap saved: ", output_pdf)

sessionInfo()

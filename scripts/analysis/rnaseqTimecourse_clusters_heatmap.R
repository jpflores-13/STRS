## RNA-seq Timecourse Heatmap Generation
## This script generates a heatmap of significantly changing genes across time
## Output: Heatmap visualization without clustering splits

library(DESeq2)
library(DEGreport)
library(tibble)
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(stringr)

## Load scripts and data
source("scripts/utils/make_norm_matrix.R")
dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

## Get time comparison coefficients
coef <- resultsNames(dds) |> 
  str_subset("Time")

## Perform variance stabilizing shrinkage on log fold changes
## This helps reduce noise in fold change estimates for genes with low counts
message("Performing log fold change shrinkage...")
resShrunk <- list()
for (i in seq_along(coef)) {
  df <- lfcShrink(dds, coef = coef[i], type = "apeglm")
  df$timepoint <- coef[i]
  resShrunk[[i]] <- df
}
shrunkenLFCs <- do.call(rbind, resShrunk)

## Identify significantly changing genes
## Filter for genes with strong fold changes (>2 or <-2) and significant p-values
message("Identifying significantly changing genes...")
sigLRT_genes <- shrunkenLFCs |> 
  data.frame() |> 
  mutate(
    timepoint = case_when(
      timepoint == "Time_1h_vs_0h" ~ "1h",
      timepoint == "Time_3h_vs_0h" ~ "3h",
      timepoint == "Time_6h_vs_0h" ~ "6h",
      timepoint == "Time_9h_vs_0h" ~ "9h",
      timepoint == "Time_12h_vs_0h" ~ "12h",
      timepoint == "Time_24h_vs_0h" ~ "24h"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 2 ~ "upregulated",
      padj < 0.05 & log2FoldChange < -2 ~ "downregulated")) |>
  filter(regulation %in% c("upregulated", "downregulated")) |>
  rownames_to_column(var = "gene") |> 
  pull(gene)

message(paste("Found", length(sigLRT_genes), "significantly changing genes"))

## Create normalized count matrix for visualization
message("Creating normalized count matrix...")
count_matrix_norm <- make_norm_matrix(dds)

# SUBSET: keep only significant time-varying genes
count_matrix_norm <- count_matrix_norm[
  rownames(count_matrix_norm) %in% sigLRT_genes,]

## Combine biological replicates by averaging
message("Combining biological replicates...")
combinedMat <- matrix(NA, nrow = nrow(count_matrix_norm), ncol = 7)
rownames(combinedMat) <- rownames(count_matrix_norm)

for (i in seq(1, (ncol(count_matrix_norm) / 2) )) {
  combinedMat[, i] <- (count_matrix_norm[, 2*i-1] + count_matrix_norm[, 2*i]) / 2
}

## Name columns with timepoints
colnames(combinedMat) <- unique(as.vector(colData(dds)$Time))

## Create metadata for visualization
meta <- data.frame(
  Time = factor(c("0h", "1h", "3h", "6h", "9h", "12h", "24h"),
                levels = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))
)
rownames(meta) <- c("0h", "1h", "3h", "6h", "9h", "12h", "24h")

# ## Determine optimal number of clusters using elbow method
# message("Performing elbow analysis for optimal cluster number...")
# wss <- numeric(15)
# for (k in 1:15) {
#   km <- kmeans(t(scale(t(combinedMat))), centers = k, nstart = 25)
#   wss[k] <- km$tot.withinss
# }
# 
# ## Create and save elbow plot
# message("Creating elbow plot...")
# pdf(file = "plots/elbow_plot.pdf")
# plot(1:15, wss, type = "b", 
#      xlab = "Number of Clusters (k)", 
#      ylab = "Total Within-cluster Sum of Squares",
#      main = "Elbow Plot for K-means Clustering")
# dev.off()
# 
# ## Perform k-means clustering
# message("Performing k-means clustering...")
# k <- 3  # Adjust based on your elbow plot
# set.seed(123)
# km <- kmeans(t(scale(t(combinedMat))), centers = k, nstart = 25)
# 
# ## Prepare data for visualization
# clustered_data <- data.frame(
#   cluster = factor(km$cluster),
#   row.names = rownames(combinedMat)
# )
# 
# ## Order rows by cluster
# row_order <- order(clustered_data$cluster)
# ordered_mat <- combinedMat[row_order,]

## Create and save unified heatmap without cluster splits
message("Creating heatmap...")
pdf(file = "plots/rnaseqTimecourse-heatmap.pdf")
Heatmap(combinedMat,
        col = colorRamp2(
          c(-2, 0, 2), 
          c("deepskyblue2", "black", "gold")
        ),
        column_order = colnames(combinedMat),
        heatmap_legend_param = list(
          title = "Expression", 
          at = c(-2, 0, 2)
        ),
        show_column_dend = FALSE,
        show_row_dend = TRUE,  # Show row dendrogram for hierarchical clustering visualization
        show_row_names = FALSE,
        column_names_rot = 0
)
dev.off()

# ## Save clustering results for downstream GO analysis
# message("Saving clustering results for downstream analysis...")
# 
# # Create list of genes for each cluster
# cluster_genes <- split(rownames(ordered_mat), km$cluster[row_order])
# 
# # Save clustering data
# clustering_results <- list(
#   kmeans_object = km,
#   cluster_genes = cluster_genes,
#   ordered_matrix = ordered_mat,
#   row_order = row_order,
#   clustered_data = clustered_data,
#   metadata = list(
#     k_clusters = k,
#     total_genes = nrow(combinedMat),
#     cluster_sizes = sapply(cluster_genes, length)
#   )
# )
# 
# saveRDS(clustering_results, 
#         file = "data/processed/rna/timecourse/output/clustering_results.rds")
# 
# ## Create summary report
# message("\n=== Clustering Analysis Summary ===")
# message(paste("Total genes analyzed:", nrow(combinedMat)))
# message(paste("Number of clusters:", k))
# for(i in seq_along(cluster_genes)) {
#   message(sprintf("Cluster %d: %d genes", i, length(cluster_genes[[i]])))
# }
# 
# ## Save cluster gene lists as separate files for easy access
# for(i in seq_along(cluster_genes)) {
#   write.csv(data.frame(gene_id = cluster_genes[[i]]),
#             file = paste0("tables/cluster_", i, "_genes.csv"),
#             row.names = FALSE)
# }

message("\nHeatmap generation completed!")
message("- Heatmap saved: plots/rnaseqTimecourse-heatmap.pdf")
## GO Enrichment Analysis and Visualization for K-means Clusters

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggrepel)
library(tidyr)

## Load clustering results from previous script
clustering_data <- readRDS("data/processed/rna/timecourse/output/clustering_results.rds")

# Extract components
cluster_genes <- clustering_data$cluster_genes
k_clusters <- clustering_data$metadata$k_clusters
cluster_sizes <- clustering_data$metadata$cluster_sizes

message("Clustering data loaded successfully!")
message(paste("Number of clusters:", k_clusters))
for(i in seq_along(cluster_genes)) {
  message(sprintf("Cluster %d: %d genes", i, length(cluster_genes[[i]])))
}

## Define color palette for clusters (matching heatmap colors)
cluster_colors <- c("1" = "#F8766D", "2" = "#619CFF", "3" = "#A6D854")

## Enhanced function to perform GO enrichment with advanced visualizations
run_cluster_go_analysis <- function(genes, cluster_number, all_genes) {
  
  # Convert gene IDs from ENSEMBL to ENTREZ
  message(sprintf("Converting gene IDs for Cluster %d...", cluster_number))
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = genes,
                       keytype = "ENSEMBL",
                       column = "ENTREZID")
  
  # Remove NA values
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Log conversion results
  message(sprintf("Cluster %d: Mapped %d out of %d genes to ENTREZ IDs",
                  cluster_number, length(entrez_ids), length(genes)))
  
  # Skip if not enough genes
  if (length(entrez_ids) < 10) {
    message(paste0("Not enough genes in cluster ", cluster_number, " for GO analysis"))
    return(NULL)
  }
  
  # Generate GO enrichment analysis
  message(paste0("Running GO analysis for cluster ", cluster_number, " genes..."))
  ego <- enrichGO(gene = entrez_ids,
                  universe = all_genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  minGSSize = 10,
                  maxGSSize = 1000,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
  
  # Check if any results were found
  if (is.null(ego) || nrow(ego@result) == 0) {
    message(paste0("No enriched GO terms found for cluster ", cluster_number, " genes"))
    return(NULL)
  }
  
  # Simplify terms to reduce redundancy
  ego_simplified <- clusterProfiler::simplify(ego,
                                              cutoff = 0.7,
                                              by = "p.adjust",
                                              select_fun = min)
  
  # Save the results
  saveRDS(ego_simplified,
          paste0("data/processed/rna/timecourse/output/GOobjs/cluster_", cluster_number, "_GO_results.rds"))
  
  # Create enhanced visualizations if results exist
  if (nrow(ego_simplified@result) > 0) {
    # Get top 20 terms or all if less than 20
    top_terms <- ego_simplified@result |>
      slice_min(order_by = p.adjust, n = 20) |>
      as.data.frame()
    
    # Create a shortened version of the Description for display on bars
    top_terms$short_desc <- str_wrap(top_terms$Description, width = 40)
    
    # Create barplot with text labels anchored at the left side of the bars
    pdf(paste0("plots/GO_barplot_cluster_", cluster_number, ".pdf"), width = 12, height = 10)
    
    p <- ggplot(top_terms, aes(x = reorder(short_desc, -p.adjust), y = -log10(p.adjust))) +
      geom_col(fill = cluster_colors[as.character(cluster_number)], alpha = 0.8) +
      geom_text(aes(y = 0.2, label = short_desc), 
                color = "white", fontface = "bold", hjust = 0, size = 3) +
      labs(
        # title = paste0("Top Enriched GO Terms in Cluster ", cluster_number),
        x = "GO Terms",
        y = "-log10(adjusted p-value)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
        plot.title = element_text(size = 12, face = "bold")
      ) +
      scale_y_continuous(expand = c(0, 0)) +  # Start bars at zero
      coord_flip()
    
    print(p)
    dev.off()
    
    
    # Also create the traditional dotplot for comparison
    pdf(paste0("plots/cluster_", cluster_number, "_GO_dotplot.pdf"))
    print(dotplot(ego_simplified, 
                  showCategory = 20) +
            theme(axis.text.y = element_text(size = 8)))
            # ggtitle(paste0("GO Enrichment - Cluster ", cluster_number)))
    dev.off()
    
    # Save detailed results
    write.csv(as.data.frame(ego_simplified),
              file = paste0("tables/cluster_", cluster_number, "_GO_enrichment.csv"),
              row.names = FALSE)
  }
  
  return(ego_simplified)
}

## Function to create focused visualizations (top N terms)
create_focused_go_plots <- function(go_result, cluster_number, top_n, plot_suffix = "") {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    message(paste0("No GO results available for focused plots - Cluster ", cluster_number))
    return(NULL)
  }
  
  # Get top N terms and convert to data frame
  top_terms <- go_result@result |>
    slice_min(order_by = p.adjust, n = top_n) |>
    as.data.frame()
  
  if (nrow(top_terms) == 0) {
    message(paste0("No terms available for top ", top_n, " plot - Cluster ", cluster_number))
    return(NULL)
  }
  
  # Create shortened descriptions for better display
  top_terms$short_desc <- str_wrap(top_terms$Description, width = 50)
  
  # Create focused barplot
  pdf(paste0("plots/GO_barplot_cluster_", cluster_number, "_top", top_n, plot_suffix, ".pdf"), 
      width = 12, height = max(6, top_n * 0.4))
  
  p <- ggplot(top_terms, aes(x = reorder(short_desc, -p.adjust), y = -log10(p.adjust))) +
    geom_col(fill = cluster_colors[as.character(cluster_number)], alpha = 0.8) +
    geom_text(aes(y = 0.2, label = short_desc), 
              color = "white", fontface = "bold", hjust = 0, size = 3.5) +
    labs(
      # title = paste0("Top ", top_n, " Enriched GO Terms in Cluster ", cluster_number),
      x = "GO Terms",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +  # Add space at top for p-value labels
    coord_flip()
  
  print(p)
  dev.off()
  
  message(paste0("Created focused plots for top ", top_n, " terms - Cluster ", cluster_number))
}

## Main GO Analysis Execution

# Get all genes from human genome for background
message("Preparing GO analysis background...")
all_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")

## Run enhanced GO analysis for each cluster
message("Starting GO enrichment analysis...")
go_results <- list()
for(i in seq_along(cluster_genes)) {
  go_results[[i]] <- run_cluster_go_analysis(cluster_genes[[i]], i, all_genes)
}

## Create enhanced summary of results
message("Creating analysis summary...")
enrichment_summary <- data.frame(
  Cluster = integer(),
  Num_Genes = integer(),
  Num_Enriched_Terms = integer(),
  Top_Terms = character(),
  stringsAsFactors = FALSE
)

for(i in seq_along(go_results)) {
  top_terms <- if(!is.null(go_results[[i]])) {
    paste(head(go_results[[i]]@result$Description, 5), collapse = "; ")
  } else "No enriched terms"
  
  enrichment_summary <- rbind(enrichment_summary,
                              data.frame(
                                Cluster = i,
                                Num_Genes = length(cluster_genes[[i]]),
                                Num_Enriched_Terms = 
                                  if(!is.null(go_results[[i]])) nrow(go_results[[i]]@result) else 0,
                                Top_Terms = top_terms
                              ))
}

## Create focused visualizations for Cluster 1 (top 5 and top 10 GO terms)
message("Creating focused visualizations for Cluster 1...")
if (!is.null(go_results[[1]])) {
  # Top 5 GO terms for Cluster 1
  create_focused_go_plots(go_results[[1]], 1, 5)
  
  # Top 10 GO terms for Cluster 1
  create_focused_go_plots(go_results[[1]], 1, 10)
  
  message("\n=== Focused Cluster 1 Visualizations Created ===")
  message("- Top 5 GO terms: GO_barplot_cluster_1_top5.pdf & GO_geneCount_cluster_1_top5.pdf")
  message("- Top 10 GO terms: GO_barplot_cluster_1_top10.pdf & GO_geneCount_cluster_1_top10.pdf")
} else {
  message("No GO results available for Cluster 1 focused plots")
}

## Print final summary to console
message("\n=== GO Enrichment Analysis Summary ===")
for(i in seq_along(go_results)) {
  message(sprintf("Cluster %d: %d genes, %d enriched GO terms", 
                  i, 
                  length(cluster_genes[[i]]),
                  if(!is.null(go_results[[i]])) nrow(go_results[[i]]@result) else 0))
}

message("\nGO enrichment analysis completed successfully!")
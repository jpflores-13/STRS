## Helper functions for exploring RNA-seq timecourse GO enrichment results

library(clusterProfiler)
library(dplyr)
library(ggplot2)

## Function to compare GO terms across timepoints ----
compare_go_across_timepoints <- function(timepoints = c("1h", "3h", "6h", "9h", "12h", "24h"),
                                         direction = "upregulated") {
  
  # Load all GO results
  go_list <- list()
  for (tp in timepoints) {
    file_path <- paste0("data/processed/rna/timecourse/output/GOobjs/", 
                        tp, "_", direction, "_GO_results.rds")
    
    if (file.exists(file_path)) {
      go_list[[tp]] <- readRDS(file_path)
    }
  }
  
  # Extract GO terms for each timepoint
  go_terms <- lapply(names(go_list), function(tp) {
    if (!is.null(go_list[[tp]])) {
      go_list[[tp]]@result |>
        mutate(timepoint = tp) |>
        select(timepoint, ID, Description, p.adjust, Count, GeneRatio)
    }
  })
  
  # Combine into single dataframe
  go_combined <- bind_rows(go_terms)
  
  return(go_combined)
}

## Function to identify persistent GO terms ----
find_persistent_terms <- function(direction = "upregulated", min_timepoints = 3) {
  
  # Get all GO terms across timepoints
  all_terms <- compare_go_across_timepoints(direction = direction)
  
  # Count how many timepoints each term appears in
  persistent <- all_terms |>
    group_by(ID, Description) |>
    summarize(
      n_timepoints = n(),
      timepoints = paste(timepoint, collapse = ", "),
      mean_padj = mean(p.adjust),
      .groups = "drop"
    ) |>
    filter(n_timepoints >= min_timepoints) |>
    arrange(desc(n_timepoints), mean_padj)
  
  return(persistent)
}

## Function to extract genes from specific GO term ----
get_genes_in_go_term <- function(timepoint, go_id, direction = "upregulated") {
  
  # Load GO results
  file_path <- paste0("data/processed/rna/timecourse/output/GOobjs/", 
                      timepoint, "_", direction, "_GO_results.rds")
  
  if (!file.exists(file_path)) {
    stop(sprintf("GO results not found for %s %s genes", timepoint, direction))
  }
  
  ego <- readRDS(file_path)
  
  # Get gene list for specified GO term
  term_data <- ego@result |>
    filter(ID == go_id)
  
  if (nrow(term_data) == 0) {
    stop(sprintf("GO term %s not found in %s %s results", go_id, timepoint, direction))
  }
  
  # Extract gene symbols
  gene_entrez <- unlist(strsplit(term_data$geneID, "/"))
  
  # Convert to symbols
  gene_symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_entrez,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )
  
  return(gene_symbols)
}

## Function to create heatmap of GO term significance across timepoints ----
plot_go_term_heatmap <- function(direction = "upregulated", top_n = 20) {
  
  library(pheatmap)
  library(RColorBrewer)
  
  # Get all GO terms
  all_terms <- compare_go_across_timepoints(direction = direction)
  
  # Identify top terms by average significance
  top_terms <- all_terms |>
    group_by(ID, Description) |>
    summarize(
      mean_padj = mean(p.adjust),
      n_timepoints = n(),
      .groups = "drop"
    ) |>
    arrange(mean_padj) |>
    head(top_n) |>
    pull(ID)
  
  # Create matrix for heatmap
  go_matrix <- all_terms |>
    filter(ID %in% top_terms) |>
    mutate(neg_log10_padj = -log10(p.adjust)) |>
    select(timepoint, Description, neg_log10_padj) |>
    tidyr::pivot_wider(
      names_from = timepoint,
      values_from = neg_log10_padj,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("Description") |>
    as.matrix()
  
  # Reorder columns
  timepoint_order <- c("1h", "3h", "6h", "9h", "12h", "24h")
  go_matrix <- go_matrix[, intersect(timepoint_order, colnames(go_matrix))]
  
  # Create heatmap
  pheatmap(
    go_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "firebrick"))(100),
    main = paste0("GO Term Enrichment - ", direction, " genes"),
    fontsize_row = 8,
    angle_col = 0
  )
}

## Function to summarize differential expression across all timepoints ----
summarize_timecourse_de <- function() {
  
  dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
  
  timepoints <- c("1h", "3h", "6h", "9h", "12h", "24h")
  coef_names <- resultsNames(dds) |>
    stringr::str_subset("Time")
  
  summary_df <- data.frame()
  
  for (i in seq_along(timepoints)) {
    res <- results(dds, name = coef_names[i])
    
    n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 2, na.rm = TRUE)
    n_down <- sum(res$padj < 0.05 & res$log2FoldChange < -2, na.rm = TRUE)
    
    summary_df <- rbind(
      summary_df,
      data.frame(
        timepoint = timepoints[i],
        upregulated = n_up,
        downregulated = n_down,
        total = n_up + n_down
      )
    )
  }
  
  return(summary_df)
}

## Function to compare GO enrichment with different background sets ----
compare_background_impact <- function(timepoint, direction = "upregulated", 
                                      gene_list_ensembl, 
                                      min_basemean = 10) {
  
  library(org.Hs.eg.db)
  
  message(sprintf("\nComparing GO enrichment backgrounds for %s %s genes", 
                  timepoint, direction))
  
  # Convert gene list to ENTREZ
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = gene_list_ensembl,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Background 1: All human genes
  bg1 <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  message(sprintf("  Background 1 (all human genes): %d genes", length(bg1)))
  
  ego1 <- enrichGO(
    gene = entrez_ids,
    universe = bg1,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # Background 2: All expressed genes (from DESeq2 object)
  dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
  expressed_ensembl <- rownames(dds)
  bg2 <- mapIds(
    org.Hs.eg.db,
    keys = expressed_ensembl,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  bg2 <- unique(bg2[!is.na(bg2)])
  message(sprintf("  Background 2 (all expressed genes): %d genes", length(bg2)))
  
  ego2 <- enrichGO(
    gene = entrez_ids,
    universe = bg2,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # Background 3: Well-expressed genes (baseMean >= threshold)
  basemeans <- rowMeans(counts(dds, normalized = TRUE))
  well_expressed_ensembl <- names(basemeans[basemeans >= min_basemean])
  bg3 <- mapIds(
    org.Hs.eg.db,
    keys = well_expressed_ensembl,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  bg3 <- unique(bg3[!is.na(bg3)])
  message(sprintf("  Background 3 (expressed with baseMean >= %.1f): %d genes", 
                  min_basemean, length(bg3)))
  
  ego3 <- enrichGO(
    gene = entrez_ids,
    universe = bg3,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # Compare results
  comparison <- list(
    all_human = list(
      n_genes = length(bg1),
      n_terms = if(!is.null(ego1)) nrow(ego1@result) else 0,
      top_terms = if(!is.null(ego1) && nrow(ego1@result) > 0) {
        head(ego1@result$Description, 10)
      } else NULL
    ),
    all_expressed = list(
      n_genes = length(bg2),
      n_terms = if(!is.null(ego2)) nrow(ego2@result) else 0,
      top_terms = if(!is.null(ego2) && nrow(ego2@result) > 0) {
        head(ego2@result$Description, 10)
      } else NULL
    ),
    well_expressed = list(
      n_genes = length(bg3),
      n_terms = if(!is.null(ego3)) nrow(ego3@result) else 0,
      top_terms = if(!is.null(ego3) && nrow(ego3@result) > 0) {
        head(ego3@result$Description, 10)
      } else NULL
    )
  )
  
  message(sprintf("\nResults:"))
  message(sprintf("  All human genes:        %d enriched terms", 
                  comparison$all_human$n_terms))
  message(sprintf("  All expressed genes:    %d enriched terms", 
                  comparison$all_expressed$n_terms))
  message(sprintf("  Well-expressed genes:   %d enriched terms", 
                  comparison$well_expressed$n_terms))
  
  # Find terms unique to each approach
  if (!is.null(ego1) && !is.null(ego2) && !is.null(ego3) &&
      nrow(ego1@result) > 0 && nrow(ego2@result) > 0 && nrow(ego3@result) > 0) {
    terms1 <- ego1@result$ID
    terms2 <- ego2@result$ID
    terms3 <- ego3@result$ID
    
    shared_all <- Reduce(intersect, list(terms1, terms2, terms3))
    
    message(sprintf("\n  Shared across all three approaches: %d terms", 
                    length(shared_all)))
    message(sprintf("  Unique to approach 1: %d", 
                    length(setdiff(terms1, union(terms2, terms3)))))
    message(sprintf("  Unique to approach 2: %d", 
                    length(setdiff(terms2, union(terms1, terms3)))))
    message(sprintf("  Unique to approach 3 (RECOMMENDED): %d", 
                    length(setdiff(terms3, union(terms1, terms2)))))
  }
  
  return(comparison)
}

## Print usage instructions ----
## Print usage instructions ----
message("Helper functions loaded! Available functions:")
message("  - compare_go_across_timepoints(): Compare GO terms across timepoints")
message("  - find_persistent_terms(): Find GO terms enriched at multiple timepoints")
message("  - get_genes_in_go_term(): Extract genes belonging to a specific GO term")
message("  - plot_go_term_heatmap(): Create heatmap of GO term significance")
message("  - summarize_timecourse_de(): Get summary of DE genes per timepoint")
message("  - compare_background_impact(): Compare GO results with different backgrounds")
message("\nExample usage:")
message('  persistent <- find_persistent_terms(direction = "upregulated", min_timepoints = 4)')
message('  genes <- get_genes_in_go_term(timepoint = "6h", go_id = "GO:0001525", direction = "upregulated")')
message('  # Compare impact of background choice:')
message('  dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")')
message('  res <- results(dds, name = "Time_6h_vs_0h")')
message('  upreg <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 2]')
message('  bg_comparison <- compare_background_impact("6h", "upregulated", upreg)')
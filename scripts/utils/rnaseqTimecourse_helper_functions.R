# rnaseqTimecourse_helper_functions.R ------------------------------------ ## Helper functions for RNA-seq timecourse GO enrichment analysis

# Libraries ----
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)
library(tidyr)
library(tibble)

# Functions ----

#' Compare GO terms across timepoints
#'
#' Loads pre-saved GO enrichment RDS files for each timepoint and returns a
#' combined data frame of enriched terms.
#'
#' @param timepoints character. Vector of timepoint labels (default: c("1h","3h","6h","9h","12h","24h"))
#' @param direction character. Either "upregulated" or "downregulated" (default: "upregulated")
#' @param go_dir character. Directory containing per-timepoint GO result RDS files
#' @return data.frame. Combined GO results with columns: timepoint, ID, Description,
#'   p.adjust, Count, GeneRatio
compare_go_across_timepoints <- function(timepoints = c("1h", "3h", "6h", "9h", "12h", "24h"),
                                         direction = "upregulated",
                                         go_dir) {

  # Load all GO results
  go_list <- list()
  for (tp in timepoints) {
    file_path <- file.path(go_dir, paste0(tp, "_", direction, "_GO_results.rds"))

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

#' Identify GO terms enriched at multiple timepoints
#'
#' @param direction character. Either "upregulated" or "downregulated" (default: "upregulated")
#' @param min_timepoints integer. Minimum number of timepoints a term must appear in (default: 3)
#' @param go_dir character. Directory containing per-timepoint GO result RDS files
#' @return data.frame. Persistent GO terms with columns: ID, Description, n_timepoints,
#'   timepoints, mean_padj; sorted by n_timepoints then mean_padj
find_persistent_terms <- function(direction = "upregulated", min_timepoints = 3, go_dir) {

  # Get all GO terms across timepoints
  all_terms <- compare_go_across_timepoints(direction = direction, go_dir = go_dir)

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

#' Extract genes belonging to a specific GO term at one timepoint
#'
#' @param timepoint character. Timepoint label (e.g., "6h")
#' @param go_id character. GO term ID (e.g., "GO:0001525")
#' @param direction character. Either "upregulated" or "downregulated" (default: "upregulated")
#' @param go_dir character. Directory containing per-timepoint GO result RDS files
#' @return data.frame. Gene symbols and ENTREZ IDs for genes in the specified GO term
get_genes_in_go_term <- function(timepoint, go_id, direction = "upregulated", go_dir) {

  # Load GO results
  file_path <- file.path(go_dir, paste0(timepoint, "_", direction, "_GO_results.rds"))

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

#' Create heatmap of GO term significance across timepoints
#'
#' @param direction character. Either "upregulated" or "downregulated" (default: "upregulated")
#' @param top_n integer. Number of top terms to include (default: 20)
#' @param go_dir character. Directory containing per-timepoint GO result RDS files
#' @return NULL (invisibly); prints a pheatmap to the active device
plot_go_term_heatmap <- function(direction = "upregulated", top_n = 20, go_dir) {

  # Get all GO terms
  all_terms <- compare_go_across_timepoints(direction = direction, go_dir = go_dir)

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

#' Summarize differential expression counts across all timepoints
#'
#' @param deseq_rds character. Path to the LRT DESeq2 RDS file
#' @return data.frame. Columns: timepoint, upregulated, downregulated, total
summarize_timecourse_de <- function(deseq_rds) {

  dds <- readRDS(deseq_rds)

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

#' Compare GO enrichment results using different background gene sets
#'
#' Runs GO enrichment three ways — all human genes, all expressed genes,
#' and well-expressed genes — and reports overlap and unique terms.
#'
#' @param timepoint character. Timepoint label (e.g., "6h")
#' @param direction character. Either "upregulated" or "downregulated" (default: "upregulated")
#' @param gene_list_ensembl character. Vector of ENSEMBL gene IDs to test
#' @param deseq_rds character. Path to the LRT DESeq2 RDS file
#' @param min_basemean numeric. baseMean threshold for the well-expressed background (default: 10)
#' @return list. Named list with elements all_human, all_expressed, well_expressed;
#'   each containing n_genes, n_terms, and top_terms
compare_background_impact <- function(timepoint, direction = "upregulated",
                                      gene_list_ensembl,
                                      deseq_rds,
                                      min_basemean = 10) {

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
  dds <- readRDS(deseq_rds)
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

# Footer ----
message("rnaseqTimecourse_helper_functions.R loaded. Available functions:")
message("  compare_go_across_timepoints(timepoints, direction, go_dir)")
message("  find_persistent_terms(direction, min_timepoints, go_dir)")
message("  get_genes_in_go_term(timepoint, go_id, direction, go_dir)")
message("  plot_go_term_heatmap(direction, top_n, go_dir)")
message("  summarize_timecourse_de(deseq_rds)")
message("  compare_background_impact(timepoint, direction, gene_list_ensembl, deseq_rds, min_basemean)")
message("Example:")
message('  go_dir <- "data/processed/rna/timecourse/output/GOobjs"')
message('  persistent <- find_persistent_terms(direction = "upregulated", min_timepoints = 4, go_dir = go_dir)')
message('  genes <- get_genes_in_go_term("6h", "GO:0001525", go_dir = go_dir)')

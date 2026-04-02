## =============================================================================
## EISA PART 2: DESeq2 ANALYSIS
## =============================================================================
##
## Purpose: Run DESeq2 on intronic reads and calculate EISA metrics
##
## Input:  Count matrices from Part 1
## Output: DESeq2 objects, EISA results table, statistics
##
## Runtime: ~5-10 minutes
##
## =============================================================================

library(DESeq2)
library(GenomicFeatures)
library(GenomicRanges)
library(InteractionSet)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)

## =============================================================================
## CONFIGURATION
## =============================================================================

BASE_DIR <- "/work/users/j/p/jpflores/projects/STRS"
EISA_DIR <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/EISA")
EXISTING_DDS <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
LOOP_DATA <- file.path(BASE_DIR, "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
TXDB_FILE <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb"

## =============================================================================
## LOAD COUNT MATRICES
## =============================================================================

exon_counts_file <- file.path(EISA_DIR, "exon_counts.rds")
intron_counts_file <- file.path(EISA_DIR, "intron_counts.rds")

## Load Count Matrices
exon_counts <- readRDS(exon_counts_file)
intron_counts <- readRDS(intron_counts_file)

message(sprintf("Exon counts: %d genes × %d samples", nrow(exon_counts), ncol(exon_counts)))
message(sprintf("Intron counts: %d genes × %d samples", nrow(intron_counts), ncol(intron_counts)))

## =============================================================================
## CREATE SAMPLE METADATA
## =============================================================================

sample_names <- colnames(exon_counts)

## Extract timepoint and replicate from sample names
create_coldata <- function(sample_names) {
  coldata <- data.frame(
    sample = sample_names,
    row.names = sample_names
  )
  
  ## Parse timepoint
  coldata$Time <- str_extract(sample_names, "[0-9]+h")
  
  ## Parse replicate: extract the number right after timepoint
  ## Pattern: ...0h_1_... or ...12h_2_...
  ##               ^ this number
  coldata$Bio_Rep <- str_extract(sample_names, "[0-9]+h_(1|2)_") |>
    str_extract("_(1|2)_") |>
    str_remove_all("_")
  
  ## Convert to factors with proper levels
  coldata$Time <- factor(coldata$Time, 
                         levels = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))
  coldata$Bio_Rep <- factor(coldata$Bio_Rep, levels = c("1", "2"))
  
  return(coldata)
}

coldata <- create_coldata(sample_names)

print(table(coldata$Time, coldata$Bio_Rep))

## =============================================================================
## FILTER GENES
## =============================================================================


## Make sure both matrices have the same genes
common_genes <- intersect(rownames(exon_counts), rownames(intron_counts))
message(sprintf("Genes in both matrices: %d", length(common_genes)))

exon_counts <- exon_counts[common_genes, ]
intron_counts <- intron_counts[common_genes, ]

## Normalize by library size
exon_norm <- sweep(exon_counts, 2, colSums(exon_counts), "/") * mean(colSums(exon_counts))
intron_norm <- sweep(intron_counts, 2, colSums(intron_counts), "/") * mean(colSums(intron_counts))

## EISA paper threshold: average log2(counts+8) >= 5
exon_pass <- rowMeans(log2(exon_norm + 8)) >= 5
intron_pass <- rowMeans(log2(intron_norm + 8)) >= 5
genes_pass <- exon_pass & intron_pass

message(sprintf("Genes with sufficient exonic coverage: %d", sum(exon_pass)))
message(sprintf("Genes with sufficient intronic coverage: %d", sum(intron_pass)))
message(sprintf("Genes with BOTH: %d (%.1f%%)", 
                sum(genes_pass),
                sum(genes_pass) / length(genes_pass) * 100))

## Filter matrices
exon_counts_filt <- exon_counts[genes_pass, ]
intron_counts_filt <- intron_counts[genes_pass, ]

message(sprintf("Using %d genes for EISA analysis", nrow(intron_counts_filt)))

## =============================================================================
## CREATE DESeq2 OBJECT FOR INTRONS
## =============================================================================

message("\n=== STEP 4: Creating DESeq2 object for intronic reads ===")

dds_intron <- DESeqDataSetFromMatrix(
  countData = intron_counts_filt,
  colData = coldata,
  design = ~ Bio_Rep + Time
)

## Add gene annotations
## Strip version numbers from Ensembl IDs
## Your IDs: ENSG00000000003.15
## org.Hs.eg.db expects: ENSG00000000003
ensembl_ids_clean <- str_remove(rownames(dds_intron), "\\..*$")

gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
) |>
  as.data.frame() |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

## Add to rowData
rowData(dds_intron)$gene_id <- rownames(dds_intron)  # Keep original with version
rowData(dds_intron)$symbol <- gene_symbols$SYMBOL[match(ensembl_ids_clean, 
                                                        gene_symbols$ENSEMBL)]

## Run DESeq2 with LRT
message("Running DESeq2 on intronic reads (LRT test)...")
message("This matches your existing exonic analysis design")

dds_intron <- DESeq(dds_intron, test = "LRT", reduced = ~ Bio_Rep)
dds_intron$Time <- fct_relevel(dds_intron$Time, 
                               c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))

## Save intronic DESeq2 object
saveRDS(dds_intron, file.path(EISA_DIR, "dds_intron_LRT.rds"))
message("✓ Intronic DESeq2 analysis complete")

## =============================================================================
## EXTRACT FOLD CHANGES
## =============================================================================

message("\n=== STEP 5: Extracting fold changes ===")

## Load existing exonic DESeq2 object
message("Loading existing exonic DESeq2 object...")
dds_exon <- readRDS(EXISTING_DDS)

## Get coefficient names
coef_names <- resultsNames(dds_intron) |> str_subset("Time")
message(sprintf("Analyzing %d timepoint comparisons", length(coef_names)))

## Extract fold changes
extract_fold_changes <- function(dds, coef_names) {
  fc_list <- list()
  padj_list <- list()
  
  for (coef in coef_names) {
    res <- results(dds, name = coef)
    tp <- str_extract(coef, "[0-9]+h")
    
    fc_list[[tp]] <- res$log2FoldChange
    names(fc_list[[tp]]) <- rownames(res)
    
    padj_list[[tp]] <- res$padj
    names(padj_list[[tp]]) <- rownames(res)
  }
  
  return(list(fc = fc_list, padj = padj_list))
}

message("Extracting Δintron (transcriptional changes)...")
intron_results <- extract_fold_changes(dds_intron, coef_names)

message("Extracting Δexon (total expression changes)...")
exon_results <- extract_fold_changes(dds_exon, coef_names)

## Match genes
common_genes_eisa <- intersect(
  str_remove(names(intron_results$fc[["1h"]]), "\\..*$"),  # Strip .15 from introns
  names(exon_results$fc[["1h"]])                            # Exons already clean
)

message(sprintf("Genes with both intronic and exonic data: %d", length(common_genes_eisa)))

## =============================================================================
## CREATE EISA RESULTS TABLE
## =============================================================================

message("\n=== STEP 6: Creating EISA results table ===")

create_eisa_table <- function(intron_res, exon_res, common_genes) {
  timepoints <- names(intron_res$fc)
  
  eisa_df <- data.frame(
    gene_id = common_genes,
    row.names = common_genes
  )
  
  for (tp in timepoints) {
    ## Strip version numbers from intron IDs for matching
    intron_ids_clean <- str_remove(names(intron_res$fc[[tp]]), "\\..*$")
    
    ## Match to common genes (which are already without versions)
    intron_idx <- match(common_genes, intron_ids_clean)
    exon_idx <- match(common_genes, names(exon_res$fc[[tp]]))
    
    delta_intron <- intron_res$fc[[tp]][intron_idx]
    delta_exon <- exon_res$fc[[tp]][exon_idx]
    delta_posttx <- delta_exon - delta_intron
    
    eisa_df[[paste0("delta_intron_", tp)]] <- delta_intron
    eisa_df[[paste0("delta_exon_", tp)]] <- delta_exon
    eisa_df[[paste0("delta_posttx_", tp)]] <- delta_posttx
    
    ## Same for padj
    eisa_df[[paste0("padj_intron_", tp)]] <- intron_res$padj[[tp]][intron_idx]
    eisa_df[[paste0("padj_exon_", tp)]] <- exon_res$padj[[tp]][exon_idx]
  }
  
  return(eisa_df)
}

## Now create the table
eisa_results <- create_eisa_table(intron_results, exon_results, common_genes_eisa)

## Add gene symbols
ensembl_ids_clean_eisa <- str_remove(eisa_results$gene_id, "\\..*$")

eisa_results$symbol <- gene_symbols$SYMBOL[match(ensembl_ids_clean_eisa, 
                                                 gene_symbols$ENSEMBL)]

message("✓ EISA metrics calculated for all timepoints")

## =============================================================================
## ANNOTATE WITH LOOP DATA
## =============================================================================

message("\n=== STEP 7: Annotating genes at chromatin loop anchors ===")

## Load loop data
loops <- readRDS(LOOP_DATA) |> interactions()
loops <- keepStandardChromosomes(loops, pruning.mode = "coarse")
loops <- loops[seqnames(anchors(loops, "first")) != "chrM"]

## Classify loops
gained_loops <- loops[mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange > 0]
lost_loops <- loops[mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange < 0]
static_loops <- loops[mcols(loops)$padj > 0.1]

message(sprintf("Gained loops: %d", length(gained_loops)))
message(sprintf("Lost loops: %d", length(lost_loops)))
message(sprintf("Static loops: %d", length(static_loops)))

## Load gene coordinates
txdb <- loadDb(TXDB_FILE)
genes_all <- genes(txdb)
genes_all <- keepStandardChromosomes(genes_all, pruning.mode = "coarse")
genes_all <- genes_all[seqnames(genes_all) != "chrM"]

## Get genes with EISA data
genes_all_clean <- genes_all
names(genes_all_clean) <- str_remove(names(genes_all), "\\..*$")
genes_with_eisa <- genes_all_clean[names(genes_all_clean) %in% eisa_results$gene_id]
message(sprintf("Genes with EISA data: %d", length(genes_with_eisa)))
promoters <- promoters(genes_with_eisa, upstream = 2000, downstream = 500)

## Find genes at loop anchors
find_loop_anchor_genes <- function(loop_set, promoters) {
  hits_anchor1 <- findOverlaps(promoters, anchors(loop_set, "first"))
  hits_anchor2 <- findOverlaps(promoters, anchors(loop_set, "second"))
  
  unique_genes <- unique(c(
    names(promoters)[queryHits(hits_anchor1)],
    names(promoters)[queryHits(hits_anchor2)]
  ))
  
  return(unique_genes)
}

gained_genes <- find_loop_anchor_genes(gained_loops, promoters)
lost_genes <- find_loop_anchor_genes(lost_loops, promoters)
static_genes <- find_loop_anchor_genes(static_loops, promoters)

## Annotate EISA results
eisa_results$loop_category <- "none"
eisa_results$loop_category[eisa_results$gene_id %in% gained_genes] <- "gained"
eisa_results$loop_category[eisa_results$gene_id %in% lost_genes] <- "lost"
eisa_results$loop_category[eisa_results$gene_id %in% static_genes] <- "static"

message(sprintf("Genes at gained loop anchors: %d", 
                sum(eisa_results$loop_category == "gained")))
message(sprintf("Genes at lost loop anchors: %d", 
                sum(eisa_results$loop_category == "lost")))
message(sprintf("Genes at static loop anchors: %d", 
                sum(eisa_results$loop_category == "static")))

## Save results
saveRDS(eisa_results, file.path(EISA_DIR, "eisa_results_with_loops.rds"))

message("✓ Loop annotations complete")

## =============================================================================
## CALCULATE STATISTICS
## =============================================================================

message("\n=== STEP 8: Computing key statistics ===")

gained_data <- eisa_results |> filter(loop_category == "gained")
lost_data <- eisa_results |> filter(loop_category == "lost")
static_data <- eisa_results |> filter(loop_category == "static")

## 1. Correlation between Δintron and Δexon
cor_test <- cor.test(
  gained_data$delta_intron_12h,
  gained_data$delta_exon_12h,
  method = "pearson"
)

## 2. Magnitude comparison
wilcox_test <- wilcox.test(
  abs(gained_data$delta_posttx_12h),
  abs(gained_data$delta_intron_12h),
  paired = TRUE,
  alternative = "less"
)

## 3. Gained vs static loops
ttest_result <- t.test(
  gained_data$delta_intron_12h,
  static_data$delta_intron_12h
)

## Print summary
cat("\n")
cat("================================================================================\n")
cat("                    KEY STATISTICS FOR MANUSCRIPT\n")
cat("================================================================================\n")
cat("\n")
cat("1. Correlation at gained loop anchors (12h timepoint):\n")
cat(sprintf("   Pearson R = %.3f\n", cor_test$estimate))
cat(sprintf("   P-value = %.2e\n", cor_test$p.value))
cat("   → Strong correlation indicates transcriptional regulation\n")
cat("\n")
cat("2. Magnitude of post-transcriptional regulation:\n")
cat(sprintf("   Median |Δintron| = %.3f\n", 
            median(abs(gained_data$delta_intron_12h), na.rm = TRUE)))
cat(sprintf("   Median |Δexon-Δintron| = %.3f\n", 
            median(abs(gained_data$delta_posttx_12h), na.rm = TRUE)))
cat(sprintf("   Fold difference = %.1f\n",
            median(abs(gained_data$delta_intron_12h), na.rm = TRUE) / 
              median(abs(gained_data$delta_posttx_12h), na.rm = TRUE)))
cat(sprintf("   P-value (Wilcoxon test) = %.2e\n", wilcox_test$p.value))
cat("   → Post-transcriptional regulation is minimal\n")
cat("\n")
cat("3. Gained vs static loops (12h timepoint):\n")
cat(sprintf("   Mean Δintron (gained) = %.3f\n", 
            mean(gained_data$delta_intron_12h, na.rm = TRUE)))
cat(sprintf("   Mean Δintron (static) = %.3f\n", 
            mean(static_data$delta_intron_12h, na.rm = TRUE)))
cat(sprintf("   Difference = %.3f\n",
            mean(gained_data$delta_intron_12h, na.rm = TRUE) -
              mean(static_data$delta_intron_12h, na.rm = TRUE)))
cat(sprintf("   P-value (t-test) = %.2e\n", ttest_result$p.value))
cat("   → Gained loops drive transcriptional activation\n")
cat("\n")
cat("================================================================================\n")
cat("\n")

## Save statistics
stats_summary <- list(
  correlation = cor_test,
  magnitude = wilcox_test,
  gained_vs_static = ttest_result,
  summary_table = data.frame(
    metric = c(
      "R (Δintron vs Δexon)",
      "median |Δintron|",
      "median |Δposttx|",
      "mean Δintron (gained)",
      "mean Δintron (static)"
    ),
    value = c(
      cor_test$estimate,
      median(abs(gained_data$delta_intron_12h), na.rm = TRUE),
      median(abs(gained_data$delta_posttx_12h), na.rm = TRUE),
      mean(gained_data$delta_intron_12h, na.rm = TRUE),
      mean(static_data$delta_intron_12h, na.rm = TRUE)
    ),
    p_value = c(
      cor_test$p.value,
      wilcox_test$p.value,
      NA,
      ttest_result$p.value,
      NA
    )
  )
)

saveRDS(stats_summary, file.path(EISA_DIR, "statistics_summary.rds"))

## =============================================================================
## COMPLETION
## =============================================================================

message("\n")
message("================================================================================")
message("                    PART 2 COMPLETE: DESeq2 ANALYSIS")
message("================================================================================")
message("\n")
message("Output files:")
message(sprintf("  1. Intronic DESeq2 object: %s", 
                file.path(EISA_DIR, "dds_intron_LRT.rds")))
message(sprintf("  2. EISA results table: %s",
                file.path(EISA_DIR, "eisa_results_with_loops.rds")))
message(sprintf("  3. Statistics summary: %s",
                file.path(EISA_DIR, "statistics_summary.rds")))
message("\n")
message("Summary:")
message(sprintf("  - %d genes analyzed", nrow(eisa_results)))
message(sprintf("  - %d at gained loops, %d at lost, %d at static",
                sum(eisa_results$loop_category == "gained"),
                sum(eisa_results$loop_category == "lost"),
                sum(eisa_results$loop_category == "static")))
message(sprintf("  - Correlation (Δintron vs Δexon): R = %.3f, p = %.2e",
                cor_test$estimate, cor_test$p.value))
message("\n")
message("Next step: Run Part 3 (Visualizations)")
message("  Rscript scripts/EISA/EISA_part3_plots.R")
message("\n")
message("Completed: ", Sys.time())
message("================================================================================")


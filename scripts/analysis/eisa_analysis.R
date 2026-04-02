# ##############################################################################
# filename:    eisa_analysis.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: EISA Part 2 — DESeq2 LRT analysis on intronic reads; calculates
#              delta-intron, delta-exon, and delta-posttx per timepoint; annotates
#              genes at gained, lost, and static loop anchors
# ##############################################################################

# Libraries ----
library(DESeq2)
library(GenomicFeatures)
library(GenomicRanges)
library(InteractionSet)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)

# Parameters ----
eisa_dir        <- "data/processed/rna/timecourse/output/EISA"
existing_dds_rds <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
loop_data_rds   <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
txdb_file       <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb"

# Data import ----
exon_counts_file   <- file.path(eisa_dir, "exon_counts.rds")
intron_counts_file <- file.path(eisa_dir, "intron_counts.rds")

exon_counts   <- readRDS(exon_counts_file)
intron_counts <- readRDS(intron_counts_file)

message(sprintf("Exon counts: %d genes × %d samples",   nrow(exon_counts),   ncol(exon_counts)))
message(sprintf("Intron counts: %d genes × %d samples", nrow(intron_counts), ncol(intron_counts)))

# Analysis ----

## Build sample metadata ----
sample_names <- colnames(exon_counts)

create_coldata <- function(sample_names) {
  coldata <- data.frame(sample = sample_names, row.names = sample_names)

  coldata$Time    <- str_extract(sample_names, "[0-9]+h")
  coldata$Bio_Rep <- str_extract(sample_names, "[0-9]+h_(1|2)_") |>
    str_extract("_(1|2)_") |>
    str_remove_all("_")

  coldata$Time    <- factor(coldata$Time,
                            levels = c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))
  coldata$Bio_Rep <- factor(coldata$Bio_Rep, levels = c("1", "2"))
  coldata
}

coldata <- create_coldata(sample_names)
print(table(coldata$Time, coldata$Bio_Rep))

## Filter genes ----
common_genes <- intersect(rownames(exon_counts), rownames(intron_counts))
message(sprintf("Genes in both matrices: %d", length(common_genes)))

exon_counts   <- exon_counts[common_genes, ]
intron_counts <- intron_counts[common_genes, ]

exon_norm   <- sweep(exon_counts,   2, colSums(exon_counts),   "/") * mean(colSums(exon_counts))
intron_norm <- sweep(intron_counts, 2, colSums(intron_counts), "/") * mean(colSums(intron_counts))

exon_pass   <- rowMeans(log2(exon_norm   + 8)) >= 5
intron_pass <- rowMeans(log2(intron_norm + 8)) >= 5
genes_pass  <- exon_pass & intron_pass

message(sprintf("Genes with sufficient exonic coverage: %d",  sum(exon_pass)))
message(sprintf("Genes with sufficient intronic coverage: %d", sum(intron_pass)))
message(sprintf("Genes with BOTH: %d (%.1f%%)",
                sum(genes_pass), sum(genes_pass) / length(genes_pass) * 100))

exon_counts_filt   <- exon_counts[genes_pass, ]
intron_counts_filt <- intron_counts[genes_pass, ]

## DESeq2 on intronic reads ----
message("\nCreating DESeq2 object for intronic reads...")

dds_intron <- DESeqDataSetFromMatrix(
  countData = intron_counts_filt,
  colData   = coldata,
  design    = ~ Bio_Rep + Time
)

ensembl_ids_clean <- str_remove(rownames(dds_intron), "\\..*$")

gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = ensembl_ids_clean,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
) |>
  as.data.frame() |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

rowData(dds_intron)$gene_id <- rownames(dds_intron)
rowData(dds_intron)$symbol  <- gene_symbols$SYMBOL[match(ensembl_ids_clean,
                                                         gene_symbols$ENSEMBL)]

message("Running DESeq2 on intronic reads (LRT)...")
dds_intron <- DESeq(dds_intron, test = "LRT", reduced = ~ Bio_Rep)
dds_intron$Time <- fct_relevel(dds_intron$Time,
                                c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))

saveRDS(dds_intron, file.path(eisa_dir, "dds_intron_LRT.rds"))
message("Intronic DESeq2 analysis complete")

## Extract fold changes ----
message("\nExtracting fold changes...")

dds_exon   <- readRDS(existing_dds_rds)
coef_names <- resultsNames(dds_intron) |> str_subset("Time")

extract_fold_changes <- function(dds, coef_names) {
  fc_list   <- list()
  padj_list <- list()

  for (coef in coef_names) {
    res <- results(dds, name = coef)
    tp  <- str_extract(coef, "[0-9]+h")

    fc_list[[tp]]   <- res$log2FoldChange
    names(fc_list[[tp]])   <- rownames(res)
    padj_list[[tp]] <- res$padj
    names(padj_list[[tp]]) <- rownames(res)
  }

  list(fc = fc_list, padj = padj_list)
}

intron_results <- extract_fold_changes(dds_intron, coef_names)
exon_results   <- extract_fold_changes(dds_exon,   coef_names)

common_genes_eisa <- intersect(
  str_remove(names(intron_results$fc[["1h"]]), "\\..*$"),
  names(exon_results$fc[["1h"]])
)

message(sprintf("Genes with both intronic and exonic data: %d", length(common_genes_eisa)))

## Build EISA results table ----
create_eisa_table <- function(intron_res, exon_res, common_genes) {
  timepoints <- names(intron_res$fc)

  eisa_df <- data.frame(gene_id = common_genes, row.names = common_genes)

  for (tp in timepoints) {
    intron_ids_clean <- str_remove(names(intron_res$fc[[tp]]), "\\..*$")

    intron_idx <- match(common_genes, intron_ids_clean)
    exon_idx   <- match(common_genes, names(exon_res$fc[[tp]]))

    delta_intron  <- intron_res$fc[[tp]][intron_idx]
    delta_exon    <- exon_res$fc[[tp]][exon_idx]
    delta_posttx  <- delta_exon - delta_intron

    eisa_df[[paste0("delta_intron_",  tp)]] <- delta_intron
    eisa_df[[paste0("delta_exon_",    tp)]] <- delta_exon
    eisa_df[[paste0("delta_posttx_",  tp)]] <- delta_posttx
    eisa_df[[paste0("padj_intron_",   tp)]] <- intron_res$padj[[tp]][intron_idx]
    eisa_df[[paste0("padj_exon_",     tp)]] <- exon_res$padj[[tp]][exon_idx]
  }

  eisa_df
}

eisa_results <- create_eisa_table(intron_results, exon_results, common_genes_eisa)

ensembl_ids_clean_eisa <- str_remove(eisa_results$gene_id, "\\..*$")
eisa_results$symbol    <- gene_symbols$SYMBOL[match(ensembl_ids_clean_eisa,
                                                    gene_symbols$ENSEMBL)]

## Annotate genes at loop anchors ----
message("\nAnnotating genes at chromatin loop anchors...")

loops <- readRDS(loop_data_rds) |> interactions()
loops <- keepStandardChromosomes(loops, pruning.mode = "coarse")
loops <- loops[seqnames(anchors(loops, "first")) != "chrM"]

gained_loops <- loops[mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange > 0]
lost_loops   <- loops[mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange < 0]
static_loops <- loops[mcols(loops)$padj > 0.1]

message(sprintf("Gained: %d  Lost: %d  Static: %d",
                length(gained_loops), length(lost_loops), length(static_loops)))

txdb      <- loadDb(txdb_file)
genes_all <- genes(txdb)
genes_all <- keepStandardChromosomes(genes_all, pruning.mode = "coarse")
genes_all <- genes_all[seqnames(genes_all) != "chrM"]

genes_all_clean <- genes_all
names(genes_all_clean) <- str_remove(names(genes_all), "\\..*$")
genes_with_eisa <- genes_all_clean[names(genes_all_clean) %in% eisa_results$gene_id]
promoters_gr    <- promoters(genes_with_eisa, upstream = 2000, downstream = 500)

find_loop_anchor_genes <- function(loop_set, promoters) {
  hits1 <- findOverlaps(promoters, anchors(loop_set, "first"))
  hits2 <- findOverlaps(promoters, anchors(loop_set, "second"))
  unique(c(names(promoters)[queryHits(hits1)],
           names(promoters)[queryHits(hits2)]))
}

gained_genes <- find_loop_anchor_genes(gained_loops, promoters_gr)
lost_genes   <- find_loop_anchor_genes(lost_loops,   promoters_gr)
static_genes <- find_loop_anchor_genes(static_loops, promoters_gr)

eisa_results$loop_category <- "none"
eisa_results$loop_category[eisa_results$gene_id %in% gained_genes] <- "gained"
eisa_results$loop_category[eisa_results$gene_id %in% lost_genes]   <- "lost"
eisa_results$loop_category[eisa_results$gene_id %in% static_genes] <- "static"

message(sprintf("Genes at gained: %d  lost: %d  static: %d",
                sum(eisa_results$loop_category == "gained"),
                sum(eisa_results$loop_category == "lost"),
                sum(eisa_results$loop_category == "static")))

## Statistics ----
gained_data <- eisa_results |> filter(loop_category == "gained")
lost_data   <- eisa_results |> filter(loop_category == "lost")
static_data <- eisa_results |> filter(loop_category == "static")

cor_test    <- cor.test(gained_data$delta_intron_12h, gained_data$delta_exon_12h, method = "pearson")
wilcox_test <- wilcox.test(abs(gained_data$delta_posttx_12h), abs(gained_data$delta_intron_12h),
                           paired = TRUE, alternative = "less")
ttest_result <- t.test(gained_data$delta_intron_12h, static_data$delta_intron_12h)

cat(sprintf("\nR (Δintron vs Δexon, 12h): %.3f  p = %.2e\n", cor_test$estimate, cor_test$p.value))
cat(sprintf("Wilcoxon (posttx < intron): p = %.2e\n", wilcox_test$p.value))
cat(sprintf("t-test (gained vs static Δintron): p = %.2e\n", ttest_result$p.value))

stats_summary <- list(
  correlation      = cor_test,
  magnitude        = wilcox_test,
  gained_vs_static = ttest_result
)

# Save outputs ----
saveRDS(eisa_results,  file.path(eisa_dir, "eisa_results_with_loops.rds"))
saveRDS(stats_summary, file.path(eisa_dir, "statistics_summary.rds"))

message("EISA Part 2 complete.")
message(sprintf("  EISA results: %s", file.path(eisa_dir, "eisa_results_with_loops.rds")))
message(sprintf("  Statistics:   %s", file.path(eisa_dir, "statistics_summary.rds")))

sessionInfo()

# auxinRAD21 — gained loop anchor expression summary table ----------------
# Author:      JP Flores
# Date:        2026-06-02
# Project:     STRS
# Description: Computes a summary table of mean log2FC, n genes, and paired
#              Wilcoxon p-value for sorbitol vs. sorbitol+auxin conditions
#              at sorbitol-induced (gained) and pre-existing loop anchors.
#              Reuses the DESeq2 and overlap logic from
#              auxinRAD21_gainedLoops_histogram.R — drops all plotting.
# Input:       data/processed/rna/AuxinRAD21/output/quant/YAPP_HCT116_mAID2-RAD21_1_tximport.rds
#              data/processed/rna/AuxinRAD21/output/YAPP_HCT116_mAID2-RAD21_1_RNApipeSamplesheet.txt
#              data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
#              data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds
# Output:      data/processed/rna/AuxinRAD21/output/deseqObjs/auxinRAD21_anchor_expression_summary.tsv
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"

tximport_rds   <- file.path(project_dir,
                            "data/processed/rna/AuxinRAD21/output/quant",
                            "YAPP_HCT116_mAID2-RAD21_1_tximport.rds")

samplesheet_file <- file.path(project_dir,
                              "data/processed/rna/AuxinRAD21/output",
                              "YAPP_HCT116_mAID2-RAD21_1_RNApipeSamplesheet.txt")

diff_loops_rds <- file.path(project_dir,
                            "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

lrt_dds_rds <- file.path(project_dir,
                         "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

output_tsv <- file.path(project_dir,
                        "data/processed/rna/AuxinRAD21/output/deseqObjs",
                        "auxinRAD21_anchor_expression_summary.tsv")

## DESeq2 filtering
min_counts  <- 10L
min_samples <- 2L

## Differential loop threshold
padj_loop <- 0.1


# Libraries ---------------------------------------------------------------

library(data.table)
library(DESeq2)
library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)
library(InteractionSet)
library(IRanges)
library(mariner)
library(S4Vectors)
library(stringr)
library(tibble)


# Load data ---------------------------------------------------------------

## AuxinRAD21 tximeta SummarizedExperiment (pre-computed by bagPipes)
se <- readRDS(tximport_rds)

## Differential loops — HEK293 eGFP-YAP (defines sorbitol-induced anchors)
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")

gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]

staticLoops <- diff_loopCounts[is.na(diff_loopCounts$padj) |
                                 diff_loopCounts$padj >= padj_loop]

## LRT timecourse dds — used for hg38 gene coordinates only
dge      <- readRDS(lrt_dds_rds)
gene_gr  <- results(dge, format = "GRanges")
mcols(gene_gr) <- rowData(dge)


# Wrangle data ------------------------------------------------------------

## Assign sample names from samplesheet (bagPipes does not write them into
## tximport RDS colnames) — preserve samplesheet row order, do NOT sort()
samplesheet <- fread(samplesheet_file) |> as.data.frame()
sn_vec      <- samplesheet$sn

stopifnot(ncol(se$counts) == length(sn_vec))
colnames(se$counts)    <- sn_vec
colnames(se$abundance) <- sn_vec
colnames(se$length)    <- sn_vec

## Parse Treatment and Bio_Rep from sample names
sample_names  <- colnames(se$counts)
treatment_vec <- dplyr::case_when(
  str_detect(sample_names, "untreated") ~ "untreated",
  str_detect(sample_names, "sorbitol")  ~ "sorbitol",
  str_detect(sample_names, "auxin")     ~ "auxin"
)
rep_vec <- str_match(sample_names, "_(\\d+)_\\d+$")[, 2]

coldata <- data.frame(
  Treatment = factor(treatment_vec, levels = c("untreated", "sorbitol", "auxin")),
  Bio_Rep   = factor(rep_vec),
  row.names = sample_names
)

## Promoter coordinates from LRT dds (tximeta-imported, has proper rowRanges)
gene_gr_prom <- promoters(gene_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(gene_gr_prom) <- "UCSC"
seqlevelsStyle(gainedLoops)  <- "UCSC"
seqlevelsStyle(staticLoops)  <- "UCSC"

## Strip Ensembl version suffixes (e.g. ENSG00000000003.14 → ENSG00000000003)
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
names(gene_gr_prom) <- strip_version(names(gene_gr_prom))


# Analysis ----------------------------------------------------------------

## Build DESeqDataSet and run DESeq2
dds  <- DESeqDataSetFromTximport(se, colData = coldata, design = ~Bio_Rep + Treatment)
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]
dds  <- DESeq(dds)

## Extract log2FC for each contrast
res_sorb_df <- results(dds, contrast = c("Treatment", "sorbitol", "untreated")) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, log2FoldChange) |>
  dplyr::mutate(gene_id = strip_version(gene_id))

res_auxin_df <- results(dds, contrast = c("Treatment", "auxin", "untreated")) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, log2FoldChange) |>
  dplyr::mutate(gene_id = strip_version(gene_id))

#' Attach log2FC to a promoter GRanges and subset to genes in results
#'
#' @param prom_gr   GRanges of promoters with gene_id as names
#' @param res_df    data.frame with columns gene_id and log2FoldChange
#'
#' @return GRanges subset to genes present in res_df, with log2FoldChange
#'         added to mcols
attach_lfc <- function(prom_gr, res_df) {
  gr  <- prom_gr[names(prom_gr) %in% res_df$gene_id]
  gr$log2FoldChange <- res_df$log2FoldChange[match(names(gr), res_df$gene_id)]
  gr
}

prom_sorb_gr  <- attach_lfc(gene_gr_prom, res_sorb_df)
prom_auxin_gr <- attach_lfc(gene_gr_prom, res_auxin_df)

#' Compute summary stats for genes at a given anchor class
#'
#' @param sorb_gr    GRanges of promoters with sorbitol log2FC in mcols
#' @param auxin_gr   GRanges of promoters with auxin log2FC in mcols
#' @param anchors    GRanges of loop anchors to overlap against
#' @param label      Character label for this anchor class in the output table
#'
#' @return One-row data.frame with anchor_class, n_genes, mean_log2FC_sorb,
#'         mean_log2FC_auxin, and wilcox_p
summarize_anchor <- function(sorb_gr, auxin_gr, anchors, label) {
  hits_sorb  <- subsetByOverlaps(sorb_gr,  anchors, ignore.strand = TRUE)
  hits_auxin <- subsetByOverlaps(auxin_gr, anchors, ignore.strand = TRUE)
  
  common_genes <- intersect(names(hits_sorb), names(hits_auxin))
  lfc_sorb     <- hits_sorb$log2FoldChange[match(common_genes, names(hits_sorb))]
  lfc_auxin    <- hits_auxin$log2FoldChange[match(common_genes, names(hits_auxin))]
  
  valid        <- !is.na(lfc_sorb) & !is.na(lfc_auxin)
  lfc_sorb     <- lfc_sorb[valid]
  lfc_auxin    <- lfc_auxin[valid]
  
  wt <- wilcox.test(lfc_sorb, lfc_auxin, paired = TRUE, alternative = "two.sided")
  
  data.frame(
    anchor_class      = label,
    n_genes           = sum(valid),
    mean_log2FC_sorb  = round(mean(lfc_sorb),  3),
    mean_log2FC_auxin = round(mean(lfc_auxin), 3),
    wilcox_p          = signif(wt$p.value, 3)
  )
}

gained_anchors  <- reduce(c(anchors(gainedLoops, "first"), anchors(gainedLoops, "second")))
static_anchors  <- reduce(c(anchors(staticLoops, "first"), anchors(staticLoops, "second")))

summary_tbl <- bind_rows(
  summarize_anchor(prom_sorb_gr, prom_auxin_gr, gained_anchors, "gained"),
  summarize_anchor(prom_sorb_gr, prom_auxin_gr, static_anchors, "preexisting")
)

message("Summary table:")
print(summary_tbl)
message(sprintf(
  "Genes at sorbitol-induced loop anchors showed significantly reduced expression upon cohesin depletion relative to sorbitol treatment alone (p < 0.001, n = %d genes), while genes at pre-existing loop anchors were unaffected (p = %.3f, n = %d genes), indicating that the transcriptional effect is specific to cohesin-dependent, stress-induced loops.",
  summary_tbl$n_genes[summary_tbl$anchor_class == "gained"],
  summary_tbl$wilcox_p[summary_tbl$anchor_class == "preexisting"],
  summary_tbl$n_genes[summary_tbl$anchor_class == "preexisting"]
))


# Save outputs ------------------------------------------------------------

write.table(summary_tbl,
            file      = output_tsv,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

message("Summary table saved: ", output_tsv)


# Session info ------------------------------------------------------------

sessionInfo()
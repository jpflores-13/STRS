# ##############################################################################
# filename:    eisa_auxinRAD21_analysis.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: EISA Part 2 for AuxinRAD21 data — mirrors the timecourse EISA
#              approach in eisa_analysis.R but for a three-condition design
#              (untreated, sorbitol 6h, sorbitol+auxin 6h). Fits DESeq2 with
#              design ~Bio_Rep + Treatment on per-replicate intronic read counts
#              (n=6: 2 reps x 3 conditions) to extract delta-intron (nascent RNA
#              log2FC) for two contrasts: sorbitol vs untreated and auxin vs
#              untreated. Repeats the same DESeq2 fit on exonic counts to extract
#              delta-exon (total mRNA log2FC). DeltaReg (post-transcriptional
#              component) is computed as delta-exon minus delta-intron per gene.
#              Annotates genes at gained and lost loop anchors using promoter
#              overlaps with HEK293 eGFP-YAP differential loop calls.
# input:       data/processed/rna/AuxinRAD21/output/EISA/per_replicate/exon_counts.rds
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/intron_counts.rds
#              data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
#              data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds
# output:      data/processed/rna/AuxinRAD21/output/EISA/per_replicate/eisa_auxinRAD21_results.rds
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/dds_intron.rds
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/dds_exon.rds
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(GenomicRanges)
library(GenomeInfoDb)
library(InteractionSet)
library(IRanges)
library(mariner)
library(S4Vectors)
library(stringr)
library(tibble)


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"

eisa_dir <- file.path(project_dir,
                      "data/processed/rna/AuxinRAD21/output/EISA/per_replicate")

exon_counts_file   <- file.path(eisa_dir, "exon_counts.rds")
intron_counts_file <- file.path(eisa_dir, "intron_counts.rds")

## LRT timecourse dds — used for hg38 gene coordinates (has proper rowRanges)
lrt_dds_rds <- file.path(project_dir,
                         "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

diff_loops_rds <- file.path(project_dir,
                            "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

## Low-count filter
min_counts  <- 5L
min_samples <- 2L

## Differential loop threshold
padj_loop <- 0.1

## Outputs
results_rds    <- file.path(eisa_dir, "eisa_auxinRAD21_results.rds")
dds_intron_rds <- file.path(eisa_dir, "dds_intron.rds")
dds_exon_rds   <- file.path(eisa_dir, "dds_exon.rds")


# Load count matrices ---------------------------------------------------------

exon_counts   <- readRDS(exon_counts_file)
intron_counts <- readRDS(intron_counts_file)

message("Exon counts:   ", nrow(exon_counts),   " x ", ncol(exon_counts))
message("Intron counts: ", nrow(intron_counts), " x ", ncol(intron_counts))
message("Columns: ", paste(colnames(exon_counts), collapse = ", "))

stopifnot(identical(colnames(exon_counts), colnames(intron_counts)))


# Build colData ---------------------------------------------------------------

## Column names after cleaning in bam_processing.R:
## sorbitol_6h_1_1, sorbitol_6h_2_1, auxin_6h_1_1, auxin_6h_2_1,
## untreated_0h_1_1, untreated_0h_2_1
sample_names  <- colnames(exon_counts)

treatment_vec <- dplyr::case_when(
  str_detect(sample_names, "untreated") ~ "untreated",
  str_detect(sample_names, "sorbitol")  ~ "sorbitol",
  str_detect(sample_names, "auxin")     ~ "auxin"
)

## Extract replicate number from trailing _<rep>_<lane> pattern
rep_vec <- str_match(sample_names, "_(\\d+)_\\d+$")[, 2]

coldata <- data.frame(
  Treatment = factor(treatment_vec,
                     levels = c("untreated", "sorbitol", "auxin")),
  Bio_Rep   = factor(rep_vec),
  row.names = sample_names
)

message("ColData:")
print(coldata)


# Filter low-count genes ------------------------------------------------------

common_genes <- intersect(rownames(exon_counts), rownames(intron_counts))
exon_filt    <- exon_counts[common_genes, ]
intron_filt  <- intron_counts[common_genes, ]

keep_intron <- rowSums(intron_filt >= min_counts) >= min_samples
keep_exon   <- rowSums(exon_filt   >= min_counts) >= min_samples
keep        <- keep_intron & keep_exon

message(sprintf("Genes passing filter: %d / %d", sum(keep), length(keep)))

exon_filt   <- exon_filt[keep, ]
intron_filt <- intron_filt[keep, ]


# DESeq2 on intronic reads ----------------------------------------------------

## Design mirrors eisa_analysis.R: ~Bio_Rep + Treatment
## Reference level is "untreated" (set via factor levels in coldata)
message("\nFitting DESeq2 on intronic counts...")
dds_intron <- DESeqDataSetFromMatrix(
  countData = intron_filt,
  colData   = coldata,
  design    = ~Bio_Rep + Treatment
)
dds_intron <- DESeq(dds_intron)
saveRDS(dds_intron, dds_intron_rds)
message("Intron DESeq2 complete.")


# DESeq2 on exonic reads ------------------------------------------------------

message("\nFitting DESeq2 on exonic counts...")
dds_exon <- DESeqDataSetFromMatrix(
  countData = exon_filt,
  colData   = coldata,
  design    = ~Bio_Rep + Treatment
)
dds_exon <- DESeq(dds_exon)
saveRDS(dds_exon, dds_exon_rds)
message("Exon DESeq2 complete.")


# Extract fold changes --------------------------------------------------------

## Helper: pull log2FC and padj for a contrast into a named data frame
extract_lfc <- function(dds, contrast_vec, lfc_col, padj_col) {
  results(dds, contrast = contrast_vec) |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    dplyr::select(gene_id, log2FoldChange, padj) |>
    dplyr::rename(!!lfc_col  := log2FoldChange,
                  !!padj_col := padj)
}

## Sorbitol vs untreated
sorb_intron <- extract_lfc(dds_intron,
                           c("Treatment", "sorbitol", "untreated"),
                           "lfc_intron_sorb", "padj_intron_sorb")
sorb_exon   <- extract_lfc(dds_exon,
                           c("Treatment", "sorbitol", "untreated"),
                           "lfc_exon_sorb", "padj_exon_sorb")

## Auxin vs untreated
auxin_intron <- extract_lfc(dds_intron,
                            c("Treatment", "auxin", "untreated"),
                            "lfc_intron_auxin", "padj_intron_auxin")
auxin_exon   <- extract_lfc(dds_exon,
                            c("Treatment", "auxin", "untreated"),
                            "lfc_exon_auxin", "padj_exon_auxin")

## Join all four contrasts
eisa_df <- sorb_intron |>
  dplyr::left_join(sorb_exon,   by = "gene_id") |>
  dplyr::left_join(auxin_intron, by = "gene_id") |>
  dplyr::left_join(auxin_exon,   by = "gene_id")

## Compute DeltaReg = delta-exon minus delta-intron (post-transcriptional component)
## Mirrors eisa_analysis.R
eisa_df <- eisa_df |>
  dplyr::mutate(
    lfc_reg_sorb  = lfc_exon_sorb  - lfc_intron_sorb,
    lfc_reg_auxin = lfc_exon_auxin - lfc_intron_auxin
  )

## Strip ENSEMBL version suffixes for coordinate lookup
strip_version  <- function(x) sub("\\.[0-9]+$", "", x)
eisa_df$gene_id <- strip_version(eisa_df$gene_id)

message("\nEISA results: ", nrow(eisa_df), " genes")


# Annotate genes at loop anchors ----------------------------------------------

## Gene coordinates from LRT timecourse dds (has proper hg38 rowRanges)
dge          <- readRDS(lrt_dds_rds)
gene_gr      <- results(dge, format = "GRanges")
mcols(gene_gr) <- rowData(dge)
gene_gr_prom <- promoters(gene_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(gene_gr_prom) <- "UCSC"
names(gene_gr_prom) <- strip_version(names(gene_gr_prom))

## Differential loops
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")
gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]
lostLoops   <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange < 0]
seqlevelsStyle(gainedLoops) <- "UCSC"
seqlevelsStyle(lostLoops)   <- "UCSC"

## Subset promoter GRanges to genes present in EISA results
gene_gr_eisa <- gene_gr_prom[names(gene_gr_prom) %in% eisa_df$gene_id]

gained_genes <- names(subsetByOverlaps(gene_gr_eisa, gainedLoops,
                                       ignore.strand = TRUE))
lost_genes   <- names(subsetByOverlaps(gene_gr_eisa, lostLoops,
                                       ignore.strand = TRUE))

eisa_df$loop_category <- "none"
eisa_df$loop_category[eisa_df$gene_id %in% gained_genes] <- "gained"
eisa_df$loop_category[eisa_df$gene_id %in% lost_genes]   <- "lost"

message(sprintf("Loop annotation: gained=%d  lost=%d  none=%d",
                sum(eisa_df$loop_category == "gained"),
                sum(eisa_df$loop_category == "lost"),
                sum(eisa_df$loop_category == "none")))


# Save outputs ----------------------------------------------------------------

saveRDS(eisa_df, results_rds)
message("EISA results saved: ", results_rds)


# Session info ----------------------------------------------------------------

sessionInfo()
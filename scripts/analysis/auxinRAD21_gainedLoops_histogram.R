# ##############################################################################
# filename:    auxinRAD21_gainedLoops_histogram.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Overlapping density histograms showing the log2FC distribution
#              (treated/untreated) for genes whose promoters overlap sorbitol-
#              induced (gained) loop anchors. Two distributions are overlaid:
#              sorbitol 6h vs untreated and sorbitol+auxin 6h vs untreated,
#              both from HCT116 mAID2-RAD21 cells. Centered on 0 with a dashed
#              reference line. Gained loop anchors are defined from HEK293
#              eGFP-YAP differential loop analysis.
# input:       AuxinRAD21 tximeta RDS (quant/), differential loop RDS
# output:      plots/auxinRAD21_gainedLoops_histogram.pdf
#              data/processed/rna/AuxinRAD21/output/deseqObjs/auxinRAD21_gg_histogram.rds
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(data.table)
library(DESeq2)
library(dplyr)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
library(InteractionSet)
library(IRanges)
library(mariner)
library(S4Vectors)
library(stringr)
library(tibble)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# Utility scripts -------------------------------------------------------------

source("scripts/utils/ggplot2_pgTheme.R")


# Parameters ------------------------------------------------------------------

project_dir    <- "/work/users/j/p/jpflores/projects/STRS"
plots_dir      <- file.path(project_dir, "plots")

tximport_rds   <- file.path(project_dir,
                            "data/processed/rna/AuxinRAD21/output/quant",
                            "YAPP_HCT116_mAID2-RAD21_1_tximport.rds")

samplesheet_file <- file.path(project_dir,
                              "data/processed/rna/AuxinRAD21/output",
                              "YAPP_HCT116_mAID2-RAD21_1_RNApipeSamplesheet.txt")

diff_loops_rds <- file.path(project_dir,
                            "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

## LRT timecourse dds — used for hg38 gene coordinates (imported with tximeta)
lrt_dds_rds <- file.path(project_dir,
                         "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

## Intermediate ggplot object saved for FigureRAD21.R assembly
gg_rds_out <- file.path(project_dir,
                        "data/processed/rna/AuxinRAD21/output/deseqObjs",
                        "auxinRAD21_gg_histogram.rds")

## DESeq2 filtering
min_counts  <- 10L
min_samples <- 2L

## Differential loop thresholds
padj_loop   <- 0.1

## Histogram bins
n_bins <- 40L

## Colors
color_sorbitol <- "#F8766D"
color_auxin    <- "#619CFF"

## Plot dimensions
plot_width  <- 5
plot_height <- 4

## Output PDF
output_pdf <- file.path(plots_dir, "auxinRAD21_gainedLoops_histogram.pdf")


# Load data -------------------------------------------------------------------

## AuxinRAD21 tximeta SummarizedExperiment (pre-computed by bagPipes)
se <- readRDS(tximport_rds)

## Differential loops — HEK293 eGFP-YAP (defines sorbitol-induced anchors)
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")

gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]


# Wrangle data ----------------------------------------------------------------

## Read samplesheet for sample metadata — same pattern as rnaseqTimecourse_LRT.R
## This is more robust than parsing colnames(se$counts) directly since the
## tximport column names may not follow a predictable substring pattern
samplesheet <- fread(samplesheet_file) |> as.data.frame()
sn_vec      <- samplesheet$sn

## colnames(se$counts) are all NA because bagPipes did not write sample names
## into the tximport RDS. Reconstruct them: bagPipes calls list.files() on the
## quant directories, which returns paths in alphabetical order — so the matrix
## columns correspond to sort(sn_vec).
stopifnot(ncol(se$counts) == length(sn_vec))
sn_ordered             <- sort(sn_vec)
colnames(se$counts)    <- sn_ordered
colnames(se$abundance) <- sn_ordered
colnames(se$length)    <- sn_ordered

message("Assigned column names (alphabetical): ",
        paste(sn_ordered, collapse = ", "))

## Parse Treatment and Bio_Rep from the now-named columns
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


# Analysis --------------------------------------------------------------------

## Build DESeqDataSet — reference level is "untreated" (set via factor levels)
dds  <- DESeqDataSetFromTximport(se, colData = coldata, design = ~Bio_Rep + Treatment)
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]
dds  <- DESeq(dds)

## Extract results as data frames — tximport-based dds has no rowRanges,
## so format = "GRanges" is not available here
res_sorb_df <- results(dds,
                       contrast = c("Treatment", "sorbitol", "untreated")) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, log2FoldChange)

res_auxin_df <- results(dds,
                        contrast = c("Treatment", "auxin", "untreated")) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  dplyr::select(gene_id, log2FoldChange)

## Get gene genomic coordinates from the LRT timecourse dds, which was
## imported with tximeta and has proper hg38 rowRanges — used here only
## for gene body positions, not expression values
dge          <- readRDS(lrt_dds_rds)
gene_gr      <- results(dge, format = "GRanges")
mcols(gene_gr) <- rowData(dge)
gene_gr_prom <- promoters(gene_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(gene_gr_prom) <- "UCSC"

## Strip ENSEMBL version suffixes before joining
## (e.g. ENSG00000000003.14 → ENSG00000000003)
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
names(gene_gr_prom)  <- strip_version(names(gene_gr_prom))
res_sorb_df$gene_id  <- strip_version(res_sorb_df$gene_id)
res_auxin_df$gene_id <- strip_version(res_auxin_df$gene_id)

## Subset to genes present in both coordinate set and AuxinRAD21 results,
## then attach log2FC by gene_id position matching
prom_sorb_gr  <- gene_gr_prom[names(gene_gr_prom) %in% res_sorb_df$gene_id]
prom_auxin_gr <- gene_gr_prom[names(gene_gr_prom) %in% res_auxin_df$gene_id]

prom_sorb_gr$log2FoldChange  <- res_sorb_df$log2FoldChange[
  match(names(prom_sorb_gr),  res_sorb_df$gene_id)]
prom_auxin_gr$log2FoldChange <- res_auxin_df$log2FoldChange[
  match(names(prom_auxin_gr), res_auxin_df$gene_id)]

## seqlevels compatibility for overlap with loop anchors
seqlevelsStyle(gainedLoops) <- "UCSC"

## Genes whose promoters overlap either anchor of a gained loop
hits_sorb  <- subsetByOverlaps(prom_sorb_gr,  gainedLoops, ignore.strand = TRUE)
hits_auxin <- subsetByOverlaps(prom_auxin_gr, gainedLoops, ignore.strand = TRUE)

## Paired Wilcoxon signed-rank test ----------------------------------------
## The same genes appear in both contrasts (same promoters overlapping the same
## gained loop anchors, just with different log2FC values), so the correct test
## is PAIRED — it tests whether the per-gene difference in log2FC between
## conditions is systematically shifted from zero.
common_genes     <- intersect(names(hits_sorb), names(hits_auxin))
lfc_sorb_paired  <- hits_sorb$log2FoldChange[match(common_genes, names(hits_sorb))]
lfc_auxin_paired <- hits_auxin$log2FoldChange[match(common_genes, names(hits_auxin))]

## Drop pairs where either value is NA
valid_pairs      <- !is.na(lfc_sorb_paired) & !is.na(lfc_auxin_paired)
lfc_sorb_paired  <- lfc_sorb_paired[valid_pairs]
lfc_auxin_paired <- lfc_auxin_paired[valid_pairs]
n_pairs          <- sum(valid_pairs)

wt <- wilcox.test(lfc_sorb_paired, lfc_auxin_paired,
                  paired      = TRUE,
                  alternative = "two.sided")

message("Paired Wilcoxon signed-rank test")
message("  n gene pairs : ", n_pairs)
message("  W statistic  : ", wt$statistic)
message("  p-value      : ", wt$p.value)

## Format p-value for plot annotation (parse = TRUE for italic p)
p_label <- if (wt$p.value < 0.001) {
  "italic(p) < 0.001"
} else {
  paste0("italic(p) == ", round(wt$p.value, 3))
}
n_label <- paste0("n = ", n_pairs, " genes")

## Tidy data frame for plotting — one row per gene per condition
hist_df <- dplyr::bind_rows(
  data.frame(log2FC    = hits_sorb$log2FoldChange,
             condition = "Sorbitol"),
  data.frame(log2FC    = hits_auxin$log2FoldChange,
             condition = "Sorbitol + Auxin")
) |>
  dplyr::filter(!is.na(log2FC))

hist_df$condition <- factor(hist_df$condition,
                            levels = c("Sorbitol", "Sorbitol + Auxin"))


# Visualization ---------------------------------------------------------------

## Compute density peaks for direct label placement
## (same approach as Figure 3's create_density_plot)
dens_sorb  <- density(hist_df$log2FC[hist_df$condition == "Sorbitol"],
                      na.rm = TRUE)
dens_auxin <- density(hist_df$log2FC[hist_df$condition == "Sorbitol + Auxin"],
                      na.rm = TRUE)

sorb_label_x  <- dens_sorb$x[which.max(dens_sorb$y)]  + 0.6
sorb_label_y  <- max(dens_sorb$y) * 0.85
auxin_label_x <- dens_auxin$x[which.max(dens_auxin$y)] - 2.8
auxin_label_y <- max(dens_auxin$y) * 0.65

p_hist <- ggplot(hist_df, aes(x = log2FC, fill = condition)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray75", linewidth = 0.3) +
  geom_density(alpha = 0.4, color = NA) +
  annotate("text",
           x        = sorb_label_x,
           y        = sorb_label_y,
           label    = "Sorbitol",
           color    = color_sorbitol,
           fontface = "bold",
           size     = 7 / .pt,
           hjust    = 0) +
  annotate("text",
           x          = auxin_label_x,
           y          = auxin_label_y,
           label      = "Sorbitol +\nAuxin",
           color      = color_auxin,
           fontface   = "bold",
           size       = 7 / .pt,
           hjust      = 0,
           lineheight = 0.8) +
  ## Paired Wilcoxon p-value — upper right where density is sparse
  annotate("text",
           x     = 4.4,
           y     = max(dens_sorb$y) * 0.95,
           label = p_label,
           parse = TRUE,
           size  = 7 / .pt,
           hjust = 1,
           color = "black") +
  annotate("text",
           x     = 4.4,
           y     = max(dens_sorb$y) * 0.83,
           label = n_label,
           size  = 6 / .pt,
           hjust = 1,
           color = "gray40") +
  scale_fill_manual(values = c(
    "Sorbitol"         = color_sorbitol,
    "Sorbitol + Auxin" = color_auxin
  )) +
  coord_cartesian(xlim = c(-4.5, 4.5)) +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  labs(
    x     = expression(log[2]~FC~(treated / untreated)),
    y     = "Density",
    title = "Genes at Sorbitol-Induced Loop Anchors"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 9),
    legend.position = "none",
    axis.text       = element_text(size = 7, color = "black"),
    axis.title      = element_text(size = 8.5),
    axis.line       = element_line(linewidth = 0.3),
    axis.ticks      = element_line(linewidth = 0.3),
    plot.margin     = margin(5.5, 5.5, 5.5, 5.5, "pt")
  )


# Save outputs ----------------------------------------------------------------

## PDF for standalone viewing
ggsave(output_pdf, p_hist,
       width  = plot_width,
       height = plot_height,
       units  = "in")

## ggplot object RDS for assembly in FigureRAD21.R
saveRDS(p_hist, file = gg_rds_out)

message("Histogram saved: ", output_pdf)
message("ggplot RDS saved: ", gg_rds_out)


# Session info ----------------------------------------------------------------

sessionInfo()
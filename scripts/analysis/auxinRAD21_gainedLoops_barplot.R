# ##############################################################################
# filename:    auxinRAD21_gainedLoops_barplot.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Grouped barplot showing mean log2FC (± SEM) for genes at
#              sorbitol-induced (gained) and pre-existing (lost) loop anchors
#              across two conditions: sorbitol 6h vs untreated and
#              sorbitol+auxin 6h vs untreated. Uses HCT116 mAID2-RAD21 RNA-seq.
#              The x-axis uses a two-level grouping label (condition within
#              loop type) produced via facet_wrap(strip.position = "bottom").
#              "Pre-existing loops" are defined as loops lost upon sorbitol
#              treatment (padj < 0.1, log2FC < 0) — i.e. canonical loops
#              present before stress.
# input:       AuxinRAD21 tximeta RDS (quant/), differential loop RDS
# output:      plots/auxinRAD21_gainedLoops_barplot.pdf
#              data/processed/rna/AuxinRAD21/output/deseqObjs/auxinRAD21_gg_barplot.rds
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
library(tidyr)


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
                        "auxinRAD21_gg_barplot.rds")

## DESeq2 filtering
min_counts  <- 10L
min_samples <- 2L

## Differential loop thresholds
padj_loop   <- 0.1

## Colors — consistent with color_sorbitol / color_auxin used elsewhere
color_sorbitol <- "#F8766D"
color_auxin    <- "#619CFF"

## Plot dimensions
plot_width  <- 6
plot_height <- 4

## Output PDF
output_pdf <- file.path(plots_dir, "auxinRAD21_gainedLoops_barplot.pdf")


# Load data -------------------------------------------------------------------

## AuxinRAD21 tximeta SummarizedExperiment
se <- readRDS(tximport_rds)

## Differential loops (HEK293 eGFP-YAP reference)
diff_loopCounts <- readRDS(diff_loops_rds) |> interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts,
                                           pruning.mode = "coarse")

gainedLoops <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange > 0]

## Pre-existing loops = lost loops (present before sorbitol, disrupted after)
lostLoops   <- diff_loopCounts[diff_loopCounts$padj < padj_loop &
                                 diff_loopCounts$log2FoldChange < 0]


# Wrangle data ----------------------------------------------------------------

samplesheet <- fread(samplesheet_file) |> as.data.frame()
sn_vec      <- samplesheet$sn

stopifnot(ncol(se$counts) == length(sn_vec))
sn_ordered             <- sort(sn_vec)
colnames(se$counts)    <- sn_ordered
colnames(se$abundance) <- sn_ordered
colnames(se$length)    <- sn_ordered

message("Assigned column names (alphabetical): ",
        paste(sn_ordered, collapse = ", "))

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

dds  <- DESeqDataSetFromTximport(se, colData = coldata, design = ~Bio_Rep + Treatment)
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]
dds  <- DESeq(dds)

## Results as data frames — no format = "GRanges" for tximport-based dds
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

## Gene coordinates from LRT timecourse dds (tximeta — has proper hg38 rowRanges)
dge          <- readRDS(lrt_dds_rds)
gene_gr      <- results(dge, format = "GRanges")
mcols(gene_gr) <- rowData(dge)
gene_gr_prom <- promoters(gene_gr) |>
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(gene_gr_prom) <- "UCSC"

## Strip ENSEMBL version suffixes before joining
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
names(gene_gr_prom)  <- strip_version(names(gene_gr_prom))
res_sorb_df$gene_id  <- strip_version(res_sorb_df$gene_id)
res_auxin_df$gene_id <- strip_version(res_auxin_df$gene_id)

## Subset to genes present in both datasets, then attach log2FC by gene_id
prom_sorb_gr  <- gene_gr_prom[names(gene_gr_prom) %in% res_sorb_df$gene_id]
prom_auxin_gr <- gene_gr_prom[names(gene_gr_prom) %in% res_auxin_df$gene_id]

prom_sorb_gr$log2FoldChange  <- res_sorb_df$log2FoldChange[
  match(names(prom_sorb_gr),  res_sorb_df$gene_id)]
prom_auxin_gr$log2FoldChange <- res_auxin_df$log2FoldChange[
  match(names(prom_auxin_gr), res_auxin_df$gene_id)]

seqlevelsStyle(gainedLoops) <- "UCSC"
seqlevelsStyle(lostLoops)   <- "UCSC"

## Genes at sorbitol-induced loop anchors
sorb_at_gained  <- subsetByOverlaps(prom_sorb_gr,  gainedLoops, ignore.strand = TRUE)
auxin_at_gained <- subsetByOverlaps(prom_auxin_gr, gainedLoops, ignore.strand = TRUE)

## Genes at pre-existing (lost) loop anchors
sorb_at_lost  <- subsetByOverlaps(prom_sorb_gr,  lostLoops, ignore.strand = TRUE)
auxin_at_lost <- subsetByOverlaps(prom_auxin_gr, lostLoops, ignore.strand = TRUE)

## Paired Wilcoxon signed-rank tests ------------------------------------------
## Same gene appears in both contrasts → paired test is appropriate

format_pval <- function(p) {
  if (p < 0.001) "p < 0.001"
  else paste0("p = ", round(p, 3))
}

## Gained loops: sorbitol vs sorbitol+auxin
common_gained   <- intersect(names(sorb_at_gained), names(auxin_at_gained))
lfc_sorb_g      <- sorb_at_gained$log2FoldChange[
  match(common_gained, names(sorb_at_gained))]
lfc_auxin_g     <- auxin_at_gained$log2FoldChange[
  match(common_gained, names(auxin_at_gained))]
valid_g         <- !is.na(lfc_sorb_g) & !is.na(lfc_auxin_g)
wt_gained       <- wilcox.test(lfc_sorb_g[valid_g], lfc_auxin_g[valid_g],
                               paired = TRUE, alternative = "two.sided")
n_gained        <- sum(valid_g)

message("Sorbitol-induced loops — paired Wilcoxon:")
message("  n = ", n_gained, " genes | W = ", wt_gained$statistic,
        " | p = ", wt_gained$p.value)

## Pre-existing loops: sorbitol vs sorbitol+auxin
common_lost     <- intersect(names(sorb_at_lost), names(auxin_at_lost))
lfc_sorb_l      <- sorb_at_lost$log2FoldChange[
  match(common_lost, names(sorb_at_lost))]
lfc_auxin_l     <- auxin_at_lost$log2FoldChange[
  match(common_lost, names(auxin_at_lost))]
valid_l         <- !is.na(lfc_sorb_l) & !is.na(lfc_auxin_l)
wt_lost         <- wilcox.test(lfc_sorb_l[valid_l], lfc_auxin_l[valid_l],
                               paired = TRUE, alternative = "two.sided")
n_lost          <- sum(valid_l)

message("Pre-existing loops — paired Wilcoxon:")
message("  n = ", n_lost, " genes | W = ", wt_lost$statistic,
        " | p = ", wt_lost$p.value)

## Build tidy data frame with all four groups
bar_raw_df <- dplyr::bind_rows(
  data.frame(
    log2FC     = sorb_at_gained$log2FoldChange,
    condition  = "Sorbitol",
    loop_group = "Sorbitol-Induced Loops"
  ),
  data.frame(
    log2FC     = auxin_at_gained$log2FoldChange,
    condition  = "Sorbitol + Auxin",
    loop_group = "Sorbitol-Induced Loops"
  ),
  data.frame(
    log2FC     = sorb_at_lost$log2FoldChange,
    condition  = "Sorbitol",
    loop_group = "Pre-Existing Loops"
  ),
  data.frame(
    log2FC     = auxin_at_lost$log2FoldChange,
    condition  = "Sorbitol + Auxin",
    loop_group = "Pre-Existing Loops"
  )
) |>
  dplyr::filter(!is.na(log2FC))

bar_raw_df$condition  <- factor(bar_raw_df$condition,
                                levels = c("Sorbitol", "Sorbitol + Auxin"))
bar_raw_df$loop_group <- factor(bar_raw_df$loop_group,
                                levels = c("Sorbitol-Induced Loops",
                                           "Pre-Existing Loops"))

## Summary statistics: mean ± SEM per group
## Used as bar heights and error bar extents
bar_summary <- bar_raw_df |>
  dplyr::group_by(loop_group, condition) |>
  dplyr::summarise(
    mean_lfc = mean(log2FC, na.rm = TRUE),
    se_lfc   = sd(log2FC, na.rm = TRUE) / sqrt(dplyr::n()),
    n        = dplyr::n(),
    .groups  = "drop"
  )

## Significance annotation data frame — one row per facet (loop_group)
## y positions computed relative to the tallest bar top in each facet
bracket_offset <- diff(range(c(bar_summary$mean_lfc - bar_summary$se_lfc,
                               bar_summary$mean_lfc + bar_summary$se_lfc),
                             na.rm = TRUE)) * 0.06

sig_df <- bar_summary |>
  dplyr::group_by(loop_group) |>
  dplyr::summarise(
    y_bracket = max(mean_lfc + se_lfc, na.rm = TRUE) + bracket_offset,
    .groups   = "drop"
  ) |>
  dplyr::mutate(
    pval_label = c(format_pval(wt_gained$p.value),
                   format_pval(wt_lost$p.value)),
    n_label    = c(paste0("n = ", n_gained, " genes"),
                   paste0("n = ", n_lost,   " genes")),
    y_text     = y_bracket + bracket_offset * 0.8,
    y_n        = y_bracket + bracket_offset * 2.2
  )


# Visualization ---------------------------------------------------------------

p_bar <- ggplot(bar_summary,
                aes(x = condition, y = mean_lfc, fill = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.3) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_lfc - se_lfc,
        ymax = mean_lfc + se_lfc),
    width     = 0.22,
    linewidth = 0.4
  ) +
  ## Significance bracket — horizontal segment spanning the two bars
  geom_segment(data = sig_df,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE,
               linewidth = 0.3, color = "black") +
  ## Short tick marks at each end of the bracket
  geom_segment(data = sig_df,
               aes(x = 1, xend = 1,
                   y = y_bracket, yend = y_bracket - bracket_offset * 0.4),
               inherit.aes = FALSE, linewidth = 0.3, color = "black") +
  geom_segment(data = sig_df,
               aes(x = 2, xend = 2,
                   y = y_bracket, yend = y_bracket - bracket_offset * 0.4),
               inherit.aes = FALSE, linewidth = 0.3, color = "black") +
  ## P-value text above bracket midpoint
  geom_text(data = sig_df,
            aes(x = 1.5, y = y_text, label = pval_label),
            inherit.aes = FALSE, size = 6 / .pt, color = "black") +
  ## n genes in smaller gray text
  geom_text(data = sig_df,
            aes(x = 1.5, y = y_n, label = n_label),
            inherit.aes = FALSE, size = 5 / .pt, color = "gray40") +
  scale_fill_manual(values = c(
    "Sorbitol"         = color_sorbitol,
    "Sorbitol + Auxin" = color_auxin
  )) +
  ## Extra top expansion to accommodate the bracket and text
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.30))) +
  ## strip.position = "bottom" + strip.placement = "outside" puts the group
  ## label ("Sorbitol-Induced Loops", "Pre-Existing Loops") below the x-axis
  ## tick labels, with a thin border acting as the visual bracket
  facet_wrap(~loop_group, strip.position = "bottom") +
  labs(
    x = NULL,
    y = expression(Mean~log[2]~FC~(treated / untreated))
  ) +
  theme_classic() +
  theme(
    strip.placement    = "outside",
    strip.background   = element_rect(fill = NA, color = "black",
                                      linewidth = 0.3),
    strip.text.x       = element_text(size = 7, face = "bold"),
    legend.position    = "none",
    axis.text          = element_text(size = 7, color = "black"),
    axis.title         = element_text(size = 8.5),
    axis.line          = element_line(linewidth = 0.3),
    axis.ticks         = element_line(linewidth = 0.3),
    plot.margin        = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    panel.spacing      = unit(1.5, "lines")
  )


# Save outputs ----------------------------------------------------------------

ggsave(output_pdf, p_bar,
       width  = plot_width,
       height = plot_height,
       units  = "in")

saveRDS(p_bar, file = gg_rds_out)

message("Barplot saved: ", output_pdf)
message("ggplot RDS saved: ", gg_rds_out)


# Session info ----------------------------------------------------------------

sessionInfo()
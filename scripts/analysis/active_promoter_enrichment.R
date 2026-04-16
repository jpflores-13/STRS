# Active Promoter Enrichment at Loop Anchors (Reviewer Response) ----------
# Author:      JP Flores
# Date:        2026-04-07
# Project:     STRS
# Description: Responds to reviewer comment requesting that Fig. 4B be redone
#              using "active" rather than "annotated" promoters. Active genes
#              are defined using a data-driven expression threshold derived from
#              baseline (0h control) RNA-seq normalized counts. Produces:
#                (1) A density plot of log2(normCount + 1) at 0h — inspect the
#                    left shoulder of the distribution to set expr_threshold.
#                    Note: the LRT DDS is pre-filtered (>= 10 counts in >= 4
#                    samples), so the distribution is unimodal; threshold of 5
#                    captures the expressed portion of the transcriptome.
#                (2) A revised stacked bar plot (Fig. 4B-style) breaking loop
#                    anchor promoter overlap into Active / Inactive / No
#                    Promoter, stratified by loop class (Gained, Lost, Static).
# Input:       - data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
#                  GInteractions object with DESeq2 mcols() (same as Fig. 4)
#              - data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds
#                  DESeqDataSet; rowData contains gene_id (Ensembl) and symbol;
#                  colData contains a "Time" column (e.g. "0h", "1h", ...)
# Output:      - plots/active_promoter_enrichment_histogram.pdf
#              - plots/active_promoter_enrichment.pdf
#              - data/processed/rna/timecourse/output/active_genes_0h.rds
#                  character vector of gene symbols classified as active at 0h
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir    <- "/work/users/j/p/jpflores/projects/STRS"
diff_loops_rds <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
dds_rds        <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
output_density <- "plots/active_promoter_enrichment_histogram.pdf"
output_barplot <- "plots/active_promoter_enrichment.pdf"
output_genes   <- "data/processed/rna/timecourse/output/active_genes_0h.rds"

## Which colData "Time" level is the untreated baseline
baseline_timepoint <- "0h"

## log2(normCount + 1) cutoff to call a gene "active" at baseline.
## The LRT DDS is pre-filtered (>= 10 counts in >= 4 samples), so the
## distribution is unimodal. A threshold of 5 sits on the left shoulder
## and captures genes expressed at a meaningful level (~96% of the
## pre-filtered set). Adjust if needed after inspecting the density plot.
expr_threshold <- 5

## Promoter window — mirrors Figure4.R exactly
promoter_upstream   <- 2000
promoter_downstream <- 500

## Loop classification thresholds — mirrors Figure4.R exactly
padj_cutoff <- 0.05
lfc_cutoff  <- 1

## Output figure dimensions (inches)
pdf_width  <- 4
pdf_height <- 4.5


# Libraries ---------------------------------------------------------------

library(DESeq2)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(glue)
library(scales)


# Utility scripts ---------------------------------------------------------

source(file.path(project_dir, "scripts/utils/ggplot2_pgTheme.R"))


# Load data ---------------------------------------------------------------

## GInteractions object — identical to what Figure4.R Panel B uses
diff_loopCounts <- readRDS(file.path(project_dir, diff_loops_rds)) |>
  interactions()
diff_loopCounts <- keepStandardChromosomes(diff_loopCounts, pruning.mode = "coarse")

## LRT DESeqDataSet — same object used in Figure5.R
## rowData(dds) has columns: gene_id (Ensembl), symbol
## colData(dds) has column:  Time (character: "0h", "1h", "3h", ...)
dds <- readRDS(file.path(project_dir, dds_rds))


# Wrangle data ------------------------------------------------------------

## --- Step 1: Annotate loop types (mirrors Figure4.R exactly) ------------

mcols(diff_loopCounts)$loop_type <- dplyr::case_when(
  mcols(diff_loopCounts)$padj < padj_cutoff &
    mcols(diff_loopCounts)$log2FoldChange >  lfc_cutoff ~ "Gained",
  mcols(diff_loopCounts)$padj < padj_cutoff &
    mcols(diff_loopCounts)$log2FoldChange < -lfc_cutoff ~ "Lost",
  mcols(diff_loopCounts)$padj > padj_cutoff              ~ "Static",
  TRUE                                                    ~ "Other"
)

## --- Step 2: Build promoter GRanges from TxDb ---------------------------

txdb        <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
seqlevels(diff_loopCounts) <- seqlevels(txdb)
seqinfo(diff_loopCounts)   <- seqinfo(txdb)

promoter_gr <- genes(txdb) |>
  promoters(upstream = promoter_upstream, downstream = promoter_downstream)

## --- Step 3: Define active genes from RNA-seq 0h baseline ---------------

## Size-factor-normalized counts (DESeq2 corrects for library depth)
norm_counts <- counts(dds, normalized = TRUE)

## Identify samples corresponding to the untreated control time point
baseline_cols <- colData(dds) |>
  as.data.frame() |>
  dplyr::filter(Time == baseline_timepoint) |>
  rownames()

stopifnot(
  "No baseline samples found — check that colData(dds) has a 'Time' column
   whose values match baseline_timepoint" = length(baseline_cols) > 0
)

## Average normalized counts across 0h replicates, then log2(x + 1)
## Pseudocount of 1 prevents log2(0) = -Inf for unexpressed genes
baseline_log2expr <- norm_counts[, baseline_cols] |>
  rowMeans() |>
  (\(x) log2(x + 1))()

## Pull the gene_id (Ensembl) <-> symbol table from rowData, as in Figure5.R
gene_info <- rowData(dds) |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_rowname") |>
  dplyr::select(gene_id, symbol) |>
  dplyr::distinct(gene_id, .keep_all = TRUE)

## Map Ensembl rownames -> Entrez IDs to match promoter_gr$gene_id
## TxDb uses Entrez IDs as its gene identifier
ensembl_ids <- names(baseline_log2expr)

ensembl_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = ensembl_ids,
  keytype = "ENSEMBL",
  columns = "ENTREZID"
) |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

## Named vector: Entrez ID -> log2 expression at baseline
expr_by_entrez <- baseline_log2expr[ensembl_to_entrez$ENSEMBL]
names(expr_by_entrez) <- ensembl_to_entrez$ENTREZID
expr_by_entrez <- expr_by_entrez[!is.na(names(expr_by_entrez))]

## Active = expressed at or above threshold at 0h
active_entrez <- names(expr_by_entrez)[expr_by_entrez >= expr_threshold]

## Store active gene symbols for output (mirrors Figure5.R style)
active_symbols <- gene_info |>
  dplyr::filter(
    gene_id %in% ensembl_to_entrez$ENSEMBL[
      ensembl_to_entrez$ENTREZID %in% active_entrez
    ]
  ) |>
  dplyr::pull(symbol)

message(glue(
  "Active genes at {baseline_timepoint} ",
  "(log2(normCount+1) >= {expr_threshold}): ",
  "{length(active_entrez)} / {length(expr_by_entrez)} ",
  "({round(100 * length(active_entrez) / length(expr_by_entrez), 1)}%)"
))

## Stamp active/inactive status onto promoter_gr so it travels with overlaps
promoter_gr$is_active <- promoter_gr$gene_id %in% active_entrez

## --- Step 4: Classify each loop anchor ----------------------------------
## Three-way: Active Promoter / Inactive Promoter / No Promoter
## An anchor is "Active" if it overlaps at least one active promoter.
## An anchor is "Inactive" if it overlaps a promoter, but none are active.

classify_anchors <- function(loop_subset, promoter_gr) {
  if (length(loop_subset) == 0) {
    return(c(active = 0L, inactive = 0L, no_promoter = 0L, total = 0L))
  }
  
  all_anchors <- c(
    anchors(loop_subset, "first"),
    anchors(loop_subset, "second")
  )
  
  hits <- findOverlaps(all_anchors, promoter_gr)
  
  has_any_promoter <- seq_along(all_anchors) %in% queryHits(hits)
  
  ## For each anchor, check whether any overlapping promoter is active
  has_active <- logical(length(all_anchors))
  active_flags <- promoter_gr$is_active[subjectHits(hits)]
  for (q in unique(queryHits(hits))) {
    has_active[q] <- any(active_flags[queryHits(hits) == q])
  }
  
  c(
    active      = sum(has_active),
    inactive    = sum(has_any_promoter & !has_active),
    no_promoter = sum(!has_any_promoter),
    total       = length(all_anchors)
  )
}

loop_classes <- c("Gained", "Lost", "Static")

anchor_counts <- purrr::map_dfr(loop_classes, function(cls) {
  subset_gi <- diff_loopCounts[mcols(diff_loopCounts)$loop_type == cls]
  cts       <- classify_anchors(subset_gi, promoter_gr)
  tibble::tibble(
    loop_type   = cls,
    active      = cts["active"],
    inactive    = cts["inactive"],
    no_promoter = cts["no_promoter"],
    total       = cts["total"]
  )
})

## Long format for ggplot stacked bar
bar_df <- anchor_counts |>
  tidyr::pivot_longer(
    cols      = c(active, inactive, no_promoter),
    names_to  = "category",
    values_to = "n"
  ) |>
  dplyr::mutate(
    pct      = 100 * n / total,
    category = dplyr::recode(category,
                             active      = "Active Promoter",
                             inactive    = "Inactive Promoter",
                             no_promoter = "No Promoter"
    ),
    category  = factor(category,
                       levels = c("No Promoter",
                                  "Inactive Promoter",
                                  "Active Promoter")),
    loop_type = factor(loop_type, levels = loop_classes)
  )

## Tidy data frame for the density plot
expr_density_df <- tibble::tibble(log2_expr = baseline_log2expr)


# Visualization -----------------------------------------------------------

color_upregulated   <- "#F8766D"
color_downregulated <- "#619CFF"
color_static        <- "#999999"

## --- Plot 1: Expression density -----------------------------------------

p_density <- ggplot(expr_density_df, aes(x = log2_expr)) +
  geom_density(
    fill      = "grey80",
    color     = "grey40",
    linewidth = 0.6,
    alpha     = 0.8
  ) +
  geom_vline(
    xintercept = expr_threshold,
    color      = "firebrick",
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  annotate(
    "text",
    x     = expr_threshold + 0.1,
    y     = Inf,
    label = glue("threshold = {expr_threshold}"),
    color = "firebrick",
    hjust = 0,
    vjust = 1.5,
    size  = 3
  ) +
  labs(
    x        = expression(log[2](normalized~count + 1)~"at 0h"),
    y        = "Density",
    title    = "Baseline expression distribution",
    subtitle = glue(
      "Active >= threshold: {length(active_entrez)} ",
      "({round(100 * length(active_entrez) / length(expr_by_entrez), 1)}%)"
    )
  )

## --- Plot 2: Revised Fig. 4B — styled to match Figure4.R Panel B -------
## Lighter hues of the loop-class colors for each category:
##   Active Promoter   -> light versions of Gained/Lost/Static
##   Inactive Promoter -> even lighter / more pastel
##   No Promoter       -> light grey (same as Fig4.R non_promoter_fraction)

## Flat colors — one per category, consistent across all loop types
fill_values <- c(
  "Active Promoter.Gained"    = "#F1948A",
  "Active Promoter.Lost"      = "#F1948A",
  "Active Promoter.Static"    = "#F1948A",
  "Inactive Promoter.Gained"  = "#AAB7B8",
  "Inactive Promoter.Lost"    = "#AAB7B8",
  "Inactive Promoter.Static"  = "#AAB7B8",
  "No Promoter.Gained"        = "#EAECEE",
  "No Promoter.Lost"          = "#EAECEE",
  "No Promoter.Static"        = "#EAECEE"
)

## Label data: show % for Active and Inactive Promoter segments
label_df <- bar_df |>
  dplyr::filter(category != "No Promoter") |>
  dplyr::mutate(label = paste0(round(pct, 1)))

## x-axis label colors: Gained = red, Lost = blue, Static = grey
axis_label_colors <- c(color_upregulated, color_downregulated, color_static)

p_bar <- ggplot(
  bar_df,
  aes(
    x    = loop_type,
    y    = pct,
    fill = interaction(category, loop_type)
  )
) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(
    data     = label_df,
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size     = 2.5
  ) +
  scale_fill_manual(
    values = fill_values,
    ## Restrict legend to one representative key per category
    breaks = c("Active Promoter.Gained",
               "Inactive Promoter.Gained",
               "No Promoter.Gained"),
    labels = c("Active Promoter", "Inactive Promoter", "No Promoter")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  ## No parentheses around n =
  scale_x_discrete(
    labels = setNames(
      glue("{anchor_counts$loop_type}\nn = {anchor_counts$total}"),
      anchor_counts$loop_type
    )
  ) +
  labs(
    x = "Loop Anchor Type",
    y = "Promoter Overlap (%)"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(fill = c("#F1948A", "#AAB7B8", "#EAECEE")),
      title        = NULL,
      nrow         = 1
    )
  ) +
  theme_minimal() +
  theme(
    ## Colored x-axis labels: Gained = red, Lost = blue, Static = grey
    axis.text.x        = element_text(size = 8.5, color = axis_label_colors),
    axis.text.y        = element_text(size = 8.5),
    axis.title         = element_text(size = 9.5),
    legend.position    = "top",
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(0, 0, 0, 0, "pt")
  )


# Save outputs ------------------------------------------------------------

## Active gene symbol list — consistent with Figure5.R's gene_info usage
saveRDS(
  active_symbols,
  file = file.path(project_dir, output_genes)
)

## Plot 1 — expression density
pdf(
  file   = file.path(project_dir, output_density),
  width  = pdf_width,
  height = pdf_height
)
print(p_density)
dev.off()

## Plot 2 — revised Fig. 4B stacked bar plot
pdf(
  file   = file.path(project_dir, output_barplot),
  width  = pdf_width,
  height = pdf_height
)
print(p_bar)
dev.off()

message("Done. Outputs written to: ", file.path(project_dir, "plots/"))


# Session info ------------------------------------------------------------

sessionInfo()
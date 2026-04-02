# ##############################################################################
# filename:    GOanalysis_timecourse.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: GO enrichment analysis (Biological Process) on up- and
#              down-regulated genes from the RNA-seq LRT timecourse; osmotic
#              stress-focused barplot of top 5 selected terms per direction
# ##############################################################################

# Libraries ----
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(enrichplot)

# Parameters ----
dds_rds               <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
tables_dir            <- "tables"
output_plot_pdf       <- "plots/go_enrichment_osmotic_stress.pdf"
output_upreg_csv      <- file.path(tables_dir, "go_enrichment_upregulated.csv")
output_downreg_csv    <- file.path(tables_dir, "go_enrichment_downregulated.csv")
output_osmotic_csv    <- file.path(tables_dir, "go_enrichment_osmotic_stress_data.csv")
output_annotated_csv  <- file.path(tables_dir, "go_enrichment_osmotic_stress_annotated.csv")

padj_threshold    <- 0.05
lfc_threshold     <- 2
go_ontology       <- "BP"
go_pvalue_cutoff  <- 0.05
go_qvalue_cutoff  <- 0.2
go_padjust_method <- "BH"
go_min_gset_size  <- 10
go_max_gset_size  <- 500
top_n_terms       <- 5
description_width <- 50

# Data import ----
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dds <- readRDS(dds_rds)

# Analysis ----

## Extract differential genes ----
time_coefficients <- resultsNames(dds) |>
  str_subset("Time")

cat("Found", length(time_coefficients), "time coefficients:\n")
print(time_coefficients)

extract_upreg_genes <- function(coefficient) {
  results(dds, name = coefficient) |>
    as.data.frame() |>
    filter(!is.na(padj), padj < padj_threshold, log2FoldChange > lfc_threshold) |>
    rownames()
}

extract_downreg_genes <- function(coefficient) {
  results(dds, name = coefficient) |>
    as.data.frame() |>
    filter(!is.na(padj), padj < padj_threshold, log2FoldChange < -lfc_threshold) |>
    rownames()
}

upreg_list   <- map(time_coefficients, extract_upreg_genes)
downreg_list <- map(time_coefficients, extract_downreg_genes)

cat("\nUpregulated genes per timepoint:\n")
walk2(time_coefficients, upreg_list,   ~cat(.x, ":", length(.y), "genes\n"))
cat("\nDownregulated genes per timepoint:\n")
walk2(time_coefficients, downreg_list, ~cat(.x, ":", length(.y), "genes\n"))

all_upreg_genes   <- unlist(upreg_list)   |> unique()
all_downreg_genes <- unlist(downreg_list) |> unique()

cat("\nTotal unique upregulated genes:",   length(all_upreg_genes),   "\n")
cat("Total unique downregulated genes:", length(all_downreg_genes), "\n")

## Convert gene IDs to Entrez IDs ----
gene_info <- rowData(dds) |> as.data.frame()

if (!"gene_id" %in% colnames(gene_info)) {
  gene_info <- gene_info |> rownames_to_column("gene_id")
}

gene_info <- gene_info |> dplyr::select(gene_id, symbol)

convert_to_entrez <- function(gene_ids) {
  gene_info |>
    filter(gene_id %in% gene_ids) |>
    pull(symbol) |>
    (\(syms) AnnotationDbi::select(
      org.Hs.eg.db, keys = syms, keytype = "SYMBOL", columns = "ENTREZID"
    ))() |>
    pull(ENTREZID) |>
    na.omit() |>
    unique()
}

upreg_entrez      <- convert_to_entrez(all_upreg_genes)
downreg_entrez    <- convert_to_entrez(all_downreg_genes)
background_entrez <- convert_to_entrez(gene_info$gene_id)

cat("\nUpregulated:   converted to", length(upreg_entrez),      "Entrez IDs\n")
cat("Downregulated: converted to", length(downreg_entrez),    "Entrez IDs\n")
cat("Background:    converted to", length(background_entrez), "Entrez IDs\n")

## GO enrichment ----
run_go <- function(gene_entrez) {
  enrichGO(
    gene          = gene_entrez,
    universe      = background_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = go_ontology,
    pAdjustMethod = go_padjust_method,
    pvalueCutoff  = go_pvalue_cutoff,
    qvalueCutoff  = go_qvalue_cutoff,
    minGSSize     = go_min_gset_size,
    maxGSSize     = go_max_gset_size,
    readable      = TRUE
  )
}

go_upreg   <- run_go(upreg_entrez)
go_downreg <- run_go(downreg_entrez)

cat("Enriched terms (up):",   nrow(go_upreg@result),   "\n")
cat("Enriched terms (down):", nrow(go_downreg@result), "\n")

## Select osmotic stress-relevant terms ----
selected_upreg_ids <- c(
  "GO:0006814",  # sodium ion transport
  "GO:0009612",  # response to mechanical stimulus
  "GO:0007159",  # leukocyte cell-cell adhesion
  "GO:0050727",  # regulation of inflammatory response
  "GO:0001525"   # angiogenesis
)

selected_downreg_ids <- c(
  "GO:0008630",  # intrinsic apoptotic signaling in response to DNA damage
  "GO:0006448",  # regulation of translational elongation
  "GO:0009451",  # RNA modification
  "GO:0006400",  # tRNA modification
  "GO:0002097"   # tRNA wobble base modification
)

prepare_go_data <- function(go_result, selected_ids, category_label) {
  go_result@result |>
    as.data.frame() |>
    filter(ID %in% selected_ids) |>
    arrange(p.adjust) |>
    head(top_n_terms) |>
    mutate(
      log_padj    = -log10(p.adjust),
      Description = str_trunc(Description, width = description_width),
      category    = category_label
    )
}

go_upreg_data   <- prepare_go_data(go_upreg,   selected_upreg_ids,   "Upregulated Genes")
go_downreg_data <- prepare_go_data(go_downreg, selected_downreg_ids, "Downregulated Genes")

cat("\n=== Top 5 osmotic stress-relevant GO terms ===\n")
cat("Upregulated:\n")
go_upreg_data   |> dplyr::select(Description, p.adjust, Count) |> print()
cat("\nDownregulated (NOTE: weak significance, n =", length(all_downreg_genes), "genes):\n")
go_downreg_data |> dplyr::select(Description, p.adjust, Count) |> print()

combined_data <- bind_rows(go_upreg_data, go_downreg_data) |>
  mutate(
    Description = fct_reorder(Description, log_padj),
    category    = factor(category, levels = c("Upregulated Genes", "Downregulated Genes"))
  )

# Visualization ----
go_plot <- ggplot(combined_data, aes(x = log_padj, y = Description, fill = category)) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(x = 0.3, label = Description),
            hjust = 0, size = 2.5, color = "black", fontface = "plain") +
  facet_wrap(~category, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Upregulated Genes"   = "#F8766D",
                               "Downregulated Genes" = "#619CFF")) +
  labs(x = expression(-log[10]("adjusted p-value")),
       y = NULL, title = "GO Enrichment Analysis") +
  theme_minimal() +
  theme(
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(size = 8, color = "black"),
    axis.title.x       = element_text(size = 9, face = "bold"),
    axis.ticks.y       = element_blank(),
    axis.line.x        = element_line(color = "black", linewidth = 0.4),
    axis.line.y        = element_blank(),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    strip.text         = element_text(size = 9, face = "bold", hjust = 0),
    strip.background   = element_blank(),
    legend.position    = "none"
  )

# Save outputs ----
write_csv(go_upreg@result,   output_upreg_csv)
write_csv(go_downreg@result, output_downreg_csv)
write_csv(combined_data,     output_osmotic_csv)

combined_data_annotated <- combined_data |>
  mutate(
    pathway_group = case_when(
      ID == "GO:0006814" ~ "Ion Homeostasis",
      ID == "GO:0009612" ~ "Mechanical Stress",
      ID == "GO:0007159" ~ "Cell Adhesion",
      ID == "GO:0050727" ~ "Inflammatory Response",
      ID == "GO:0001525" ~ "Angiogenesis",
      ID == "GO:0008630" ~ "DNA Damage Response",
      ID == "GO:0006448" ~ "Translation Regulation",
      ID %in% c("GO:0009451", "GO:0006400", "GO:0002097") ~ "RNA Metabolism",
      TRUE ~ "Other"
    )
  )

write_csv(combined_data_annotated, output_annotated_csv)

ggsave(output_plot_pdf, go_plot, width = 8, height = 6, units = "in")

sessionInfo()

## RNA-seq Timecourse Analysis

## Load required libraries
library(data.table)
library(tximeta)
library(DESeq2)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BiocManager)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
library(DEGreport)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(AnnotationDbi)
library(ComplexHeatmap)
library(dplyr)

## Load utility scripts
source("scripts/utils/make_norm_matrix.R")

## Read in sample sheet as colData
colData <-
  fread("data/processed/rna/timecourse/output/STRS_HEK293_WT_1_RNApipeSamplesheet.txt") |> 
  as.data.frame()

## Edit quant paths 
## *note: necessary for `tximeta`
colData$files <- list.files(
  "data/processed/rna/timecourse/output/quant",
  full.names = T,
  recursive = T,
  pattern = "quant.sf")

## rename data.table column `sn` to `names` 
## *note: necessary for `tximeta`
setnames(colData, "sn", "names")

## Check that quant paths are correct
file.exists(colData$files)

## Import data with tximeta & summarize to gene
se <- tximeta(colData)
gse <- summarizeToGene(se)

## Convert all colData to factors
colData(gse)[] <- lapply(colData(gse), as.factor)

## Build DESeq object 
## *note: Left out Treatment so that matrix is in full rank (no linear combos)
dds <- DESeqDataSet(gse, design = ~Bio_Rep + Time)

## Filter out lowly expressed genes (at least 10 counts in at least 4 samples)
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]

## Fit model *note: Using LRT test because this is a timecourse analysis
dds <- DESeq(dds, test = "LRT", reduced = ~Bio_Rep)
res <- results(dds)
dds$Time <- fct_relevel(dds$Time, c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))

## Define coefs
coef <- resultsNames(dds) |> 
  str_subset("Time")

## Find number of differential + significant genes
numDiff <- list()
for (i in seq_along(coef)) {
  diffGenes <- results(dds,
                       name = coef[i]) |>
    as.data.frame() |>
    filter(padj < 0.05 & log2FoldChange > 2 |
             padj < 0.05 & log2FoldChange < -2) |>
    mutate(timepoint = coef[i])
  
  numDiff[[i]] <- diffGenes
}
numDiff_comb <- do.call(rbind, numDiff) |>
  rownames_to_column(var = "ENSEMBL")

## Map ENSEMBL IDs to SYMBOL
geneIDs <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = numDiff_comb$ENSEMBL,
                                 keytype = "ENSEMBL",
                                 columns = c("SYMBOL")) |>
  distinct(ENSEMBL, .keep_all = T)

numDiff_comb <- numDiff_comb |>
  left_join(geneIDs)

## How many differential genes per timepoint?
numDiff_comb <- numDiff_comb |>
  mutate(timepoint = case_when(
    timepoint == "Time_1h_vs_0h" ~ "1h",
    timepoint == "Time_3h_vs_0h" ~ "3h",
    timepoint == "Time_6h_vs_0h" ~ "6h",
    timepoint == "Time_9h_vs_0h" ~ "9h",
    timepoint == "Time_12h_vs_0h" ~ "12h",
    timepoint == "Time_24h_vs_0h" ~ "24h"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 2 ~ "upregulated",
      padj < 0.05 & log2FoldChange < -2 ~ "downregulated"))

pdf("plots/rnaseqTimecourse_LRT_sigGenes_barplot.pdf")
(bp <- numDiff_comb |> 
  dplyr::count(timepoint, regulation) |>
  ggplot(aes(x = timepoint, y = n, fill = fct_rev(regulation))) +
  geom_col() +
  geom_text(aes(label = n),
            vjust = -2.55
  ) +
  labs(x = "",
       y = "",
       title = "# of Differential Genes",
       fill = "") +
  scale_x_discrete(limits = c("1h", "3h", "6h",
                              "9h", "12h", "24h")) +
  theme(
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10,
                               face = "bold"),
    title = element_text(), 
    legend.position = "top",
    axis.ticks.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")))
dev.off()

sig_genes <- dds[(which(res$padj < 0.05))] |> 
  results() |> 
  as.data.frame() |> 
  arrange(padj) |> 
  rownames_to_column(var = "ENSEMBL") |> 
  pull(ENSEMBL) |> 
  head(100)

## QC
plotPCA(vst(dds),
        intgroup = c("Treatment", "Time")) + 
  ggplot2::theme(aspect.ratio = 1)

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

dds <- dds |> 
  keepStandardChromosomes(pruning.mode = "coarse")

seqlevelsStyle(dds) <- "UCSC"
seqlevels(dds) <- seqlevels(txdb)
seqinfo(dds) <- seqinfo(txdb)

## Save DESeq object as .rds
saveRDS(dds, "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

## Calculate key statistics for paper

# Total genes retained after filtering
total_genes_analyzed <- nrow(dds)
cat("Total genes analyzed after filtering:", total_genes_analyzed, "\n")

# Genes with significant time-dependent variation (LRT test)
lrt_sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
cat("Genes with significant time-dependent variation (LRT, FDR < 0.05):", lrt_sig_genes, "\n")

# Peak response timepoint
peak_counts <- numDiff_comb |>
  dplyr::count(timepoint) |>
  arrange(desc(n))
peak_timepoint <- peak_counts$timepoint[1]
peak_n <- peak_counts$n[1]
cat("Peak response timepoint:", peak_timepoint, "with", peak_n, "differential genes\n")

# Earliest significant response
early_response <- numDiff_comb |>
  filter(timepoint == "1h") |>
  nrow()
cat("Differential genes at earliest timepoint (1h):", early_response, "\n")

# Sustained response (genes differential at 24h)
sustained_response <- numDiff_comb |>
  filter(timepoint == "24h") |>
  nrow()
cat("Differential genes at final timepoint (24h):", sustained_response, "\n")

# Genes differential at multiple timepoints (indicating temporal dynamics)
genes_multiple_timepoints <- numDiff_comb |>
  dplyr::count(ENSEMBL) |>
  filter(n > 1) |>
  nrow()
cat("Genes differential at multiple timepoints:", genes_multiple_timepoints, "\n")

# Breakdown by regulation direction
regulation_summary <- numDiff_comb |>
  dplyr::count(regulation)
cat("\nRegulation breakdown across all timepoints:\n")
print(regulation_summary)

# Ratio of up vs down
ratio_up_down <- regulation_summary$n[regulation_summary$regulation == "upregulated"] / 
  regulation_summary$n[regulation_summary$regulation == "downregulated"]
cat("Ratio of upregulated to downregulated genes:", round(ratio_up_down, 1), ":1\n")

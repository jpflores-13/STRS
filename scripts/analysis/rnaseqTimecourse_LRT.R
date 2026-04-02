# ##############################################################################
# filename:    rnaseqTimecourse_LRT.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: DESeq2 LRT timecourse analysis on RNA-seq data; imports salmon
#              quantifications via tximeta, fits ~Bio_Rep + Time model, extracts
#              per-timepoint differentially expressed genes, and saves the dds
# ##############################################################################

# Libraries ----
library(data.table)
library(tximeta)
library(DESeq2)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
library(DEGreport)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(AnnotationDbi)
library(ComplexHeatmap)

# Parameters ----
samplesheet_file <- "data/processed/rna/timecourse/output/STRS_HEK293_WT_1_RNApipeSamplesheet.txt"
quant_dir        <- "data/processed/rna/timecourse/output/quant"
output_dds_rds   <- "data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds"
output_barplot   <- "plots/rnaseqTimecourse_LRT_sigGenes_barplot.pdf"
output_pca       <- "plots/rnaseqTimecourse_LRT_PCA.pdf"
min_counts       <- 50
min_samples      <- 4
padj_cutoff      <- 0.05
lfc_cutoff       <- 2

# Data import ----
source("scripts/utils/make_norm_matrix.R")

colData <- fread(samplesheet_file) |> as.data.frame()

quant_files <- list.files(quant_dir, full.names = TRUE, recursive = TRUE,
                          pattern = "quant.sf")
tp_num      <- as.numeric(sub("h", "", stringr::str_extract(quant_files, "[0-9]+h")))
colData$files <- quant_files[order(tp_num)]

setnames(colData, "sn", "names")
stopifnot(all(file.exists(colData$files)))

se  <- tximeta(colData)
gse <- summarizeToGene(se)
colData(gse)[] <- lapply(colData(gse), as.factor)

# Analysis ----
dds  <- DESeqDataSet(gse, design = ~Bio_Rep + Time)
keep <- rowSums(counts(dds) >= min_counts) >= min_samples
dds  <- dds[keep, ]

dds  <- DESeq(dds, test = "LRT", reduced = ~Bio_Rep)
res  <- results(dds)
dds$Time <- fct_relevel(dds$Time, c("0h", "1h", "3h", "6h", "9h", "12h", "24h"))

coef <- resultsNames(dds) |> str_subset("Time")

numDiff_list <- lapply(coef, function(c) {
  results(dds, name = c) |>
    as.data.frame() |>
    filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) |>
    mutate(timepoint = c)
})

numDiff_comb <- do.call(rbind, numDiff_list) |>
  rownames_to_column(var = "ENSEMBL")

geneIDs <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = numDiff_comb$ENSEMBL,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) |> dplyr::distinct(ENSEMBL, .keep_all = TRUE)

numDiff_comb <- numDiff_comb |>
  left_join(geneIDs) |>
  mutate(
    timepoint = case_when(
      timepoint == "Time_1h_vs_0h"  ~ "1h",
      timepoint == "Time_3h_vs_0h"  ~ "3h",
      timepoint == "Time_6h_vs_0h"  ~ "6h",
      timepoint == "Time_9h_vs_0h"  ~ "9h",
      timepoint == "Time_12h_vs_0h" ~ "12h",
      timepoint == "Time_24h_vs_0h" ~ "24h"
    ),
    regulation = case_when(
      padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "upregulated",
      padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "downregulated"
    )
  )

## Summary statistics ----
cat("Total genes analyzed after filtering:", nrow(dds), "\n")
cat("Genes with significant time-dependent variation (LRT, FDR < 0.05):",
    sum(res$padj < 0.05, na.rm = TRUE), "\n")

peak_counts    <- numDiff_comb |> dplyr::count(timepoint) |> arrange(desc(n))
cat("Peak response timepoint:", peak_counts$timepoint[1],
    "with", peak_counts$n[1], "differential genes\n")

cat("Differential at 1h:", nrow(filter(numDiff_comb, timepoint == "1h")), "\n")
cat("Differential at 24h:", nrow(filter(numDiff_comb, timepoint == "24h")), "\n")

genes_multi_tp <- numDiff_comb |> dplyr::count(ENSEMBL) |> filter(n > 1) |> nrow()
cat("Genes differential at multiple timepoints:", genes_multi_tp, "\n")

reg_summary  <- numDiff_comb |> dplyr::count(regulation)
print(reg_summary)

n_up   <- reg_summary$n[reg_summary$regulation == "upregulated"]
n_down <- reg_summary$n[reg_summary$regulation == "downregulated"]
cat("Up:down ratio:", round(n_up / n_down, 1), ":1\n")

## Add seqinfo ----
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene |> keepStandardChromosomes()
dds  <- keepStandardChromosomes(dds, pruning.mode = "coarse")
seqlevelsStyle(dds) <- "UCSC"
seqlevels(dds) <- seqlevels(txdb)
seqinfo(dds)   <- seqinfo(txdb)

# Visualization ----
bp <- numDiff_comb |>
  dplyr::count(timepoint, regulation) |>
  ggplot(aes(x = timepoint, y = n, fill = fct_rev(regulation))) +
  geom_col() +
  geom_text(aes(label = n), vjust = -2.55) +
  scale_x_discrete(limits = c("1h", "3h", "6h", "9h", "12h", "24h")) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  labs(x = "", y = "", title = "# of Differential Genes", fill = "") +
  theme(axis.text.x     = element_text(size = 10, face = "bold"),
        legend.position = "top",
        axis.ticks.x    = element_blank(),
        panel.background = element_blank())

pca <- plotPCA(vst(dds), intgroup = c("Treatment", "Time")) +
  ggplot2::theme(aspect.ratio = 1)

# Save outputs ----
saveRDS(dds, output_dds_rds)

pdf(output_barplot)
print(bp)
dev.off()

pdf(output_pca)
print(pca)
dev.off()

sessionInfo()

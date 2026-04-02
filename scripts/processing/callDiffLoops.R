# ##############################################################################
# filename:    callDiffLoops.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Perform differential loop analysis with DESeq2 using Wald test
#              and apeglm LFC shrinkage on eGFP-YAP HEK293 Hi-C data
# ##############################################################################

# Libraries ----
library(DESeq2)
library(InteractionSet)
library(glue)
library(stringr)
library(purrr)
library(plyranges)
library(mariner)
library(tidyverse)

# Parameters ----
loop_counts_rds <- "data/processed/hic/h5/sipLoops_eGFP-YAP_noDroso_10kb.rds"
ma_plot_pdf     <- "plots/diffLoops_eGFP-YAP_noDroso_10kb_MA.pdf"
output_rds      <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"

# Data import ----
loopCounts <- readRDS(loop_counts_rds)

# Analysis ----

## Extract count matrix for countData
m <- counts(loopCounts)

## Construct colData/metadata
colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")),
                         stringsAsFactors = TRUE)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Genotype",
                       "Lab", "Treatment", "Replicate")
colData <- colData[c(1:6)]

## Rename Replicate
colData <- colData |>
  mutate(Replicate = factor(rep(c(1:4), 2)))

## Make sure sample names match
all(colnames(m) == rownames(colData))

## Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## Disable DESeq's default normalization
sizeFactors(dds) <- rep(1, ncol(dds))

## Hypothesis testing with Wald with betaPrior = FALSE
dds <- DESeq(dds)

# Visualization ----

## QC plots before saving
plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

## Results with LFC shrinkage
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef = "Treatment_sorbitol_vs_control", type = "apeglm")

pdf(file = ma_plot_pdf)
plotMA(res, ylim = c(-4, 4), main = "Differential Loop Analysis",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Save outputs ----

## Concatenate loopCounts and DESeqResults
rowData(loopCounts) <- cbind(rowData(loopCounts), res)
diff_loopCounts <- loopCounts

base::saveRDS(diff_loopCounts, file = output_rds)

sessionInfo()

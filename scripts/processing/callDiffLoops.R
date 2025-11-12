## Perform differential loops analysis with DESeq2

library(DESeq2)
library(InteractionSet)
library(glue)
library(stringr)
library(purrr)
library(plyranges)
library(mariner)
library(tidyverse)

# Load extracted loops counts ---------------------------------------------
loopCounts <- readRDS("data/processed/hic/h5/sipLoops_eGFP-YAP_noDroso_10kb.rds")

# Extract count matrix for countData -------------------------------------------
m <- counts(loopCounts)

# Construct colData/metadata ----------------------------------------------

## String split the colnames 
colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")),
                         stringsAsFactors = T)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Genotype",
                       "Lab", "Treatment", "Replicate")
colData <- colData[c(1:6)]

## Rename Replicate
colData <- colData |> 
  mutate(Replicate = factor(rep(c(1:4), 2)))

## Make sure sample names match
all(colnames(m) == rownames(colData))

# Run DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## disable DESeq's default normalization 
sizeFactors(dds) <- rep(1, ncol(dds))

## Hypothesis testing with Wald with `betaPrior = F`
dds <- DESeq(dds)

# QC via data visualization before moving on ------------------------------

## plot PCA
plotPCA(vst(dds), intgroup = "Treatment") + ggplot2::theme(aspect.ratio = 1)
plotPCA(vst(dds), intgroup = "Treatment", returnData = TRUE)

## outputs
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")

pdf(file = "plots/diffLoops_eGFP-YAP_noDroso_10kb_MA.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential Loop Analysis",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Concatenate loopCounts and DESeqResults ---------------------------------
rowData(loopCounts) <- cbind(rowData(loopCounts), res)
diff_loopCounts <- loopCounts

## save as .rds
base::saveRDS(diff_loopCounts, file = "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")

## print sessionInfo
sessionInfo()
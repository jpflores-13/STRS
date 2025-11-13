## Create eGFP-YAP lost loop APA plots comparing 
## HEK293 eGFP-YAP/WT, HCT116, & T47D genotypes

## load packages
library(mariner)
library(plotgardener)
library(glue)
library(RColorBrewer)
library(tidyverse)
library(InteractionSet)

# Load data ---------------------------------------------------------------

## Load HiC files
## Load HEK293 eGFP-YAP overexpression files
oe_hic <- list.files("data/processed/hic/maps",
                     full.names = T,
                     pattern = "eGFP-YAP_Cai") |> 
  str_subset("megaMap", negate = T)

## Load HEK293 WT files  
wt_hic <- list.files("data/processed/hic/maps",
                     full.names = T,
                     pattern = "WT") |> 
  str_subset("megaMap", negate = T)

## Load HCT116 WT files  
hct_hic <- list.files("data/processed/hic/maps",
                      full.names = T,
                      pattern = "STRS_HCT116") |> 
  str_subset("megaMap", negate = T)

## Load HEK293 eGFP-YAPdTAD files
dtad_hic <- list.files("data/processed/hic/maps",
                       full.names = T,
                       pattern = "eGFP-YAPdTAD_Cai")

## Set a vector up for `glue` & load Amat et al 2019 HiC data
cond <- c("cont", "nacl")
amat_hic <- list.files(glue("/users/j/p/jpflores/projects/YAPP/MYAP/external/HYPE/data/raw/hic/hg38/220717_dietJuicerCore/{cond}"),
                       full.names = T)

## Load differential loop calls
diffLoops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
lostLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                               rowData(diffLoops)$log2FoldChange < 0)]
lostLoops <- interactions(lostLoops)

# calculate APA matrices --------------------------------------------------
## Parameters
buffer <- 10
resolution <- 10e3

## WT genotype
wt_contAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = wt_hic[str_detect(wt_hic, "HEK293_WT_Phanstiel_control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

wt_sorbAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = wt_hic[str_detect(wt_hic, "HEK293_WT_Phanstiel_sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## EGFP-YAP genotype
oe_contAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = oe_hic[str_detect(oe_hic, "control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

oe_sorbAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = oe_hic[str_detect(oe_hic, "sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## HCT116 Cells
hct_contAPA <- lostLoops |>
  pixelsToMatrices(buffer = 10) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = 10e3,
                  files = hct_hic[str_detect(hct_hic, "control")], ##cont
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

hct_sorbAPA <- lostLoops |>
  pixelsToMatrices(buffer = 10) |>
  removeShortPairs() |> 
  pullHicMatrices(binSize = 10e3,
                  files = hct_hic[str_detect(hct_hic, "sorbitol")], ##sorb
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## T47D cells
amat_contAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = amat_hic[str_detect(amat_hic, "None")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

amat_sorbAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = amat_hic[str_detect(amat_hic, "NaCl")], 
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## eGFP-YAPdTAD cells
dtad_contAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = dtad_hic[str_detect(dtad_hic, "control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

dtad_sorbAPA <- lostLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = dtad_hic[str_detect(dtad_hic, "sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## Divide each genotype/condition by nLoops
nLoops <- length(lostLoops)

wt_contAPA <- (wt_contAPA/nLoops)
wt_sorbAPA <- (wt_sorbAPA/nLoops)
oe_contAPA <- (oe_contAPA/nLoops)
oe_sorbAPA <- (oe_sorbAPA/nLoops)
hct_contAPA <- (hct_contAPA/nLoops)
hct_sorbAPA <- (hct_sorbAPA/nLoops)
amat_contAPA <- (amat_contAPA/nLoops)
amat_sorbAPA <- (amat_sorbAPA/nLoops)
dtad_contAPA <- (dtad_contAPA/nLoops)
dtad_sorbAPA <- (dtad_sorbAPA/nLoops)

## Save normalized APAs as `.rds` files
saveRDS(wt_contAPA, file = "data/processed/hic/normalizedAPA/wt_lostLoops_contAPA_normalized.rds")
saveRDS(wt_sorbAPA, file = "data/processed/hic/normalizedAPA/wt_lostLoops_sorbAPA_normalized.rds")
saveRDS(oe_contAPA, file = "data/processed/hic/normalizedAPA/oe_lostLoops_contAPA_normalized.rds")
saveRDS(oe_sorbAPA, file = "data/processed/hic/normalizedAPA/oe_lostLoops_sorbAPA_normalized.rds")
saveRDS(hct_contAPA, file = "data/processed/hic/normalizedAPA/hct_lostLoops_contAPA_normalized.rds")
saveRDS(hct_sorbAPA, file = "data/processed/hic/normalizedAPA/hct_lostLoops_sorbAPA_normalized.rds")
saveRDS(amat_contAPA, file = "data/processed/hic/normalizedAPA/amat_lostLoops_contAPA_normalized.rds")
saveRDS(amat_sorbAPA, file = "data/processed/hic/normalizedAPA/amat_lostLoops_sorbAPA_normalized.rds")
saveRDS(dtad_contAPA, file = "data/processed/hic/normalizedAPA/dtad_lostLoops_contAPA_normalized.rds")
saveRDS(dtad_sorbAPA, file = "data/processed/hic/normalizedAPA/dtad_lostLoops_sorbAPA_normalized.rds")

## print sessionInfo
sessionInfo()
## Create eGFP-YAP gained loop APA plots comparing 
## HEK293 eGFP-YAP/WT/dTAD, HCT116 WT/mAID2-CTCF/mAID2-RAD21, & T47D genotypes

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
                     full.names = TRUE,
                     pattern = "eGFP-YAP_Cai") |> 
  str_subset("megaMap", negate = TRUE)

## Load HEK293 WT files  
wt_hic <- list.files("data/processed/hic/maps",
                     full.names = TRUE,
                     pattern = "WT") |> 
  str_subset("megaMap", negate = TRUE)

## Load HCT116 WT files  
hct_hic <- list.files("data/processed/hic/maps",
                      full.names = TRUE,
                      pattern = "STRS_HCT116") |> 
  str_subset("megaMap", negate = TRUE)

## Load HCT116 mAID2-CTCF files
ctcf_hic <- list.files("data/processed/hic/maps",
                       full.names = TRUE,
                       pattern = "mAID2-CTCF") |> 
  str_subset("megaMap", negate = TRUE)

## Load HCT116 mAID2-RAD21 files
rad21_hic <- list.files("data/processed/hic/maps",
                        full.names = TRUE,
                        pattern = "mAID2-RAD21") |> 
  str_subset("megaMap", negate = TRUE)

## Load HEK293 eGFP-YAPdTAD files
dtad_hic <- list.files("data/processed/hic/maps",
                       full.names = TRUE,
                       pattern = "eGFP-YAPdTAD_Cai")

## Set a vector up for `glue` & load Amat et al 2019 HiC data
cond <- c("cont", "nacl")
amat_hic <- list.files(glue("/users/j/p/jpflores/projects/YAPP/MYAP/external/HYPE/data/raw/hic/hg38/220717_dietJuicerCore/{cond}"),
                       full.names = TRUE)

## Load differential loop calls
diffLoops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange > 0)]
gainedLoops <- interactions(gainedLoops)

# calculate APA matrices --------------------------------------------------
## Parameters
buffer <- 10
resolution <- 10e3

## WT genotype
wt_contAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = wt_hic[str_detect(wt_hic, "HEK293_WT_Phanstiel_control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

wt_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = wt_hic[str_detect(wt_hic, "HEK293_WT_Phanstiel_sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## EGFP-YAP genotype
oe_contAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = oe_hic[str_detect(oe_hic, "control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

oe_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = oe_hic[str_detect(oe_hic, "sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## HCT116 Cells
hct_contAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = resolution,
                  files = hct_hic[str_detect(hct_hic, "control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

hct_sorbAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = hct_hic[str_detect(hct_hic, "sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## HCT116 mAID2-CTCF cells
ctcf_sorbAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = resolution,
                  files = ctcf_hic[str_detect(ctcf_hic, "sorbitol") & 
                                     !str_detect(ctcf_hic, "auxin")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

ctcf_sorb_auxAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = resolution,
                  files = ctcf_hic[str_detect(ctcf_hic, "sorbitol") & 
                                     str_detect(ctcf_hic, "auxin")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## HCT116 mAID2-RAD21 cells
rad21_sorbAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = resolution,
                  files = rad21_hic[str_detect(rad21_hic, "sorbitol") & 
                                      !str_detect(rad21_hic, "auxin")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

rad21_sorb_auxAPA <- gainedLoops |>
  pixelsToMatrices(buffer = buffer) |>
  removeShortPairs() |>
  pullHicMatrices(binSize = resolution,
                  files = rad21_hic[str_detect(rad21_hic, "sorbitol") & 
                                      str_detect(rad21_hic, "auxin")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## T47D cells
amat_contAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = amat_hic[str_detect(amat_hic, "None")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

amat_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = amat_hic[str_detect(amat_hic, "NaCl")], 
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## eGFP-YAPdTAD cells
dtad_contAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = dtad_hic[str_detect(dtad_hic, "control")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

dtad_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = buffer) |> 
  removeShortPairs() |> 
  pullHicMatrices(binSize = resolution,
                  files = dtad_hic[str_detect(dtad_hic, "sorbitol")],
                  half = "upper",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum)

## Divide each genotype/condition by nLoops
nLoops <- length(gainedLoops)

wt_contAPA <- (wt_contAPA / nLoops)
wt_sorbAPA <- (wt_sorbAPA / nLoops)
oe_contAPA <- (oe_contAPA / nLoops)
oe_sorbAPA <- (oe_sorbAPA / nLoops)
hct_contAPA <- (hct_contAPA / nLoops)
hct_sorbAPA <- (hct_sorbAPA / nLoops)
ctcf_sorbAPA <- (ctcf_sorbAPA / nLoops)
ctcf_sorb_auxAPA <- (ctcf_sorb_auxAPA / nLoops)
rad21_sorbAPA <- (rad21_sorbAPA / nLoops)
rad21_sorb_auxAPA <- (rad21_sorb_auxAPA / nLoops)
amat_contAPA <- (amat_contAPA / nLoops)
amat_sorbAPA <- (amat_sorbAPA / nLoops)
dtad_contAPA <- (dtad_contAPA / nLoops)
dtad_sorbAPA <- (dtad_sorbAPA / nLoops)

## Save normalized APAs as `.rds` files
saveRDS(wt_contAPA, file = "data/processed/hic/normalizedAPA/wt_gainedLoops_contAPA_normalized.rds")
saveRDS(wt_sorbAPA, file = "data/processed/hic/normalizedAPA/wt_gainedLoops_sorbAPA_normalized.rds")
saveRDS(oe_contAPA, file = "data/processed/hic/normalizedAPA/oe_gainedLoops_contAPA_normalized.rds")
saveRDS(oe_sorbAPA, file = "data/processed/hic/normalizedAPA/oe_gainedLoops_sorbAPA_normalized.rds")
saveRDS(hct_contAPA, file = "data/processed/hic/normalizedAPA/hct_gainedLoops_contAPA_normalized.rds")
saveRDS(hct_sorbAPA, file = "data/processed/hic/normalizedAPA/hct_gainedLoops_sorbAPA_normalized.rds")
saveRDS(ctcf_sorbAPA, file = "data/processed/hic/normalizedAPA/ctcf_gainedLoops_sorbAPA_normalized.rds")
saveRDS(ctcf_sorb_auxAPA, file = "data/processed/hic/normalizedAPA/ctcf_gainedLoops_sorb_auxAPA_normalized.rds")
saveRDS(rad21_sorbAPA, file = "data/processed/hic/normalizedAPA/rad21_gainedLoops_sorbAPA_normalized.rds")
saveRDS(rad21_sorb_auxAPA, file = "data/processed/hic/normalizedAPA/rad21_gainedLoops_sorb_auxAPA_normalized.rds")
saveRDS(amat_contAPA, file = "data/processed/hic/normalizedAPA/amat_gainedLoops_contAPA_normalized.rds")
saveRDS(amat_sorbAPA, file = "data/processed/hic/normalizedAPA/amat_gainedLoops_sorbAPA_normalized.rds")
saveRDS(dtad_contAPA, file = "data/processed/hic/normalizedAPA/dtad_gainedLoops_contAPA_normalized.rds")
saveRDS(dtad_sorbAPA, file = "data/processed/hic/normalizedAPA/dtad_gainedLoops_sorbAPA_normalized.rds")

## print sessionInfo
sessionInfo()
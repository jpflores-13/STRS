# ##############################################################################
# filename:    calcAPA_gained.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Calculate normalized APA matrices for eGFP-YAP gained loops
#              across HEK293 eGFP-YAP/WT/dTAD, HCT116 WT/mAID2-CTCF/mAID2-RAD21,
#              and T47D genotypes under control and sorbitol conditions
# ##############################################################################

# Libraries ----
library(mariner)
library(plotgardener)
library(glue)
library(RColorBrewer)
library(tidyverse)
library(InteractionSet)

# Parameters ----
hic_dir          <- "data/processed/hic/maps"
amat_hic_base    <- "/users/j/p/jpflores/projects/YAPP/MYAP/external/HYPE/data/raw/hic/hg38/220717_dietJuicerCore"
diff_loops_rds   <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_dir       <- "data/processed/hic/normalizedAPA"
buffer           <- 10
resolution       <- 10e3

# Data import ----

## Load HEK293 eGFP-YAP overexpression files
oe_hic <- list.files(hic_dir, full.names = TRUE, pattern = "eGFP-YAP_Cai") |>
  str_subset("megaMap", negate = TRUE)

## Load HEK293 WT files
wt_hic <- list.files(hic_dir, full.names = TRUE, pattern = "WT") |>
  str_subset("megaMap", negate = TRUE)

## Load HCT116 WT files
hct_hic <- list.files(hic_dir, full.names = TRUE, pattern = "STRS_HCT116") |>
  str_subset("megaMap", negate = TRUE)

## Load HCT116 mAID2-CTCF files
ctcf_hic <- list.files(hic_dir, full.names = TRUE, pattern = "mAID2-CTCF") |>
  str_subset("megaMap", negate = TRUE)

## Load HCT116 mAID2-RAD21 files
rad21_hic <- list.files(hic_dir, full.names = TRUE, pattern = "mAID2-RAD21") |>
  str_subset("megaMap", negate = TRUE)

## Load HEK293 eGFP-YAPdTAD files
dtad_hic <- list.files(hic_dir, full.names = TRUE, pattern = "eGFP-YAPdTAD_Cai")

## Load Amat et al 2019 HiC data (T47D)
cond     <- c("cont", "nacl")
amat_hic <- list.files(glue("{amat_hic_base}/{cond}"), full.names = TRUE)

## Load differential loop calls
diffLoops <- readRDS(diff_loops_rds)
gainedLoops <- diffLoops[which(rowData(diffLoops)$padj < 0.1 &
                                 rowData(diffLoops)$log2FoldChange > 0)]
gainedLoops <- interactions(gainedLoops)

# Analysis ----

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

## eGFP-YAP genotype
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

## HCT116 cells
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

## Normalize by loop count
nLoops <- length(gainedLoops)

wt_contAPA        <- (wt_contAPA / nLoops)
wt_sorbAPA        <- (wt_sorbAPA / nLoops)
oe_contAPA        <- (oe_contAPA / nLoops)
oe_sorbAPA        <- (oe_sorbAPA / nLoops)
hct_contAPA       <- (hct_contAPA / nLoops)
hct_sorbAPA       <- (hct_sorbAPA / nLoops)
ctcf_sorbAPA      <- (ctcf_sorbAPA / nLoops)
ctcf_sorb_auxAPA  <- (ctcf_sorb_auxAPA / nLoops)
rad21_sorbAPA     <- (rad21_sorbAPA / nLoops)
rad21_sorb_auxAPA <- (rad21_sorb_auxAPA / nLoops)
amat_contAPA      <- (amat_contAPA / nLoops)
amat_sorbAPA      <- (amat_sorbAPA / nLoops)
dtad_contAPA      <- (dtad_contAPA / nLoops)
dtad_sorbAPA      <- (dtad_sorbAPA / nLoops)

# Save outputs ----
saveRDS(wt_contAPA,        file.path(output_dir, "wt_gainedLoops_contAPA_normalized.rds"))
saveRDS(wt_sorbAPA,        file.path(output_dir, "wt_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(oe_contAPA,        file.path(output_dir, "oe_gainedLoops_contAPA_normalized.rds"))
saveRDS(oe_sorbAPA,        file.path(output_dir, "oe_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(hct_contAPA,       file.path(output_dir, "hct_gainedLoops_contAPA_normalized.rds"))
saveRDS(hct_sorbAPA,       file.path(output_dir, "hct_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(ctcf_sorbAPA,      file.path(output_dir, "ctcf_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(ctcf_sorb_auxAPA,  file.path(output_dir, "ctcf_gainedLoops_sorb_auxAPA_normalized.rds"))
saveRDS(rad21_sorbAPA,     file.path(output_dir, "rad21_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(rad21_sorb_auxAPA, file.path(output_dir, "rad21_gainedLoops_sorb_auxAPA_normalized.rds"))
saveRDS(amat_contAPA,      file.path(output_dir, "amat_gainedLoops_contAPA_normalized.rds"))
saveRDS(amat_sorbAPA,      file.path(output_dir, "amat_gainedLoops_sorbAPA_normalized.rds"))
saveRDS(dtad_contAPA,      file.path(output_dir, "dtad_gainedLoops_contAPA_normalized.rds"))
saveRDS(dtad_sorbAPA,      file.path(output_dir, "dtad_gainedLoops_sorbAPA_normalized.rds"))

sessionInfo()

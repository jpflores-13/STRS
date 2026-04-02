# ##############################################################################
# filename:    extractLoops.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Convert SIP loop calls to GInteractions, merge pairs, assign to
#              bins, and extract Hi-C pixel counts for eGFP-YAP HEK293 cells
#              (noDroso SIP setting)
# ##############################################################################

# Libraries ----
library(data.table)
library(mariner)
library(glue)
library(stringr)
library(dbscan)
library(GenomeInfoDb)
library(InteractionSet)

# Parameters ----
loop_dir   <- "data/processed/hic/loops"
hic_dir    <- "data/processed/hic/maps"
h5_out     <- "data/processed/hic/h5/sipLoops_eGFP-YAP_noDroso_10kb.h5"
rds_out    <- gsub(".h5$", ".rds", h5_out)
resolution <- 10e3

# Data import ----

## "noDroso" loops are loops called with SIP setting `-isDroso false`
noDroso_loops <- list.files(path = loop_dir, full.names = TRUE, pattern = "Cai")

hic_files <- list.files(path = hic_dir, full.names = TRUE, pattern = "eGFP-YAP_Cai") |>
  str_subset("megaMap", negate = TRUE)

# Analysis ----

## Convert SIP loop calls to GInteractions
noDroso_loops <- lapply(noDroso_loops, fread) |>
  lapply(as_ginteractions)

## Extract counts
loopCounts_10kb <- noDroso_loops |>
  mergePairs(radius = resolution,
             column = "APScoreAvg") |>
  assignToBins(binSize = resolution) |>
  pullHicPixels(binSize = resolution,
                files = hic_files,
                h5File = h5_out,
                norm = "NONE",
                matrix = "observed")

# Save outputs ----
base::saveRDS(loopCounts_10kb, file = rds_out)

sessionInfo()

## Generate loop counts for control &sorbitol-treated 
## eGFP-YAP HEK293 cells processed with SIP parameter `-isDroso` false

library(data.table)
library(mariner)
library(glue)
library(stringr)
library(dbscan)
library(GenomeInfoDb)
library(InteractionSet)

# Compile loop files of interest ------------------------------------------

## "noDroso" loops are loops called with SIP setting `-isDroso false`
noDroso_loops <- list.files(path = "data/processed/hic/loops",
                            full.names = T,
                            pattern = "Cai")

# Provide .hic files to extract from --------------------------------------

hic_files <- list.files(path = "data/processed/hic/maps",
                        full.names = T,
                        pattern = "eGFP-YAP_Cai") |> 
  str_subset("megaMap", negate = T)

# Merge bedpe files into one and extract counts  --------------------------

## parameters
resolution <- 10e3

## Convert SIP loop calls to GInteractions
noDroso_loops <- lapply(noDroso_loops, fread) |> 
  lapply(as_ginteractions) 

## Extract counts
h5 <- "data/processed/hic/h5/sipLoops_eGFP-YAP_noDroso_10kb.h5"
rds <- gsub(".h5$", ".rds", h5)
loopCounts_10kb <- noDroso_loops |>
  mergePairs(radius = resolution,
             column = "APScoreAvg") |>
  assignToBins(binSize = resolution) |>
  pullHicPixels(binSize = resolution,
                files = hic_files,
                h5File = h5,
                norm = "NONE",
                matrix = "observed")

## Save files
base::saveRDS(loopCounts_10kb, file = rds)

## print sessionInfo
sessionInfo()
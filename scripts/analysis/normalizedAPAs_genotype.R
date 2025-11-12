## Figure 1 — Panel C only
## Output: figures/Figure1_panelC.pdf (8.5" wide x 2.8" tall)

## ---- Packages ----
library(plotgardener)
library(RColorBrewer)
library(readr)
library(glue)
library(tidyverse)

## ---- I/O ----
out_file <- "plots/normalizedAPAs_genotype.pdf"
depth_csv <- "data/raw/hic/sequencing_depth_results/eight_condition_averages.csv"

## Load normalized APA matrices (already normalized in your pipeline)
normalizedAPAs <- list.files("data/processed/hic/normalizedAPA", full.names = TRUE) |>
  lapply(readRDS)
names(normalizedAPAs) <- list.files("data/processed/hic/normalizedAPA/", recursive = TRUE)

## ---- Read depth helpers ----
depth_data <- read_csv(depth_csv, show_col_types = FALSE)

get_depth_millions <- function(cell_type, condition) {
  condition_standard <- dplyr::case_when(
    condition %in% c("Control", "Untreated", "untreated") ~ "Control",
    condition %in% c("Sorbitol", "NaCl_Treated", "Treated", "+ sorbitol", "+ NaCl") ~ "Treated",
    TRUE ~ condition
  )
  
  depth_row <- depth_data |>
    dplyr::filter(Cell_Type == cell_type, condition_standard == condition_standard)
  
  if (nrow(depth_row) > 0) {
    depth_millions <- depth_row$average_depth_per_replicate
    return(sprintf("%.0fM", depth_millions))
  } else {
    warning(glue("Could not find depth for {cell_type} - {condition_standard}"))
    return("N/A")
  }
}

depths <- list(
  hek_egfp = get_depth_millions("HEK293_eGFP_YAP", "Control"),
  hek_wt   = get_depth_millions("HEK293_WT", "Control"),
  hct116   = get_depth_millions("HCT116", "Control"),
  t47d     = get_depth_millions("T47D", "Control")
)

## ---- Plot params ----
apaParams <- pgParams(
  assembly = "hg38",
  height   = 0.75,
  width    = 0.75,
  fontsize = 5,
  palette  = colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
)

## Precompute zrange maxima per panel column (gained vs lost, per genotype)
wt_gained_comb   <- c(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
                      normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds)
oe_gained_comb   <- c(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
                      normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds)
hct_gained_comb  <- c(normalizedAPAs$hct_gainedLoops_contAPA_normalized.rds,
                      normalizedAPAs$hct_gainedLoops_sorbAPA_normalized.rds)
amat_gained_comb <- c(normalizedAPAs$amat_gainedLoops_contAPA_normalized.rds,
                      normalizedAPAs$amat_gainedLoops_sorbAPA_normalized.rds)

wt_lost_comb     <- c(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
                      normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds)
oe_lost_comb     <- c(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
                      normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds)
hct_lost_comb    <- c(normalizedAPAs$hct_lostLoops_contAPA_normalized.rds,
                      normalizedAPAs$hct_lostLoops_sorbAPA_normalized.rds)
amat_lost_comb   <- c(normalizedAPAs$amat_lostLoops_contAPA_normalized.rds,
                      normalizedAPAs$amat_lostLoops_sorbAPA_normalized.rds)

## ---- Device + page ----
pdf(out_file, width = 8.5, height = 2.8)
pageCreate(width = 8.5, height = 2.8, showGuides = FALSE)

## Y positions for the two rows (untreated / osmotic)
y_top  <- 0.10 + 0.25      # first row of APAs
y_bot  <- y_top + 0.80     # second row of APAs
y_hdr  <- y_top - 0.10     # header line/labels
y_badge<- y_bot + 0.75 - 0.05

## ---- GAINED (left block, ~0–3.75 in) ----
# eGFP-YAP
plotMatrix(normalizedAPAs$oe_gainedLoops_contAPA_normalized.rds,
           x = 0.00, y = y_top,
           zrange = c(0, max(oe_gained_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 0.80, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$oe_gainedLoops_sorbAPA_normalized.rds,
           x = 0.00, y = y_bot,
           zrange = c(0, max(oe_gained_comb)),
           params = apaParams)

plotText(label = depths$hek_egfp, fontcolor = "black", fontsize = 4,
         x = 0.75 - 0.05, y = y_badge, just = c("right", "bottom"))

# HEK WT
plotMatrix(normalizedAPAs$wt_gainedLoops_contAPA_normalized.rds,
           x = 0.95, y = y_top,
           zrange = c(0, max(wt_gained_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 1.75, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$wt_gainedLoops_sorbAPA_normalized.rds,
           x = 0.95, y = y_bot,
           zrange = c(0, max(wt_gained_comb)),
           params = apaParams)

plotText(label = depths$hek_wt, fontcolor = "black", fontsize = 4,
         x = 1.70 - 0.05, y = y_badge, just = c("right", "bottom"))

# HCT116
plotMatrix(normalizedAPAs$hct_gainedLoops_contAPA_normalized.rds,
           x = 1.90, y = y_top,
           zrange = c(0, max(hct_gained_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 2.70, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$hct_gainedLoops_sorbAPA_normalized.rds,
           x = 1.90, y = y_bot,
           zrange = c(0, max(hct_gained_comb)),
           params = apaParams)

plotText(label = depths$hct116, fontcolor = "black", fontsize = 4,
         x = 2.65 - 0.05, y = y_badge, just = c("right", "bottom"))

# T47D (Amat)
plotMatrix(normalizedAPAs$amat_gainedLoops_contAPA_normalized.rds,
           x = 2.85, y = y_top,
           zrange = c(0, max(amat_gained_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 3.65, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$amat_gainedLoops_sorbAPA_normalized.rds,
           x = 2.85, y = y_bot,
           zrange = c(0, max(amat_gained_comb)),
           params = apaParams)

plotText(label = depths$t47d, fontcolor = "black", fontsize = 4,
         x = 3.60 - 0.05, y = y_badge, just = c("right", "bottom"))

## Column labels (top)
plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 0.05, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "HEK293 WT", fontcolor = "black",
         x = 1.125, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "HCT116", fontcolor = "black",
         x = 2.15, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "T47D (Amat et al 2019)", fontcolor = "black",
         x = 2.85, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))

## Row labels for gained block
plotText(label = "untreated", fontcolor = "black",
         x = 0.025, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 0.025, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 0.95 + 0.025, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 0.95 + 0.025, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 1.925, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 1.925, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 2.875, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ NaCl", fontcolor = "black",
         x = 2.875, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))

## Orange header line + title for gained
plotSegments(x0 = 0.00, y0 = y_hdr, x1 = 3.75, y1 = y_hdr,
             default.units = "inches", linecolor = "#F8766D", lwd = 1, lty = 1)
plotText(label = "gained loops", x = 1.80, y = y_hdr - 0.10,
         fontcolor = "#F8766D", fontsize = 10)

## ---- LOST (right block, ~3.8–7.6 in) ----
# eGFP-YAP
plotMatrix(normalizedAPAs$oe_lostLoops_contAPA_normalized.rds,
           x = 3.80, y = y_top,
           zrange = c(0, max(oe_lost_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 4.60, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$oe_lostLoops_sorbAPA_normalized.rds,
           x = 3.80, y = y_bot,
           zrange = c(0, max(oe_lost_comb)),
           params = apaParams)

plotText(label = depths$hek_egfp, fontcolor = "black", fontsize = 4,
         x = 4.55 - 0.05, y = y_badge, just = c("right", "bottom"))

# HEK WT
plotMatrix(normalizedAPAs$wt_lostLoops_contAPA_normalized.rds,
           x = 4.75, y = y_top,
           zrange = c(0, max(wt_lost_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 5.55, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$wt_lostLoops_sorbAPA_normalized.rds,
           x = 4.75, y = y_bot,
           zrange = c(0, max(wt_lost_comb)),
           params = apaParams)

plotText(label = depths$hek_wt, fontcolor = "black", fontsize = 4,
         x = 5.50 - 0.05, y = y_badge, just = c("right", "bottom"))

# HCT116
plotMatrix(normalizedAPAs$hct_lostLoops_contAPA_normalized.rds,
           x = 5.70, y = y_top,
           zrange = c(0, max(hct_lost_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 6.50, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$hct_lostLoops_sorbAPA_normalized.rds,
           x = 5.70, y = y_bot,
           zrange = c(0, max(hct_lost_comb)),
           params = apaParams)

plotText(label = depths$hct116, fontcolor = "black", fontsize = 4,
         x = 6.45 - 0.05, y = y_badge, just = c("right", "bottom"))

# T47D (Amat)
plotMatrix(normalizedAPAs$amat_lostLoops_contAPA_normalized.rds,
           x = 6.65, y = y_top,
           zrange = c(0, max(amat_lost_comb)),
           params = apaParams) |>
  annoHeatmapLegend(x = 7.45, y = y_top, width = 0.05, height = 0.75,
                    fontcolor = "black", fontsize = 5)

plotMatrix(normalizedAPAs$amat_lostLoops_sorbAPA_normalized.rds,
           x = 6.65, y = y_bot,
           zrange = c(0, max(amat_lost_comb)),
           params = apaParams)

plotText(label = depths$t47d, fontcolor = "black", fontsize = 4,
         x = 7.40 - 0.05, y = y_badge, just = c("right", "bottom"))

## Column labels (top) for lost block
plotText(label = "HEK293 eGFP-YAP", fontcolor = "black",
         x = 3.85, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "HEK293 WT", fontcolor = "black",
         x = 4.90, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "HCT116", fontcolor = "black",
         x = 5.90, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))
plotText(label = "T47D (Amat et al 2019)", fontcolor = "black",
         x = 6.65, y = y_badge + 0.05, params = apaParams, just = c("top", "left"))

## Row labels for lost block
plotText(label = "untreated", fontcolor = "black",
         x = 3.825, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 3.825, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 4.75 + 0.025, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 4.75 + 0.025, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 5.725, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ sorbitol", fontcolor = "black",
         x = 5.725, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "untreated", fontcolor = "black",
         x = 6.675, y = y_top + 0.025, params = apaParams, just = c("top", "left"))
plotText(label = "+ NaCl", fontcolor = "black",
         x = 6.675, y = y_bot + 0.025, params = apaParams, just = c("top", "left"))

## Blue header line + title for lost
plotSegments(x0 = 3.80, y0 = y_hdr, x1 = 7.60, y1 = y_hdr,
             default.units = "inches", linecolor = "#619CFF", lwd = 1, lty = 1)
plotText(label = "lost loops", x = 5.55, y = y_hdr - 0.10,
         fontcolor = "#619CFF", fontsize = 10)

dev.off()
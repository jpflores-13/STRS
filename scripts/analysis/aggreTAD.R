# ##############################################################################
# filename:    aggreTAD.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Aggregate the contact domain (i.e. the block of enriched contacts
#              inside a loop's anchors, not an insulation-derived TAD) for
#              gained vs matched pre-existing loops and produce domain/loop
#              enrichment plots; uses matchRanges for size- and
#              contact-frequency-matched null set
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
library(data.table)
library(mariner)
library(raster)
library(reshape2)

source("scripts/utils/aggregateLoops.R")
source("scripts/utils/aggregateTAD.R")
source("scripts/utils/plotAggTAD.R")

# Parameters ----
hic_replicates_dir  <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220715_dietJuicerCore/output/"
hic_merged_dir      <- "/proj/phanstiel_lab/Data/processed/YAPP/hic/hg38/220716_dietJuicerMerge_condition/output/"
diff_loops_rds      <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
plots_dir           <- "plots"

# Data import ----

## Load HiC replicate files
hicFiles <- list.files(hic_replicates_dir,
                       full.names = TRUE,
                       recursive = TRUE,
                       pattern = "inter_30.hic")

## Load merged HiC files for human
merged_hicFiles <- list.files(hic_merged_dir,
                              full.names = TRUE,
                              recursive = TRUE,
                              pattern = "inter_30.hic")

## Load loops and extract pixel counts at 10 kb
noDroso_loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "VC_SQRT",
                matrix = "observed")

# Analysis ----

## Create loop size mcol
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## Classify loop types
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## Aggregate contacts per condition
mcols(noDroso_loops)$sorb_contacts <- counts(noDroso_loops)[,"YAPP_HEK_sorbitol_4_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_sorbitol_5_2_inter_30.hic"] + counts(noDroso_loops)[,"YAPP_HEK_sorbitol_6_2_inter_30.hic"]

mcols(noDroso_loops)$cont_contacts <- counts(noDroso_loops)[,"YAPP_HEK_control_1_2_inter_30.hic"] +
  counts(noDroso_loops)[,"YAPP_HEK_control_2_2_inter_30.hic"] + counts(noDroso_loops)[,"YAPP_HEK_control_3_2_inter_30.hic"]

mcols(noDroso_loops)$agg_contacts <- mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

## Log-transform for matchRanges normality assumption
mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)
hist(mcols(noDroso_loops)$loop_size)

mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA
noDroso_loops <- interactions(noDroso_loops) |>
  as.data.frame() |>
  na.omit() |>
  as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))
hist(mcols(noDroso_loops)$agg_contacts)

## matchRanges: select null set matched for size and contact frequency
focal <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost","other")]
pool  <- noDroso_loops[noDroso_loops$loop_type  %in% c("static", "lost","other")]

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size,
                       method = 'stratified',
                       replace = FALSE)

gainedLoops <- focal
ctcfLoops   <- pool

gained_bed <- gainedLoops |>
  as.data.frame() |>
  dplyr::select(c(1,2,3,6,7,8))

n_gained <- nrow(gained_bed)
cat(glue("Number of sorbitol-induced loops: {n_gained}\n"))

nullSet <- nullSet |>
  as.data.frame() |>
  dplyr::select(c(1,2,3,6,7,8))

# Visualization ----

aggtad_gain <- aggregateTAD(loops = gained_bed,
                            hic = merged_hicFiles[2],
                            res = 10e3,
                            buffer = 0.5,
                            norm = "VC_SQRT")

aggtad_match <- aggregateTAD(loops = nullSet,
                              hic = merged_hicFiles[1],
                              res = 10e3,
                              buffer = 0.5,
                              norm = "VC_SQRT")

# Save outputs ----

## matchRanges QC
pdf(file.path(plots_dir, "matchRanges_covariate_agg_contacts.pdf"), width = 6, height = 5)
plotCovariate(nullSet)
dev.off()

pdf(file.path(plots_dir, "matchRanges_covariate_loop_size.pdf"), width = 6, height = 5)
plotCovariate(nullSet, covar = "loop_size")
dev.off()

pdf(file.path(plots_dir, "matchRanges_propensity.pdf"), width = 6, height = 5)
plotPropensity(nullSet, sets = c('f', 'p', 'm'), log = 'x')
dev.off()

## Aggregated TAD heatmaps
pdf(file.path(plots_dir, "aggregatedTADs_gained.pdf"), width = 8, height = 6)
(aggTAD_gainedPlot <- plotAggTAD(aggtad_gain, maxval = aggtad_gain[26,75]*1.2))
dev.off()

pdf(file.path(plots_dir, "aggregatedTADs_matched.pdf"), width = 8, height = 6)
(aggTAD_matchedPlot <- plotAggTAD(aggtad_match, maxval = aggtad_match[26,75]*1.2))
dev.off()

## TAD enrichment line plot
pdf(file.path(plots_dir, "tad_enrichment.pdf"), width = 8, height = 6)
par(mgp=c(3,.3,0))
plot(1, type="n", xaxt="n", yaxt="n", ann=FALSE, xlim=c(1,76), ylim=c(.9, 1.9))
abline(v=c(26,51), lty=2, col="grey90")

data = aggtad_match[cbind(1:76, 25:100)]
bg   = median(data[c(1:20,56:76)])
lines(data/bg, type="l", col="#619CFF")

data = aggtad_gain[cbind(1:76, 25:100)]
bg   = median(data[c(1:20,56:76)])
lines(data/bg, type="l", col="#F8766D")

axis(side=2, las=2, tcl=0.2, at=seq(1, 2.5, 0.5))
axis(side=1, at = c(26,51), labels=c("left","right"))
legend("topright", bty="n", col=c("#619CFF","#F8766D"), legend=c("Pre-Existing","Gained"), lty=1)
dev.off()

## Loop enrichment line plot
pdf(file.path(plots_dir, "loop_enrichment.pdf"), width = 8, height = 6)
plot(1, type="n", xaxt="n", yaxt="n", ann=FALSE, xlim=c(1,50), ylim=c(.9, 4.0))
abline(v=26, lty=2, col="grey90")

data = aggtad_match[cbind(1:51, 50:100)]
bg   = median(data[c(1:10,40:50)])
lines(data/bg, type="l", col="#619CFF")

data = aggtad_gain[cbind(1:50, 50:99)]
bg   = median(data[c(1:10,40:50)])
lines(data/bg, type="l", col="#F8766D")

axis(side=2, las=2, tcl=0.2, at=seq(1, 5.0, 1.0))
axis(side=1, at = 26, labels="loop pixel")
legend("topright", bty="n", col=c("#619CFF","#F8766D"), legend=c("Pre-Existing", "Gained"), lty=1)
dev.off()

## ggplot2 versions
data_list  <- list("Pre-Existing" = aggtad_match, "Gained" = aggtad_gain)
plot_colors <- c("Pre-Existing" = "#619CFF", "Gained" = "#F8766D")

(tad_plot <- create_enrichment_plot(
  data_list = data_list,
  plot_type = "tad",
  colors = plot_colors,
  title = "",
  legend_position = "top"
))

(loop_plot <- create_enrichment_plot(
  data_list = data_list,
  plot_type = "loop",
  colors = plot_colors,
  title = "",
  legend_position = "top"
))

sessionInfo()

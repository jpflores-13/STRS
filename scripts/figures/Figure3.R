## Figure 3 

library(DESeq2)
library(InteractionSet)
library(plotgardener)
library(tidyverse)
library(RColorBrewer)
library(mariner)
library(plyranges)
library(ggpubr)
library(bamsignals)

# Setup -------------------------------------------------------------------

page_width <- 11
page_height <- 7.7

pdf("figures/Figure3.pdf", width = page_width, height = page_height)
pageCreate(width = page_width, height = page_height, showGuides = FALSE)

gray_color <- "#666666"

# Define colors for protein labels in Panel D
ctcf_color <- "#253494"
rad21_color <- "#41B6C4"
yap1_color <- "#238B45"

# Load Data ---------------------------------------------------------------

diff_CTCF <- readRDS("data/processed/cutntag/deseq2/diff_CTCF_counts.rds")
diff_RAD21 <- readRDS("data/processed/cutntag/deseq2/diff_RAD21_counts.rds")
diff_YAP1 <- readRDS("data/processed/cutntag/deseq2/diff_YAP1_counts.rds")

loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |>
  interactions() |> as.data.frame() |> as_ginteractions()

ctcf_control_bw <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw"
ctcf_sorb_bw    <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw"
rad21_control_bw<- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw"
rad21_sorb_bw   <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw"
yap1_control_bw <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_YAP1_cont_0h.bw"
yap1_sorb_bw    <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_YAP1_sorbitol_1h.bw"

# Categorize loops --------------------------------------------------------

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0] |> as.data.frame() |> as_ginteractions()
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0] |> as.data.frame() |> as_ginteractions()

## Region for survey (loop #105) with buffer
loopRegions_gained <- GRanges(
  seqnames = as.character(seqnames(anchors(gainedLoops, "first"))),
  ranges   = IRanges(start = start(anchors(gainedLoops, "first")),
                     end   = end(anchors(gainedLoops, "second"))),
  mcols    = mcols(gainedLoops)
)
buffer <- 150e3
loopRegions_gained_buffed <- loopRegions_gained + buffer

# Bar plot data -----------------------------------------------------------

cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE, pattern = ".narrowPeak") |>
  lapply(read_narrowpeaks)

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

extractAnchors <- function(gi) unique(c(anchors(gi, "first"), anchors(gi, "second")))
calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  n_anchors <- length(anchors)
  boot_props <- replicate(n_bootstrap, {
    boot_idx <- sample(n_anchors, replace = TRUE)
    sum(countOverlaps(anchors[boot_idx], peaks) > 0) / n_anchors
  })
  qs <- quantile(boot_props, c((1-conf_level)/2, 1-(1-conf_level)/2))
  data.frame(lower = qs[1]*100, upper = qs[2]*100)
}
calculateOverlaps <- function(anchors, peaks_list, category) {
  props <- lapply(peaks_list, function(p) sum(countOverlaps(anchors, p) > 0) / length(anchors))
  cis   <- lapply(peaks_list, function(p) calculateCI(anchors, p))
  data.frame(
    Category  = category,
    Target    = sub("_.*", "", names(props)),
    Condition = sub(".*_", "", names(props)),
    Percentage= unlist(props)*100,
    Lower = sapply(cis, \(x) x$lower),
    Upper = sapply(cis, \(x) x$upper)
  )
}
gained_anchors_bar <- extractAnchors(loops[loops$padj < 0.1 & loops$log2FoldChange > 0])
lost_anchors_bar   <- extractAnchors(loops[loops$padj < 0.1 & loops$log2FoldChange < 0])

if (!exists("bar_results")) {
  bar_results <- rbind(
    calculateOverlaps(gained_anchors_bar, cutntag, "Gained"),
    calculateOverlaps(lost_anchors_bar,  cutntag, "Lost")
  )
}
bar_plot_data <- bar_results |>
  mutate(Category = factor(Category, levels = c("Lost","Gained")),
         Target   = factor(Target,   levels = c("CTCF","RAD21","H3K27ac","YAP1")))

# Density analysis data ---------------------------------------------------

peak_list <- list()
for (t in target) for (cond in condition) {
  pattern <- paste0(t, "_", ifelse(cond=="control","cont",cond))
  fp <- list.files("data/processed/cutntag/output/peaks/", full.names = TRUE, pattern = pattern)
  fp <- fp[grepl("\\.narrowPeak$", fp)]
  if (length(fp) == 1) peak_list[[paste0(t,"_",cond)]] <- read_narrowpeaks(fp)
}
bam_files <- character()
for (t in target) for (cond in condition) {
  pattern <- paste0(t, "_", ifelse(cond=="control","cont",cond))
  fp <- list.files("data/processed/cutntag/output/mergeAlign/", full.names = TRUE, pattern = pattern)
  fp <- fp[grepl("\\.bam$", fp)]
  if (length(fp) == 1) bam_files[paste0(t,"_",cond)] <- fp
}
get_between_regions <- function(loops) {
  a1 <- anchors(loops, "first"); a2 <- anchors(loops, "second")
  st <- end(a1) + 1; en <- start(a2) - 1; ok <- en >= st
  if (!any(ok)) return(GRanges())
  GRanges(seqnames = seqnames(a1)[ok],
          ranges   = IRanges(start = st[ok], end = en[ok]),
          loop_id  = paste0(which(ok), "_between"))
}
analyze_regions <- function(regions, peaks_control, peaks_treat, bam_control, bam_treat, region_type, loop_category) {
  ov <- findOverlaps(peaks_control, regions); if (length(ov) == 0) return(data.frame())
  p  <- peaks_control[queryHits(ov)]
  mcols(p)$region_id <- mcols(regions)$loop_id[subjectHits(ov)]
  pid <- paste0(as.character(seqnames(p)), ":", start(p), "-", end(p))
  cs <- bamCount(bam_control, p, paired.end = "midpoint")
  ts <- bamCount(bam_treat,   p, paired.end = "midpoint")
  l2 <- dplyr::case_when(cs==0 & ts==0 ~ 0, cs==0 ~ log2(ts+1), ts==0 ~ -log2(cs+1),
                         TRUE ~ log2((ts+1)/(cs+1)))
  data.frame(loop_category = loop_category, region_type = region_type, region_id = mcols(p)$region_id,
             peak_id = pid, log2FC = l2, control_signal = cs, treat_signal = ts,
             peak_score = mcols(p)$signalValue) |>
    dplyr::group_by(peak_id) |> dplyr::slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}
results_list <- list()
for (protein in c("CTCF","RAD21","YAP1")) {
  cp <- peak_list[[paste0(protein,"_control")]]
  tp <- peak_list[[paste0(protein,"_sorbitol")]]
  cb <- bam_files[paste0(protein,"_control")]
  tb <- bam_files[paste0(protein,"_sorbitol")]
  for (lt in c("Gained","Lost")) {
    sub <- if (lt=="Gained") gainedLoops else lostLoops
    ar <- c(anchors(sub,"first"), anchors(sub,"second"))
    mcols(ar)$loop_id <- c(paste0(seq_along(anchors(sub,"first")),"_anchor1"),
                           paste0(seq_along(anchors(sub,"second")),"_anchor2"))
    br <- get_between_regions(sub)
    r1 <- analyze_regions(ar, cp, tp, cb, tb, "At Anchors", lt)
    r2 <- analyze_regions(br, cp, tp, cb, tb, "Between Anchors", lt)
    results_list[[paste0(protein,"_",lt)]] <- bind_rows(r1, r2) |> mutate(protein = protein)
  }
}
density_data <- bind_rows(results_list)

# Identify retained vs lost peaks -----------------------------------------

identify_peak_changes <- function(control_peaks, sorbitol_peaks, buffer = 1000) {
  # Expand peaks slightly for overlap detection
  control_expanded <- control_peaks + buffer
  sorbitol_expanded <- sorbitol_peaks + buffer
  
  # Find overlaps
  control_overlaps <- countOverlaps(control_peaks, sorbitol_expanded) > 0
  sorbitol_overlaps <- countOverlaps(sorbitol_peaks, control_expanded) > 0
  
  # Retained peaks = present in both conditions
  retained_control <- control_peaks[control_overlaps]
  retained_sorbitol <- sorbitol_peaks[sorbitol_overlaps]
  
  # Lost peaks = present in control but not sorbitol
  lost_peaks <- control_peaks[!control_overlaps]
  
  # Gained peaks = present in sorbitol but not control
  gained_peaks <- sorbitol_peaks[!sorbitol_overlaps]
  
  list(
    retained_control = retained_control,
    retained_sorbitol = retained_sorbitol,
    lost = lost_peaks,
    gained = gained_peaks
  )
}

# Apply to each protein
ctcf_changes <- identify_peak_changes(
  peak_list[["CTCF_control"]], 
  peak_list[["CTCF_sorbitol"]]
)

rad21_changes <- identify_peak_changes(
  peak_list[["RAD21_control"]], 
  peak_list[["RAD21_sorbitol"]]
)

yap1_changes <- identify_peak_changes(
  peak_list[["YAP1_control"]], 
  peak_list[["YAP1_sorbitol"]]
)

# Plot helpers ------------------------------------------------------------

create_ma_data <- function(diff_obj, protein_name) {
  as.data.frame(mcols(diff_obj)) |>
    dplyr::select(baseMean, log2FoldChange, padj) |>
    mutate(isDE = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "Increased",
      log2FoldChange < -1 & padj < 0.05 ~ "Decreased",
      TRUE ~ "Not significant")) |>
    arrange(isDE)
}

create_ma_plot <- function(res_df, protein_name) {
  n_up <- sum(res_df$isDE=="Increased", na.rm=TRUE)
  n_dn <- sum(res_df$isDE=="Decreased", na.rm=TRUE)
  ggplot(res_df, aes(x=baseMean, y=log2FoldChange, color=isDE)) +
    geom_point(alpha=.7, size=.6) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=.3) +
    scale_color_manual(values=c("Increased"="#F8766D","Decreased"="#619CFF","Not significant"="grey80")) +
    ylim(c(-4,4)) + scale_x_log10(breaks=c(1,50,500)) +
    labs(y=paste0(protein_name, " log2\n(sorbitol/control)"), x="mean of normalized counts") +
    theme_classic() +
    theme(legend.position="none", 
          plot.title=element_text(hjust=.5, face="bold", size=9),
          axis.text=element_text(size=7), 
          axis.title=element_text(size=8.5),
          axis.title.y=element_text(angle=90, vjust=.5), 
          axis.line=element_line(linewidth=.3),
          axis.ticks=element_line(linewidth=.3),
          plot.margin=margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio=1) +
    annotate("text", label="Increased", x=max(res_df$baseMean,na.rm=TRUE)*.1, y=3.6, 
             color="#F8766D", fontface="bold", size=7/.pt) +
    annotate("text", label=paste0("n = ", n_up), x=max(res_df$baseMean,na.rm=TRUE)*.1, y=3.1, 
             color="#F8766D", size=6/.pt) +
    annotate("text", label="Decreased", x=max(res_df$baseMean,na.rm=TRUE)*.1, y=-3.4, 
             color="#619CFF", fontface="bold", size=7/.pt) +
    annotate("text", label=paste0("n = ", n_dn), x=max(res_df$baseMean,na.rm=TRUE)*.1, y=-3.9, 
             color="#619CFF", size=6/.pt)
}

create_ma_density_plot <- function(res_df) {
  ggplot(res_df, aes(y=log2FoldChange)) +
    geom_density(color="#9370DB", fill="#9370DB", alpha=.25) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=.3) +
    ylim(c(-4,4)) +  
    theme_classic() +
    theme(legend.position="none", axis.text=element_blank(), axis.title=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          plot.margin=margin(0,0,0,0,"pt"), panel.spacing=unit(0,"pt"))
}

create_bar_plot <- function(d) {
  protein_name <- unique(d$Target)[1]
  ggplot(d, aes(x=Category, y=Percentage, fill=Condition)) +
    geom_hline(yintercept=seq(0,100,25), color="gray90", linetype="dashed") +
    geom_bar(stat="identity", position=position_dodge(.7), width=.6) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), position=position_dodge(.7), width=.25, linewidth=.3) +
    geom_text(aes(label=sprintf("%.1f%%", Percentage)), position=position_dodge(.7), vjust=-1.2, size=6/.pt) +
    scale_fill_manual(values=c(control="#619CFF", sorbitol="#F8766D")) +
    theme_classic() +
    theme(legend.position="none",
          axis.text=element_text(size=7),
          axis.title=element_text(size=8.5),
          axis.title.x=element_text(),
          axis.line=element_line(linewidth=.3),
          axis.ticks=element_line(linewidth=.3),
          plot.margin=margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio=1) +
    labs(y=paste0("% of anchors bound\nby ", protein_name), x="Loop Anchor Type") +
    scale_y_continuous(limits=c(0,100), expand=expansion(mult=c(0,.05)), breaks=seq(0,100,25)) +
    # Add direct text labels inside plot
    annotate("text", x=0.5, y=95, label="control", color="#619CFF", fontface="bold", size=7/.pt, hjust=0) +
    annotate("text", x=0.5, y=88, label="sorbitol", color="#F8766D", fontface="bold", size=7/.pt, hjust=0)
}

create_density_plot <- function(df, protein_name_input) {
  protein_name <- protein_name_input
  df <- df |> mutate(region_type=factor(region_type, levels=c("At Anchors","Between Anchors")))
  
  # Define colors
  at_anchors_color <- "#5DA5DA"
  between_anchors_color <- "#FAA43A"
  
  # Calculate density to find maxima
  dens_at <- density(df$log2FC[df$region_type == "At Anchors"])
  dens_between <- density(df$log2FC[df$region_type == "Between Anchors"])
  
  max_at_x <- dens_at$x[which.max(dens_at$y)]
  max_between_x <- dens_between$x[which.max(dens_between$y)]
  max_at_y <- max(dens_at$y)
  max_between_y <- max(dens_between$y)
  
  # Position labels based on protein
  if (protein_name == "CTCF") {
    at_label_x <- max_at_x + 0.8
    at_label_y <- max_at_y * 0.8
    between_label_x <- -4.2
    between_label_y <- 0.4
  } else if (protein_name == "RAD21") {
    at_label_x <- max_at_x + 0.8
    at_label_y <- max_at_y * 0.8
    between_label_x <- -4.2
    between_label_y <- 0.45  # Moved up from 0.4
  } else {  # YAP1
    at_label_x <- max_at_x + 0.9
    at_label_y <- max_at_y * 0.8
    between_label_x <- max_between_x + 0.9
    between_label_y <- max_between_y * 0.75
  }
  
  ggplot(df, aes(x=log2FC, fill=region_type)) +
    geom_vline(xintercept=0, linetype="dashed", color="gray75", linewidth=.3) +  # Lighter gray
    geom_density(alpha=.4, color=NA) +
    scale_fill_manual(values=c("At Anchors"=at_anchors_color, "Between Anchors"=between_anchors_color), name="") +
    labs(x=paste0(protein_name, " log2 (sorbitol/control)"), y="Density") +
    theme_classic() +
    theme(plot.title=element_text(hjust=.5, face="bold", size=9),
          legend.position="none",
          axis.text=element_text(size=7),
          axis.title=element_text(size=8.5),
          axis.line=element_line(linewidth=.3),
          axis.ticks=element_line(linewidth=.3),
          plot.margin=margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio=1) +
    coord_cartesian(xlim=c(-4.5,4.5)) +
    scale_x_continuous(limits=c(-4.5,4.5)) +
    # Add direct labels with two lines - placed AFTER geom_density so they appear on top
    annotate("text", x=at_label_x, y=at_label_y, label="At\nAnchors", 
             color=at_anchors_color, fontface="bold", size=7/.pt, vjust=0.5, hjust=0, lineheight=0.8) +
    annotate("text", x=between_label_x, y=between_label_y, label="Between\nAnchors", 
             color=between_anchors_color, fontface="bold", size=7/.pt, vjust=0.5, hjust=0, lineheight=0.8)
}

# Panel positioning -------------------------------------------------------

panel_width <- 1.8
panel_height <- 1.8
panel_spacing <- 0.2

ma_panel_width  <- panel_width * 1.07
ma_panel_height <- panel_height * 1.07
ma_dx <- -0.08
ma_dy <- -0.06

density_panel_width <- 0.13
density_panel_height <- ma_panel_height * 0.77
density_spacing <- -0.05
density_y_offset <- ma_panel_height * 0.03

x_start <- 0.25
y_start_row1 <- 0.5
y_start_row2 <- y_start_row1 + panel_height + 0.65
y_start_row3 <- y_start_row2 + panel_height + 0.65

x_col1 <- x_start
x_col2 <- x_col1 + panel_width + panel_spacing + 0.05  # Moved slightly right
x_col3 <- x_col2 + panel_width + panel_spacing

survey_width <- 3.5

x_survey <- x_col3 + panel_width + 0.35

## Survey stack dims - adjust to align genome label bottom with Panel A x-axis label
# Target: genome label bottom should align with Panel A x-axis label
target_genome_label_bottom <- y_start_row3 + ma_panel_height

hic_height    <- 1.60
gene_height   <- 0.30
label_height  <- 0.10
spacing_hic   <- 0.03
spacing_after_hic_to_signal <- 0.10

genes_top_target <- target_genome_label_bottom - gene_height - 0.05 - label_height

spacing_before_genes <- 0.025

signal_height <- 0.35
peak_track_height <- 0.04
signal_gap <- 0.015
spacing_between_proteins <- 0.025
label_to_signal_gap <- 0.005

# Panels A–C --------------------------------------------------------------

plotText(label="A", x=x_col1-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
ctcf_ma_data <- create_ma_data(diff_CTCF, "CTCF")
plotGG(create_ma_plot(ctcf_ma_data, "CTCF"), x=x_col1+ma_dx, y=y_start_row1+ma_dy,
       width=ma_panel_width, height=ma_panel_height)
plotGG(create_ma_density_plot(ctcf_ma_data),
       x=(x_col1+ma_dx)+ma_panel_width+density_spacing,
       y=(y_start_row1+ma_dy)+density_y_offset,
       width=density_panel_width, height=density_panel_height)

plotText(label="B", x=x_col2-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
rad21_ma_data <- create_ma_data(diff_RAD21, "RAD21")
plotGG(create_ma_plot(rad21_ma_data, "RAD21"), x=x_col1+ma_dx, y=y_start_row2+ma_dy,
       width=ma_panel_width, height=ma_panel_height)
plotGG(create_ma_density_plot(rad21_ma_data),
       x=(x_col1+ma_dx)+ma_panel_width+density_spacing,
       y=(y_start_row2+ma_dy)+density_y_offset,
       width=density_panel_width, height=density_panel_height)

yap1_ma_data <- create_ma_data(diff_YAP1, "YAP1")
plotGG(create_ma_plot(yap1_ma_data, "YAP1"), x=x_col1+ma_dx, y=y_start_row3+ma_dy,
       width=ma_panel_width, height=ma_panel_height)
plotGG(create_ma_density_plot(yap1_ma_data),
       x=(x_col1+ma_dx)+ma_panel_width+density_spacing,
       y=(y_start_row3+ma_dy)+density_y_offset,
       width=density_panel_width, height=density_panel_height)

# Bar plots - use exact same positioning as MA plots
ctcf_bar_plot  <- create_bar_plot(bar_plot_data |> filter(Target=="CTCF"))
rad21_bar_plot <- create_bar_plot(bar_plot_data |> filter(Target=="RAD21"))
yap1_bar_plot  <- create_bar_plot(bar_plot_data |> filter(Target=="YAP1"))
plotGG(ctcf_bar_plot,  x=x_col2+ma_dx, y=y_start_row1+ma_dy, width=ma_panel_width, height=ma_panel_height)
plotGG(rad21_bar_plot, x=x_col2+ma_dx, y=y_start_row2+ma_dy, width=ma_panel_width, height=ma_panel_height)
plotGG(yap1_bar_plot,  x=x_col2+ma_dx, y=y_start_row3+ma_dy, width=ma_panel_width, height=ma_panel_height)

plotText(label="C", x=x_col3-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
# Density plots - use exact same positioning as MA plots
plotGG(create_density_plot(density_data |> filter(protein=="CTCF", loop_category=="Gained"), "CTCF"),
       x=x_col3+ma_dx, y=y_start_row1+ma_dy, width=ma_panel_width, height=ma_panel_height)
plotGG(create_density_plot(density_data |> filter(protein=="RAD21", loop_category=="Gained"), "RAD21"),
       x=x_col3+ma_dx, y=y_start_row2+ma_dy, width=ma_panel_width, height=ma_panel_height)
plotGG(create_density_plot(density_data |> filter(protein=="YAP1", loop_category=="Gained"), "YAP1"),
       x=x_col3+ma_dx, y=y_start_row3+ma_dy, width=ma_panel_width, height=ma_panel_height)

## row labels - vertical on left side
ctcf_row_mid <- y_start_row1 + panel_height/2
rad21_row_mid <- y_start_row2 + panel_height/2
yap1_row_mid <- y_start_row3 + panel_height/2

plotText("CTCF", x=x_col1 - 0.25, y=ctcf_row_mid,
         fontcolor=gray_color, fontsize=11, fontface="bold",
         just=c("right","center"), rot=0)
plotText("RAD21", x=x_col1 - 0.25, y=rad21_row_mid,
         fontcolor=gray_color, fontsize=11, fontface="bold",
         just=c("right","center"), rot=0)
plotText("YAP1", x=x_col1 - 0.25, y=yap1_row_mid,
         fontcolor=gray_color, fontsize=11, fontface="bold",
         just=c("right","center"), rot=0)


# Panel D -----------------------------------------------------------------

plotText(label="D", x=x_survey-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")

surveyParams <- pgParams(
  assembly="hg38", resolution=10e3,
  chrom = as.character(seqnames(loopRegions_gained_buffed))[105],
  chromstart = start(loopRegions_gained_buffed)[105],
  chromend   = end(loopRegions_gained_buffed)[105],
  zrange = c(0,100), norm = "SCALE"
)

current_y <- y_start_row1 + 0.03
highlight_y_start <- current_y

survey_control <- plotHicRectangle(
  data="data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic",
  params=surveyParams, x=x_survey, y=current_y, width=survey_width, height=hic_height
)

# Add heatmap legend to the right of control Hi-C map
annoHeatmapLegend(survey_control,
                  params=surveyParams,
                  orientation="v",
                  fontcolor="black",
                  digits=2,
                  x=x_survey + survey_width + 0.05,
                  y=current_y,
                  width=0.05,
                  height=hic_height/2,  # Half the height of Hi-C map
                  just=c("left","top"),
                  default.units="inches")

plotText("control", x=x_survey+0.05, y=current_y+0.05, fontsize=5,
         fontcolor="black", just=c("left","top"))
annoPixels(survey_control, data=lostLoops,  type="arrow", col="#619CFF", lwd=.001, lty=.001, shift=4)
annoPixels(survey_control, data=gainedLoops, type="arrow", col="#F8766D", lwd=.001, lty=.001, shift=4)

current_y <- current_y + hic_height + spacing_hic

survey_sorb <- plotHicRectangle(
  data="data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic",
  params=surveyParams, x=x_survey, y=current_y, width=survey_width, height=hic_height
)
plotText("200mM sorbitol", x=x_survey+0.05, y=current_y+0.05, fontsize=5,
         fontcolor="black", just=c("left","top"))
annoPixels(survey_sorb, data=lostLoops,  type="arrow", col="#619CFF", lwd=.001, lty=.001, shift=4)
annoPixels(survey_sorb, data=gainedLoops, type="arrow", col="#F8766D", lwd=.001, lty=.001, shift=4)

current_y <- current_y + hic_height + spacing_hic

## Signal ranges
ctcf_signalRange <- calcSignalRange(c(ctcf_control_bw, ctcf_sorb_bw),
                                    chrom=surveyParams$chrom,
                                    chromstart=surveyParams$chromstart,
                                    chromend=surveyParams$chromend,
                                    assembly="hg38", negData=FALSE)
rad21_signalRange <- calcSignalRange(c(rad21_control_bw, rad21_sorb_bw),
                                     chrom=surveyParams$chrom,
                                     chromstart=surveyParams$chromstart,
                                     chromend=surveyParams$chromend,
                                     assembly="hg38", negData=FALSE)
yap1_signalRange <- calcSignalRange(c(yap1_control_bw, yap1_sorb_bw),
                                    chrom=surveyParams$chrom,
                                    chromstart=surveyParams$chromstart,
                                    chromend=surveyParams$chromend,
                                    assembly="hg38", negData=FALSE)

retained_peak_color <- "#4DAF4A"
lost_peak_color <- "#FF7F0E"

# Add space between Hi-C map and first CTCF track
current_y <- current_y + spacing_after_hic_to_signal

ctcf_tracks_y_start <- current_y

# CTCF control signal
plotSignal(ctcf_control_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=ctcf_color, linecolor=ctcf_color, range=ctcf_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# CTCF control peaks - show retained (green) and lost (red)
plotRanges(
  data = ctcf_changes$retained_control, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)
plotRanges(
  data = ctcf_changes$lost, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=lost_peak_color, linecolor=lost_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add range label and "control" label at same y-position below peaks
ctcf_max <- ceiling(ctcf_signalRange[2])
range_label_x <- x_survey + survey_width
plotText(paste0("0-", ctcf_max), x=range_label_x, y=current_y, 
         fontsize=5, fontcolor=gray_color, just=c("right","top"))
plotText("control", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y + 0.03

# CTCF sorbitol signal
plotSignal(ctcf_sorb_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=ctcf_color, linecolor=ctcf_color, range=ctcf_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# CTCF sorbitol peaks - show retained only (green)
plotRanges(
  data = ctcf_changes$retained_sorbitol, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add "+ sorbitol" label at same y-position as range labels would be
plotText("+ sorbitol", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y

ctcf_tracks_y_end <- current_y
ctcf_mid <- ctcf_tracks_y_start + (ctcf_tracks_y_end - ctcf_tracks_y_start)/2

plotText("CTCF", x=x_survey - 0.12, y=ctcf_mid,
         fontsize=9, fontface="bold", fontcolor=ctcf_color,
         just=c("right","center"), rot=90)

current_y <- current_y + spacing_between_proteins

rad21_tracks_y_start <- current_y

# RAD21 control signal
plotSignal(rad21_control_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=rad21_color, linecolor=rad21_color, range=rad21_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# RAD21 control peaks - show retained (green) and lost (red)
plotRanges(
  data = rad21_changes$retained_control, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)
plotRanges(
  data = rad21_changes$lost, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=lost_peak_color, linecolor=lost_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add range label and "control" label at same y-position
rad21_max <- ceiling(rad21_signalRange[2])
plotText(paste0("0-", rad21_max), x=range_label_x, y=current_y,
         fontsize=5, fontcolor=gray_color, just=c("right","top"))
plotText("control", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y + 0.03

# RAD21 sorbitol signal
plotSignal(rad21_sorb_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=rad21_color, linecolor=rad21_color, range=rad21_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# RAD21 sorbitol peaks - show retained only (green)
plotRanges(
  data = rad21_changes$retained_sorbitol, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add "+ sorbitol" label
plotText("+ sorbitol", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y

rad21_tracks_y_end <- current_y
rad21_mid <- rad21_tracks_y_start + (rad21_tracks_y_end - rad21_tracks_y_start)/2

plotText("RAD21", x=x_survey - 0.12, y=rad21_mid,
         fontsize=9, fontface="bold", fontcolor=rad21_color,
         just=c("right","center"), rot=90)

current_y <- current_y + spacing_between_proteins

# YAP1 tracks
yap1_tracks_y_start <- current_y

# YAP1 control signal
plotSignal(yap1_control_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=yap1_color, linecolor=yap1_color, range=yap1_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# YAP1 control peaks - show retained (green) and lost (red)
plotRanges(
  data = yap1_changes$retained_control, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)
plotRanges(
  data = yap1_changes$lost, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=lost_peak_color, linecolor=lost_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add range label and "control" label at same y-position
yap1_max <- ceiling(yap1_signalRange[2])
plotText(paste0("0-", yap1_max), x=range_label_x, y=current_y,
         fontsize=5, fontcolor=gray_color, just=c("right","top"))
plotText("control", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y + 0.03

# YAP1 sorbitol signal
plotSignal(yap1_sorb_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill=yap1_color, linecolor=yap1_color, range=yap1_signalRange, scale=FALSE)

current_y <- current_y + signal_height + label_to_signal_gap

current_y <- current_y + signal_gap

# YAP1 sorbitol peaks - show retained only (green)
plotRanges(
  data = yap1_changes$retained_sorbitol, params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + 0.02

# Add "+ sorbitol" label
plotText("+ sorbitol", x=x_survey, y=current_y, fontsize=6, 
         fontcolor=gray_color, just=c("left","top"))

current_y <- current_y

yap1_tracks_y_end <- current_y
yap1_mid <- yap1_tracks_y_start + (yap1_tracks_y_end - yap1_tracks_y_start)/2

plotText("YAP1", x=x_survey - 0.12, y=yap1_mid,
         fontsize=9, fontface="bold", fontcolor=yap1_color,
         just=c("right","center"), rot=90)

current_y <- current_y + spacing_before_genes

genes_y_start <- current_y

plotGenes(param=surveyParams, chrom=surveyParams$chrom,
          x=x_survey, y=current_y, width=survey_width, height=gene_height, fontsize=7)

current_y <- current_y + gene_height + 0.05

plotGenomeLabel(params=surveyParams, x=x_survey, y=current_y,
                length=survey_width, scale="bp", fontsize=7)

# Set highlight to end at the top of genome label (where genes end)
highlight_y_end <- current_y

## Highlights for loop #105
loop_105_anchor1 <- anchors(gainedLoops, "first")[105]
loop_105_anchor2 <- anchors(gainedLoops, "second")[105]
total_highlight_height <- highlight_y_end - highlight_y_start

buffer_10kb <- 10000
anchor1_start_expanded <- max(1, start(loop_105_anchor1) - buffer_10kb)
anchor1_end_expanded <- end(loop_105_anchor1) + buffer_10kb
anchor2_start_expanded <- max(1, start(loop_105_anchor2) - buffer_10kb)
anchor2_end_expanded <- end(loop_105_anchor2) + buffer_10kb

annoHighlight(plot=survey_control,
              chrom=as.character(seqnames(loop_105_anchor1)),
              chromstart=anchor1_start_expanded, chromend=anchor1_end_expanded,
              y=highlight_y_start, height=total_highlight_height,
              just=c("left","top"), fill="grey85", alpha=0.4)
annoHighlight(plot=survey_control,
              chrom=as.character(seqnames(loop_105_anchor2)),
              chromstart=anchor2_start_expanded, chromend=anchor2_end_expanded,
              y=highlight_y_start, height=total_highlight_height,
              just=c("left","top"), fill="grey85", alpha=0.4)

anchor1_center_x <- x_survey + ((anchor1_start_expanded + anchor1_end_expanded) / 2 - 
                                  surveyParams$chromstart) / 
  (surveyParams$chromend - surveyParams$chromstart) * survey_width
anchor2_center_x <- x_survey + ((anchor2_start_expanded + anchor2_end_expanded) / 2 - 
                                  surveyParams$chromstart) / 
  (surveyParams$chromend - surveyParams$chromstart) * survey_width

text_y_position <- ctcf_tracks_y_start - 0.05
plotText("Retained Peak", x=anchor1_center_x, y=text_y_position,
         fontsize=7, fontcolor=retained_peak_color, just=c("center","center"))
plotText("Retained Peak", x=anchor2_center_x, y=text_y_position,
         fontsize=7, fontcolor=retained_peak_color, just=c("center","center"))

middle_x <- x_survey + survey_width / 2
lost_peaks_y <- ctcf_tracks_y_start - 0.05
plotText("Lost peaks", x=middle_x, y=lost_peaks_y,
         fontsize=7, fontcolor=lost_peak_color, just=c("center","center"))

dev.off()
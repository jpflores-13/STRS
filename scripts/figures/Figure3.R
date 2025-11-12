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
page_height <- 8.5

pdf("figures/Figure3.pdf", width = page_width, height = page_height)
pageCreate(width = page_width, height = page_height, showGuides = FALSE)

gray_color <- "#666666"

# Load Data ---------------------------------------------------------------

diff_CTCF <- readRDS("data/processed/cutntag/deseq2/diff_CTCF_counts.rds")
diff_RAD21 <- readRDS("data/processed/cutntag/deseq2/diff_RAD21_counts.rds")

loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |>
  interactions() |> as.data.frame() |> as_ginteractions()

ctcf_control_bw <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_cont_0h.bw"
ctcf_sorb_bw    <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h.bw"
rad21_control_bw<- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_cont_0h.bw"
rad21_sorb_bw   <- "data/processed/cutntag/output/mergeSignal/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h.bw"

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
for (protein in c("CTCF","RAD21")) {
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

# Plot helpers ------------------------------------------------------------

create_ma_data <- function(diff_obj, protein_name) {
  as.data.frame(mcols(diff_obj)) |>
    dplyr::select(baseMean, log2FoldChange, padj) |>
    mutate(isDE = case_when(
      log2FoldChange > 0 & padj < 0.1 ~ "Increased",
      log2FoldChange < 0 & padj < 0.1 ~ "Decreased",
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
    labs(y="log2FC(sorbitol/control)", x="mean of normalized counts") +
    theme_classic() +
    theme(legend.position="none", plot.title=element_text(hjust=.5, face="bold", size=9),
          axis.text=element_text(size=7.5), axis.title=element_text(size=8.5),
          axis.title.y=element_text(angle=90, vjust=.5), aspect.ratio=1) +
    annotate("text", label="Increased", x=max(res_df$baseMean,na.rm=TRUE)*.1, y=3.6, color="#F8766D", fontface="bold", size=7/.pt) +
    annotate("text", label=paste0("n = ", n_up), x=max(res_df$baseMean,na.rm=TRUE)*.1, y=3.1, color="#F8766D", size=6/.pt) +
    annotate("text", label="Decreased", x=max(res_df$baseMean,na.rm=TRUE)*.1, y=-3.4, color="#619CFF", fontface="bold", size=7/.pt) +
    annotate("text", label=paste0("n = ", n_dn), x=max(res_df$baseMean,na.rm=TRUE)*.1, y=-3.9, color="#619CFF", size=6/.pt)
}
create_ma_density_plot <- function(res_df) {
  ggplot(res_df, aes(y=log2FoldChange)) +
    geom_density(color="#9370DB", fill="#9370DB", alpha=.25) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=.3) +
    ylim(c(-4,4)) + xlim(c(0,1)) +
    theme_classic() +
    theme(legend.position="none", axis.text=element_blank(), axis.title=element_blank(),
          axis.ticks=element_blank(), axis.line=element_blank(),
          plot.margin=margin(0,0,0,0,"pt"), panel.spacing=unit(0,"pt"))
}
create_bar_plot <- function(d) {
  ggplot(d, aes(x=Category, y=Percentage, fill=Condition)) +
    geom_hline(yintercept=seq(0,100,25), color="gray90", linetype="dashed") +
    geom_bar(stat="identity", position=position_dodge(.7), width=.6) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), position=position_dodge(.7), width=.25, linewidth=.3) +
    geom_text(aes(label=sprintf("%.1f%%", Percentage)), position=position_dodge(.7), vjust=-.5, size=6/.pt) +
    scale_fill_manual(values=c(control="#619CFF", sorbitol="#F8766D")) +
    theme_bw() +
    theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_blank(), axis.title.x=element_blank(),
          legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank(),
          legend.key.size=unit(.3,"cm"), legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,0,0,0), axis.text=element_text(size=7),
          axis.title=element_text(size=8), aspect.ratio=1) +
    labs(y="% of anchors bound") +
    scale_y_continuous(limits=c(0,100), expand=expansion(mult=c(0,.05)), breaks=seq(0,100,25))
}
create_density_plot <- function(df) {
  df <- df |> mutate(region_type=factor(region_type, levels=c("At Anchors","Between Anchors")))
  ggplot(df, aes(x=log2FC, fill=region_type, color=region_type)) +
    geom_vline(xintercept=0, linetype="dashed", color="gray30", linewidth=.3) +
    geom_density(alpha=.4, linewidth=.8) +
    scale_fill_manual(values=c("At Anchors"="#5DA5DA","Between Anchors"="#FAA43A"), name="") +
    scale_color_manual(values=c("At Anchors"="#5DA5DA","Between Anchors"="#FAA43A"), name="") +
    labs(title="Gained Loops", x="log2FC(sorbitol/control)", y="Density") +
    theme_classic() +
    theme(plot.title=element_text(hjust=.5, face="bold", size=9),
          legend.position="bottom", legend.text=element_text(size=7),
          legend.key.size=unit(.3,"cm"), legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-5,0,0,0), axis.text=element_text(size=7),
          axis.title=element_text(size=8), aspect.ratio=1) +
    coord_cartesian(xlim=c(-3,3))
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

x_col1 <- x_start
x_col2 <- x_col1 + panel_width + panel_spacing
x_col3 <- x_col2 + panel_width + panel_spacing

survey_width <- 3.5

x_survey <- x_col3 + panel_width + 0.35
apa_x <- x_survey + survey_width - 0.05

apa_bottom_target <- y_start_row2 + panel_height
apa_height <- 0.65
apa_label_space <- 0.1
apa_y_start <- apa_bottom_target - apa_height

## Survey stack dims
hic_height    <- 1.20
signal_height <- 0.18
gene_height   <- 0.28
label_height  <- 0.08

signal_gap    <- 0.05
spacing       <- 0.04

# Panels A–C --------------------------------------------------------------

plotText(label="A", x=x_col1-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
ctcf_ma_data <- create_ma_data(diff_CTCF, "CTCF")
plotGG(create_ma_plot(ctcf_ma_data, "CTCF"), x=x_col1+ma_dx, y=y_start_row1+ma_dy,
       width=ma_panel_width, height=ma_panel_height)
plotGG(create_ma_density_plot(ctcf_ma_data),
       x=(x_col1+ma_dx)+ma_panel_width+density_spacing,
       y=(y_start_row1+ma_dy)+density_y_offset,
       width=density_panel_width, height=density_panel_height)

rad21_ma_data <- create_ma_data(diff_RAD21, "RAD21")
plotText(label="B", x=x_col2-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
plotGG(create_ma_plot(rad21_ma_data, "RAD21"), x=x_col1+ma_dx, y=y_start_row2+ma_dy,
       width=ma_panel_width, height=ma_panel_height)
plotGG(create_ma_density_plot(rad21_ma_data),
       x=(x_col1+ma_dx)+ma_panel_width+density_spacing,
       y=(y_start_row2+ma_dy)+density_y_offset,
       width=density_panel_width, height=density_panel_height)

ctcf_bar_plot  <- create_bar_plot(bar_plot_data |> filter(Target=="CTCF"))
rad21_bar_plot <- create_bar_plot(bar_plot_data |> filter(Target=="RAD21"))
plotGG(ctcf_bar_plot,  x=x_col2, y=y_start_row1, width=panel_width, height=panel_height)
plotGG(rad21_bar_plot, x=x_col2, y=y_start_row2, width=panel_width, height=panel_height)

plotText(label="C", x=x_col3-0.15, y=y_start_row1-0.2, fontsize=12, fontface="bold")
plotGG(create_density_plot(density_data |> filter(protein=="CTCF", loop_category=="Gained")),
       x=x_col3, y=y_start_row1, width=panel_width, height=panel_height)
plotGG(create_density_plot(density_data |> filter(protein=="RAD21", loop_category=="Gained")),
       x=x_col3, y=y_start_row2, width=panel_width, height=panel_height)

## row labels/segments
row1_y <- y_start_row1 + panel_height + 0.05
row2_y <- y_start_row2 + panel_height + 0.05
plotSegments(x0=x_col1, y0=row1_y, x1=x_col3+panel_width, y1=row1_y,
             default.units="inches", linecolor=gray_color, lwd=1)
plotText("CTCF", x=(x_col1 + (x_col3+panel_width))/2, y=row1_y+0.1,
         fontcolor=gray_color, fontsize=10)
plotSegments(x0=x_col1, y0=row2_y, x1=x_col3+panel_width, y1=row2_y,
             default.units="inches", linecolor=gray_color, lwd=1)
plotText("RAD21", x=(x_col1 + (x_col3+panel_width))/2, y=row2_y+0.1,
         fontcolor=gray_color, fontsize=10)

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
plotText("untreated", x=x_survey+0.05, y=current_y+0.05, fontsize=5,
         fontcolor="black", just=c("left","top"))
annoPixels(survey_control, data=lostLoops,  type="arrow", col="#619CFF", lwd=.001, lty=.001, shift=4)
annoPixels(survey_control, data=gainedLoops, type="arrow", col="#F8766D", lwd=.001, lty=.001, shift=4)

current_y <- current_y + hic_height + spacing

survey_sorb <- plotHicRectangle(
  data="data/processed/hic/maps/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic",
  params=surveyParams, x=x_survey, y=current_y, width=survey_width, height=hic_height
)
plotText("200mM sorbitol", x=x_survey+0.05, y=current_y+0.05, fontsize=5,
         fontcolor="black", just=c("left","top"))
annoPixels(survey_sorb, data=lostLoops,  type="arrow", col="#619CFF", lwd=.001, lty=.001, shift=4)
annoPixels(survey_sorb, data=gainedLoops, type="arrow", col="#F8766D", lwd=.001, lty=.001, shift=4)

current_y <- current_y + hic_height + spacing

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

retained_peak_color <- "#4DAF4A"
lost_peak_color <- "#FF7F0E"
peak_track_height <- 0.08

# Add extra space between Hi-C map and first CTCF track for text labels
current_y <- current_y + 0.15

ctcf_tracks_y_start <- current_y

plotSignal(ctcf_control_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill="#253494", linecolor="#253494", range=ctcf_signalRange, scale=TRUE)
plotText("untreated", x=x_survey+survey_width, y=current_y - 0.15,
         fontsize=7, fontcolor=gray_color, just=c("right","top"))

current_y <- current_y + signal_height + 0.01

plotRanges(
  data = peak_list[["CTCF_control"]], params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + signal_gap

plotSignal(ctcf_sorb_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill="#253494", linecolor="#253494", range=ctcf_signalRange, scale=TRUE)
plotText("+ sorbitol", x=x_survey+survey_width, y=current_y - 0.15,
         fontsize=7, fontcolor=gray_color, just=c("right","top"))

current_y <- current_y + signal_height + 0.01

plotRanges(
  data = peak_list[["CTCF_sorbitol"]], params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=lost_peak_color, linecolor=lost_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height

ctcf_tracks_y_end <- current_y
ctcf_mid <- ctcf_tracks_y_start + (ctcf_tracks_y_end - ctcf_tracks_y_start)/2

plotText("CTCF", x=x_survey - 0.14, y=ctcf_mid - 0.20,
         fontsize=9, fontface="bold", fontcolor=gray_color,
         just=c("right","center"), rot=90)

current_y <- current_y + signal_gap

rad21_tracks_y_start <- current_y

plotSignal(rad21_control_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill="#41B6C4", linecolor="#41B6C4", range=rad21_signalRange, scale=TRUE)
plotText("untreated", x=x_survey+survey_width, y=current_y - 0.15,
         fontsize=7, fontcolor=gray_color, just=c("right","top"))

current_y <- current_y + signal_height + 0.01

plotRanges(
  data = peak_list[["RAD21_control"]], params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=retained_peak_color, linecolor=retained_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height + signal_gap

plotSignal(rad21_sorb_bw, params=surveyParams,
           x=x_survey, y=current_y, width=survey_width, height=signal_height,
           fill="#41B6C4", linecolor="#41B6C4", range=rad21_signalRange, scale=TRUE)
plotText("+ sorbitol", x=x_survey+survey_width, y=current_y - 0.15,
         fontsize=7, fontcolor=gray_color, just=c("right","top"))

current_y <- current_y + signal_height + 0.01

plotRanges(
  data = peak_list[["RAD21_sorbitol"]], params=surveyParams,
  x=x_survey, y=current_y, width=survey_width, height=peak_track_height,
  fill=lost_peak_color, linecolor=lost_peak_color, collapse=TRUE,
  stroke=0.08, strokecolor="black",
  just=c("left","top"), default.units="inches"
)

current_y <- current_y + peak_track_height

rad21_tracks_y_end <- current_y
rad21_mid <- rad21_tracks_y_start + (rad21_tracks_y_end - rad21_tracks_y_start)/2

plotText("RAD21", x=x_survey - 0.14, y=rad21_mid - 0.20,
         fontsize=9, fontface="bold", fontcolor=gray_color,
         just=c("right","center"), rot=90)

current_y <- current_y + 0.12

plotGenes(param=surveyParams, chrom=surveyParams$chrom,
          x=x_survey, y=current_y, width=survey_width, height=gene_height, fontsize=7)

current_y <- current_y + gene_height + 0.03

plotGenomeLabel(params=surveyParams, x=x_survey, y=current_y,
                length=survey_width, scale="bp", fontsize=7)

highlight_y_end <- current_y + label_height

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

text_y_position <- ctcf_tracks_y_start - 0.08
plotText("Retained Peak", x=anchor1_center_x, y=text_y_position,
         fontsize=7, fontcolor=retained_peak_color, just=c("center","center"))
plotText("Retained Peak", x=anchor2_center_x, y=text_y_position,
         fontsize=7, fontcolor=retained_peak_color, just=c("center","center"))

middle_x <- x_survey + survey_width / 2
lost_peaks_y <- ctcf_tracks_y_start - 0.08
plotText("Lost peaks", x=middle_x, y=lost_peaks_y,
         fontsize=7, fontcolor=lost_peak_color, just=c("center","center"))

dev.off()

cat("\nFigure created successfully!\n")
cat("Saved to: figures/Figure3.pdf\n")

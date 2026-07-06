# ##############################################################################
# filename:    FigureSX.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Supplementary Figure SX — CUT&Tag anchor vs between-anchor
#              density analysis; CTCF and RAD21 log2FC distributions for
#              gained, static, and lost loop categories
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(bamsignals)
library(GenomicRanges)
library(plotgardener)

# Parameters ----
diff_loops_rds   <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
cutntag_bam_dir   <- "data/processed/cutntag/output/mergeAlign/"
output_pdf        <- "figures/FigureSX.pdf"
page_width        <- 8.5
page_height       <- 8

# Helper functions ----
get_loop_anchors <- function(loops) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  all_anchors    <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_id <- c(
    paste0(seq_along(first_anchors),  "_anchor1"),
    paste0(seq_along(second_anchors), "_anchor2")
  )
  all_anchors
}

get_between_regions <- function(loops) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  starts     <- end(first_anchors) + 1
  ends       <- start(second_anchors) - 1
  valid_idx  <- ends >= starts
  if (sum(valid_idx) == 0) return(GRanges())
  GRanges(
    seqnames = seqnames(first_anchors)[valid_idx],
    ranges   = IRanges(start = starts[valid_idx], end = ends[valid_idx]),
    loop_id  = paste0(which(valid_idx), "_between")
  )
}

analyze_regions <- function(regions, peaks_control, peaks_treat,
                             bam_control, bam_treat, region_type, loop_category) {
  overlaps <- findOverlaps(peaks_control, regions)
  if (length(overlaps) == 0) return(data.frame())

  overlapping_peaks <- peaks_control[queryHits(overlaps)]
  mcols(overlapping_peaks)$region_id <- mcols(regions)$loop_id[subjectHits(overlaps)]

  peak_ids <- paste0(as.character(seqnames(overlapping_peaks)), ":",
                     start(overlapping_peaks), "-", end(overlapping_peaks))

  control_signal <- bamCount(bam_control, overlapping_peaks, paired.end = "midpoint")
  treat_signal   <- bamCount(bam_treat,   overlapping_peaks, paired.end = "midpoint")

  log2FC <- case_when(
    control_signal == 0 & treat_signal == 0 ~ 0,
    control_signal == 0                     ~ log2(treat_signal + 1),
    treat_signal   == 0                     ~ -log2(control_signal + 1),
    TRUE                                    ~ log2((treat_signal + 1) / (control_signal + 1))
  )

  data.frame(
    loop_category  = loop_category,
    region_type    = region_type,
    region_id      = mcols(overlapping_peaks)$region_id,
    peak_id        = peak_ids,
    log2FC         = log2FC,
    control_signal = control_signal,
    treat_signal   = treat_signal,
    peak_score     = mcols(overlapping_peaks)$signalValue,
    stringsAsFactors = FALSE
  ) |>
    group_by(peak_id) |>
    slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    ungroup()
}

create_density_plot_with_labels <- function(data, protein_name) {
  at_anchors_color      <- "#5DA5DA"
  between_anchors_color <- "#FAA43A"

  ggplot(data, aes(x = log2FC, fill = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray75", linewidth = 0.3) +
    geom_density(alpha = 0.4, color = NA) +
    facet_wrap(~loop_category, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("At Anchors" = at_anchors_color,
                 "Between Anchors" = between_anchors_color),
      name = ""
    ) +
    labs(x = paste0(protein_name, " log2 (sorbitol/control)"),
         y = "Density") +
    annotate("text", x = 0, y = Inf, label = "At\nAnchors",
             color = at_anchors_color, fontface = "bold",
             size = 7 / .pt, lineheight = 0.8, vjust = 1.5) +
    annotate("text", x = -3, y = Inf, label = "Between\nAnchors",
             color = between_anchors_color, fontface = "bold",
             size = 7 / .pt, lineheight = 0.8, vjust = 1.5) +
    theme_classic() +
    theme(
      strip.text         = element_text(face = "bold", size = 8),
      strip.background   = element_blank(),
      legend.position    = "none",
      panel.spacing      = unit(0.5, "lines"),
      axis.line          = element_line(linewidth = 0.3),
      axis.ticks         = element_line(linewidth = 0.3),
      axis.title         = element_text(size = 8.5),
      axis.text          = element_text(size = 7),
      plot.margin        = margin(5.5, 5.5, 5.5, 5.5, "pt")
    ) +
    coord_cartesian(xlim = c(-4.5, 4.5), clip = "off") +
    scale_x_continuous(limits = c(-4.5, 4.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

# Data import ----
loops <- readRDS(diff_loops_rds) |> interactions()

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]
staticLoops <- loops[loops$padj > 0.1]

target    <- c("CTCF", "RAD21")
condition <- c("control", "sorbitol")

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern   <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = pattern)
    file_path <- file_path[grepl("\\.narrowPeak$", file_path)]
    if (length(file_path) == 1)
      peak_list[[paste0(t, "_", cond)]] <- read_narrowpeaks(file_path)
  }
}

bam_files <- character()
for (t in target) {
  for (cond in condition) {
    pattern   <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(cutntag_bam_dir, full.names = TRUE, pattern = pattern)
    file_path <- file_path[grepl("\\.bam$", file_path)]
    if (length(file_path) == 1)
      bam_files[paste0(t, "_", cond)] <- file_path
  }
}

# Analysis ----
loop_categories <- list(Gained = gainedLoops, Static = staticLoops, Lost = lostLoops)
results_list    <- list()

for (protein in target) {
  control_peaks  <- peak_list[[paste0(protein, "_control")]]
  sorbitol_peaks <- peak_list[[paste0(protein, "_sorbitol")]]
  control_bam    <- bam_files[paste0(protein, "_control")]
  sorbitol_bam   <- bam_files[paste0(protein, "_sorbitol")]

  for (category in names(loop_categories)) {
    loops_subset <- loop_categories[[category]]
    anch         <- get_loop_anchors(loops_subset)
    between      <- get_between_regions(loops_subset)

    anchor_results  <- analyze_regions(anch,    control_peaks, sorbitol_peaks,
                                       control_bam, sorbitol_bam, "At Anchors", category)
    between_results <- analyze_regions(between, control_peaks, sorbitol_peaks,
                                       control_bam, sorbitol_bam, "Between Anchors", category)

    results_list[[paste0(protein, "_", category)]] <- rbind(anchor_results, between_results) |>
      mutate(protein = protein)
  }
}

all_results <- bind_rows(results_list)

ctcf_data <- all_results |>
  filter(protein == "CTCF") |>
  mutate(loop_category = factor(loop_category, levels = c("Gained", "Static", "Lost")),
         region_type   = factor(region_type,   levels = c("At Anchors", "Between Anchors")))

rad21_data <- all_results |>
  filter(protein == "RAD21") |>
  mutate(loop_category = factor(loop_category, levels = c("Gained", "Static", "Lost")),
         region_type   = factor(region_type,   levels = c("At Anchors", "Between Anchors")))

ctcf_plot  <- create_density_plot_with_labels(ctcf_data,  "CTCF")
rad21_plot <- create_density_plot_with_labels(rad21_data, "RAD21")

# Visualization ----
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

panel_width   <- 4
panel_height  <- 3.5
panel_spacing <- 0.25
x_col1 <- 0.25
x_col2 <- x_col1 + panel_width + panel_spacing
y_start <- 0.5

plotText(label = "A", fontsize = 12, fontface = "bold",
         x = x_col1 - 0.15, y = y_start - 0.2)
plotGG(plot = ctcf_plot,  x = x_col1, y = y_start,
       width = panel_width, height = panel_height)

plotText(label = "B", fontsize = 12, fontface = "bold",
         x = x_col2 - 0.15, y = y_start - 0.2)
plotGG(plot = rad21_plot, x = x_col2, y = y_start,
       width = panel_width, height = panel_height)

# Save outputs ----
dev.off()

sessionInfo()

# ##############################################################################
# filename:    FigureS9.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Supplementary Figure S9 — H3K27ac differential analysis;
#              MA plot, anchor overlap bar plot, and anchor vs between-anchor
#              density plot for gained and lost loop categories
# ##############################################################################

# Libraries ----
library(DESeq2)
library(InteractionSet)
library(plotgardener)
library(tidyverse)
library(RColorBrewer)
library(mariner)
library(plyranges)
library(ggpubr)
library(bamsignals)

# Parameters ----
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
diff_H3K27ac_rds  <- "data/processed/cutntag/deseq2/diff_H3K27ac_counts.rds"
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
cutntag_bam_dir   <- "data/processed/cutntag/output/mergeAlign/"
output_pdf        <- "figures/FigureS9.pdf"
page_width        <- 7.5
page_height       <- 3.0

# Data import ----
loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

diff_H3K27ac <- readRDS(diff_H3K27ac_rds)

# Analysis ----
gainedLoops <- loops[loops$padj < 0.05 & loops$log2FoldChange > 0] |>
  as.data.frame() |>
  as_ginteractions()
lostLoops <- loops[loops$padj < 0.05 & loops$log2FoldChange < 0] |>
  as.data.frame() |>
  as_ginteractions()

target    <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")

cutntag <- list.files(cutntag_peaks_dir,
                      full.names = TRUE,
                      pattern    = ".narrowPeak") |>
  lapply(read_narrowpeaks)
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

bam_files <- character()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    fp      <- list.files(cutntag_bam_dir, full.names = TRUE, pattern = pattern)
    fp      <- fp[grepl("\\.bam$", fp)]
    if (length(fp) == 1) bam_files[paste0(t, "_", cond)] <- fp
  }
}

## MA plot data ----
create_ma_data <- function(diff_obj, protein_name) {
  as.data.frame(mcols(diff_obj)) |>
    dplyr::select(baseMean, log2FoldChange, padj) |>
    mutate(isDE = case_when(
      log2FoldChange >  1 & padj < 0.05 ~ "Increased",
      log2FoldChange < -1 & padj < 0.05 ~ "Decreased",
      TRUE                               ~ "Not significant")) |>
    arrange(isDE)
}

h3k27ac_ma_data <- create_ma_data(diff_H3K27ac, "H3K27ac")

## Bar plot data ----
extractAnchors <- function(gi) {
  unique(c(anchors(gi, "first"), anchors(gi, "second")))
}

calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  n_anchors  <- length(anchors)
  boot_props <- replicate(n_bootstrap, {
    boot_idx <- sample(n_anchors, replace = TRUE)
    sum(countOverlaps(anchors[boot_idx], peaks) > 0) / n_anchors
  })
  qs <- quantile(boot_props, c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
  data.frame(lower = qs[1] * 100, upper = qs[2] * 100)
}

calculateOverlaps <- function(anchors, peaks_list, category) {
  props <- lapply(peaks_list, function(p) {
    sum(countOverlaps(anchors, p) > 0) / length(anchors)
  })
  cis <- lapply(peaks_list, function(p) calculateCI(anchors, p))
  data.frame(
    Category   = category,
    Target     = sub("_.*", "", names(props)),
    Condition  = sub(".*_", "", names(props)),
    Percentage = unlist(props) * 100,
    Lower      = sapply(cis, \(x) x$lower),
    Upper      = sapply(cis, \(x) x$upper)
  )
}

gained_anchors_bar <- extractAnchors(loops[loops$padj < 0.05 & loops$log2FoldChange > 0])
lost_anchors_bar   <- extractAnchors(loops[loops$padj < 0.05 & loops$log2FoldChange < 0])

bar_results <- rbind(
  calculateOverlaps(gained_anchors_bar, cutntag, "Gained"),
  calculateOverlaps(lost_anchors_bar,   cutntag, "Lost")
)

bar_plot_data <- bar_results |>
  filter(Target == "H3K27ac") |>
  mutate(Category = factor(Category, levels = c("Lost", "Gained")))

## Density analysis data ----
peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    fp      <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = pattern)
    fp      <- fp[grepl("\\.narrowPeak$", fp)]
    if (length(fp) == 1) peak_list[[paste0(t, "_", cond)]] <- read_narrowpeaks(fp)
  }
}

get_between_regions <- function(loops) {
  a1 <- anchors(loops, "first")
  a2 <- anchors(loops, "second")
  st <- end(a1) + 1
  en <- start(a2) - 1
  ok <- en >= st
  if (!any(ok)) return(GRanges())
  GRanges(
    seqnames = seqnames(a1)[ok],
    ranges   = IRanges(start = st[ok], end = en[ok]),
    loop_id  = paste0(which(ok), "_between")
  )
}

analyze_regions <- function(regions, peaks_control, peaks_treat,
                             bam_control, bam_treat, region_type, loop_category) {
  ov <- findOverlaps(peaks_control, regions)
  if (length(ov) == 0) return(data.frame())

  p              <- peaks_control[queryHits(ov)]
  mcols(p)$region_id <- mcols(regions)$loop_id[subjectHits(ov)]
  pid            <- paste0(as.character(seqnames(p)), ":", start(p), "-", end(p))

  cs  <- bamCount(bam_control, p, paired.end = "midpoint")
  ts  <- bamCount(bam_treat,   p, paired.end = "midpoint")

  l2  <- dplyr::case_when(
    cs == 0 & ts == 0 ~ 0,
    cs == 0            ~ log2(ts + 1),
    ts == 0            ~ -log2(cs + 1),
    TRUE               ~ log2((ts + 1) / (cs + 1))
  )

  data.frame(
    loop_category  = loop_category,
    region_type    = region_type,
    region_id      = mcols(p)$region_id,
    peak_id        = pid,
    log2FC         = l2,
    control_signal = cs,
    treat_signal   = ts,
    peak_score     = mcols(p)$signalValue
  ) |>
    dplyr::group_by(peak_id) |>
    dplyr::slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

protein <- "H3K27ac"
cp      <- peak_list[[paste0(protein, "_control")]]
tp      <- peak_list[[paste0(protein, "_sorbitol")]]
cb      <- bam_files[paste0(protein, "_control")]
tb      <- bam_files[paste0(protein, "_sorbitol")]

results_list <- list()
for (lt in c("Gained", "Lost")) {
  sub <- if (lt == "Gained") gainedLoops else lostLoops
  ar  <- c(anchors(sub, "first"), anchors(sub, "second"))
  mcols(ar)$loop_id <- c(
    paste0(seq_along(anchors(sub, "first")),  "_anchor1"),
    paste0(seq_along(anchors(sub, "second")), "_anchor2")
  )
  br  <- get_between_regions(sub)
  r1  <- analyze_regions(ar, cp, tp, cb, tb, "At Anchors",      lt)
  r2  <- analyze_regions(br, cp, tp, cb, tb, "Between Anchors", lt)
  results_list[[paste0(protein, "_", lt)]] <- bind_rows(r1, r2) |>
    mutate(protein = protein)
}

density_data <- bind_rows(results_list)

## Plot helpers ----
create_ma_plot <- function(res_df, protein_name) {
  n_up <- sum(res_df$isDE == "Increased", na.rm = TRUE)
  n_dn <- sum(res_df$isDE == "Decreased", na.rm = TRUE)

  ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = isDE)) +
    geom_point(alpha = 0.7, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    scale_color_manual(
      values = c("Increased" = "#F8766D", "Decreased" = "#619CFF",
                 "Not significant" = "grey80")
    ) +
    ylim(c(-4, 4)) +
    scale_x_log10(breaks = c(1, 50, 500)) +
    labs(y = paste0(protein_name, " log2\n(sorbitol/control)"),
         x = "mean of normalized counts") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5, face = "bold", size = 9),
          axis.text       = element_text(size = 7),
          axis.title      = element_text(size = 8.5),
          axis.title.y    = element_text(angle = 90, vjust = 0.5),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          plot.margin     = margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio    = 1) +
    annotate("text", label = "Increased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1, y = 3.6,
             color = "#F8766D", fontface = "bold", size = 7 / .pt) +
    annotate("text", label = paste0("n = ", n_up),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1, y = 3.1,
             color = "#F8766D", size = 6 / .pt) +
    annotate("text", label = "Decreased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1, y = -3.4,
             color = "#619CFF", fontface = "bold", size = 7 / .pt) +
    annotate("text", label = paste0("n = ", n_dn),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1, y = -3.9,
             color = "#619CFF", size = 6 / .pt)
}

create_ma_density_plot <- function(res_df) {
  ggplot(res_df, aes(y = log2FoldChange)) +
    geom_density(color = "#9370DB", fill = "#9370DB", alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    ylim(c(-4, 4)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text       = element_blank(),
          axis.title      = element_blank(),
          axis.ticks      = element_blank(),
          axis.line       = element_blank(),
          plot.margin     = margin(0, 0, 0, 0, "pt"),
          panel.spacing   = unit(0, "pt"))
}

create_bar_plot <- function(d) {
  protein_name <- unique(d$Target)[1]
  ggplot(d, aes(x = Category, y = Percentage, fill = Condition)) +
    geom_hline(yintercept = seq(0, 100, 25),
               color = "gray90", linetype = "dashed") +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position  = position_dodge(0.7),
                  width     = 0.25,
                  linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_dodge(0.7),
              vjust    = -1.2,
              size     = 6 / .pt) +
    scale_fill_manual(values = c(control = "#619CFF", sorbitol = "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text       = element_text(size = 7),
          axis.title      = element_text(size = 8.5),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          plot.margin     = margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio    = 1) +
    labs(y = paste0("% of anchors bound\nby ", protein_name),
         x = "Loop Anchor Type") +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.05)),
                       breaks = seq(0, 100, 25)) +
    annotate("text", x = 0.5, y = 95, label = "control",
             color = "#619CFF", fontface = "bold", size = 7 / .pt, hjust = 0) +
    annotate("text", x = 0.5, y = 88, label = "sorbitol",
             color = "#F8766D", fontface = "bold", size = 7 / .pt, hjust = 0)
}

create_density_plot <- function(df) {
  protein_name <- unique(df$protein)[1]
  df <- df |>
    mutate(region_type = factor(region_type,
                                levels = c("At Anchors", "Between Anchors")))

  at_anchors_color      <- "#5DA5DA"
  between_anchors_color <- "#FAA43A"

  dens_at      <- density(df$log2FC[df$region_type == "At Anchors"])
  dens_between <- density(df$log2FC[df$region_type == "Between Anchors"])

  max_at_x      <- dens_at$x[which.max(dens_at$y)]
  max_between_x <- dens_between$x[which.max(dens_between$y)]
  max_at_y      <- max(dens_at$y)
  max_between_y <- max(dens_between$y)

  at_label_x      <- max_at_x + 0.8
  at_label_y      <- max_at_y * 0.8
  between_label_x <- max_between_x + 0.9
  between_label_y <- max_between_y * 0.75

  ggplot(df, aes(x = log2FC, fill = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray75", linewidth = 0.3) +
    geom_density(alpha = 0.4, color = NA) +
    scale_fill_manual(
      values = c("At Anchors" = at_anchors_color,
                 "Between Anchors" = between_anchors_color),
      name = ""
    ) +
    labs(x = paste0(protein_name, " log2 (sorbitol/control)"),
         y = "Density") +
    theme_classic() +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = "none",
          axis.text       = element_text(size = 7),
          axis.title      = element_text(size = 8.5),
          axis.line       = element_line(linewidth = 0.3),
          axis.ticks      = element_line(linewidth = 0.3),
          plot.margin     = margin(5.5, 5.5, 5.5, 5.5, "pt"),
          aspect.ratio    = 1) +
    coord_cartesian(xlim = c(-4.5, 4.5)) +
    scale_x_continuous(limits = c(-4.5, 4.5)) +
    annotate("text", x = at_label_x, y = at_label_y, label = "At\nAnchors",
             color = at_anchors_color, fontface = "bold", size = 7 / .pt,
             vjust = 0.5, hjust = 0, lineheight = 0.8) +
    annotate("text", x = between_label_x, y = between_label_y,
             label = "Between\nAnchors",
             color = between_anchors_color, fontface = "bold", size = 7 / .pt,
             vjust = 0.5, hjust = 0, lineheight = 0.8)
}

# Visualization ----
panel_width  <- 1.8
panel_height <- 1.8
panel_spacing <- 0.2

ma_panel_width  <- panel_width  * 1.07
ma_panel_height <- panel_height * 1.07
ma_dx <- -0.08
ma_dy <- -0.06

density_panel_width  <- 0.13
density_panel_height <- ma_panel_height * 0.77
density_spacing      <- -0.05
density_y_offset     <- ma_panel_height * 0.03

x_start <- 0.25
y_start <- 0.5

x_col1 <- x_start
x_col2 <- x_col1 + panel_width + panel_spacing + 0.05
x_col3 <- x_col2 + panel_width + panel_spacing

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

plotText(label = "A", x = x_col1 - 0.15, y = y_start - 0.2,
         fontsize = 12, fontface = "bold")
plotGG(create_ma_plot(h3k27ac_ma_data, "H3K27ac"),
       x = x_col1 + ma_dx, y = y_start + ma_dy,
       width = ma_panel_width, height = ma_panel_height)
plotGG(create_ma_density_plot(h3k27ac_ma_data),
       x = (x_col1 + ma_dx) + ma_panel_width + density_spacing,
       y = (y_start + ma_dy) + density_y_offset,
       width = density_panel_width, height = density_panel_height)

plotText(label = "B", x = x_col2 - 0.15, y = y_start - 0.2,
         fontsize = 12, fontface = "bold")
plotGG(create_bar_plot(bar_plot_data),
       x = x_col2 + ma_dx, y = y_start + ma_dy,
       width = ma_panel_width, height = ma_panel_height)

plotText(label = "C", x = x_col3 - 0.15, y = y_start - 0.2,
         fontsize = 12, fontface = "bold")
plotGG(create_density_plot(density_data |>
                             filter(protein == "H3K27ac", loop_category == "Gained")),
       x = x_col3 + ma_dx, y = y_start + ma_dy,
       width = ma_panel_width, height = ma_panel_height)

# Save outputs ----
dev.off()

sessionInfo()

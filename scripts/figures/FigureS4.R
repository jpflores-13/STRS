## Figure S4 - H3K27ac and YAP1

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

page_width <- 7.5
page_height <- 5.5

pdf("figures/FigureS4.pdf", width = page_width, height = page_height)
pageCreate(width = page_width, height = page_height, showGuides = FALSE)

gray_color <- "#666666"

# Load Data ---------------------------------------------------------------

loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |>
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

# Categorize loops --------------------------------------------------------

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0] |> 
  as.data.frame() |> 
  as_ginteractions()
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0] |> 
  as.data.frame() |> 
  as_ginteractions()

# Load peaks and BAM files ------------------------------------------------

cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE, pattern = ".narrowPeak") |>
  lapply(read_narrowpeaks)

target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
names(cutntag) <- paste0(rep(target, each = 2), "_", condition)

# Set up BAM files
bam_files <- character()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    fp <- list.files("data/processed/cutntag/output/mergeAlign/", 
                     full.names = TRUE, 
                     pattern = pattern)
    fp <- fp[grepl("\\.bam$", fp)]
    if (length(fp) == 1) {
      bam_files[paste0(t, "_", cond)] <- fp
    }
  }
}

# Create MA plot data for H3K27ac and YAP1 --------------------------------

## Function to merge peaks for a single protein
merge_protein_peaks <- function(protein, peak_list) {
  control_name <- paste0(protein, "_control")
  treatment_name <- paste0(protein, "_sorbitol")
  
  control_idx <- which(names(peak_list) == control_name)
  treatment_idx <- which(names(peak_list) == treatment_name)
  
  combined_peaks <- GenomicRanges::reduce(
    c(peak_list[[control_idx]], peak_list[[treatment_idx]])
  )
  
  return(combined_peaks)
}

## Function to create MA data
create_ma_data_from_peaks <- function(protein, peak_list, bam_files, 
                                      min_mean_count = 5) {
  merged_peaks <- merge_protein_peaks(protein, peak_list)
  
  control_bam <- bam_files[paste0(protein, "_control")]
  treatment_bam <- bam_files[paste0(protein, "_sorbitol")]
  
  control_counts <- bamCount(control_bam, merged_peaks, paired.end = "midpoint")
  treatment_counts <- bamCount(treatment_bam, merged_peaks, paired.end = "midpoint")
  
  mean_counts <- (control_counts + treatment_counts) / 2
  keep_peaks <- mean_counts >= min_mean_count
  
  control_counts <- control_counts[keep_peaks]
  treatment_counts <- treatment_counts[keep_peaks]
  mean_counts <- mean_counts[keep_peaks]
  
  log2fc <- log2((treatment_counts + 1) / (control_counts + 1))
  
  isDE <- case_when(
    log2fc > 1 ~ "Increased",
    log2fc < -1 ~ "Decreased",
    TRUE ~ "Not significant"
  )
  
  data.frame(
    baseMean = mean_counts,
    log2FoldChange = log2fc,
    padj = NA,  # Not calculating p-values for this approach
    isDE = isDE
  ) |>
    arrange(isDE)
}

## Generate MA data for H3K27ac and YAP1
h3k27ac_ma_data <- create_ma_data_from_peaks("H3K27ac", cutntag, bam_files)
yap1_ma_data <- create_ma_data_from_peaks("YAP1", cutntag, bam_files)

# Bar plot data -----------------------------------------------------------

extractAnchors <- function(gi) {
  unique(c(anchors(gi, "first"), anchors(gi, "second")))
}

calculateCI <- function(anchors, peaks, n_bootstrap = 1000, conf_level = 0.95) {
  n_anchors <- length(anchors)
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
    Category = category,
    Target = sub("_.*", "", names(props)),
    Condition = sub(".*_", "", names(props)),
    Percentage = unlist(props) * 100,
    Lower = sapply(cis, \(x) x$lower),
    Upper = sapply(cis, \(x) x$upper)
  )
}

gained_anchors_bar <- extractAnchors(loops[loops$padj < 0.1 & loops$log2FoldChange > 0])
lost_anchors_bar <- extractAnchors(loops[loops$padj < 0.1 & loops$log2FoldChange < 0])

if (!exists("bar_results")) {
  bar_results <- rbind(
    calculateOverlaps(gained_anchors_bar, cutntag, "Gained"),
    calculateOverlaps(lost_anchors_bar, cutntag, "Lost")
  )
}

bar_plot_data <- bar_results |>
  mutate(
    Category = factor(Category, levels = c("Lost", "Gained")),
    Target = factor(Target, levels = c("CTCF", "RAD21", "H3K27ac", "YAP1"))
  )

# Density analysis data ---------------------------------------------------

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    fp <- list.files("data/processed/cutntag/output/peaks/", 
                     full.names = TRUE, 
                     pattern = pattern)
    fp <- fp[grepl("\\.narrowPeak$", fp)]
    if (length(fp) == 1) {
      peak_list[[paste0(t, "_", cond)]] <- read_narrowpeaks(fp)
    }
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
    ranges = IRanges(start = st[ok], end = en[ok]),
    loop_id = paste0(which(ok), "_between")
  )
}

analyze_regions <- function(regions, peaks_control, peaks_treat, 
                            bam_control, bam_treat, region_type, loop_category) {
  ov <- findOverlaps(peaks_control, regions)
  if (length(ov) == 0) return(data.frame())
  
  p <- peaks_control[queryHits(ov)]
  mcols(p)$region_id <- mcols(regions)$loop_id[subjectHits(ov)]
  pid <- paste0(as.character(seqnames(p)), ":", start(p), "-", end(p))
  
  cs <- bamCount(bam_control, p, paired.end = "midpoint")
  ts <- bamCount(bam_treat, p, paired.end = "midpoint")
  
  l2 <- dplyr::case_when(
    cs == 0 & ts == 0 ~ 0,
    cs == 0 ~ log2(ts + 1),
    ts == 0 ~ -log2(cs + 1),
    TRUE ~ log2((ts + 1) / (cs + 1))
  )
  
  data.frame(
    loop_category = loop_category,
    region_type = region_type,
    region_id = mcols(p)$region_id,
    peak_id = pid,
    log2FC = l2,
    control_signal = cs,
    treat_signal = ts,
    peak_score = mcols(p)$signalValue
  ) |>
    dplyr::group_by(peak_id) |>
    dplyr::slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

results_list <- list()
for (protein in c("H3K27ac", "YAP1")) {
  cp <- peak_list[[paste0(protein, "_control")]]
  tp <- peak_list[[paste0(protein, "_sorbitol")]]
  cb <- bam_files[paste0(protein, "_control")]
  tb <- bam_files[paste0(protein, "_sorbitol")]
  
  for (lt in c("Gained", "Lost")) {
    sub <- if (lt == "Gained") gainedLoops else lostLoops
    ar <- c(anchors(sub, "first"), anchors(sub, "second"))
    mcols(ar)$loop_id <- c(
      paste0(seq_along(anchors(sub, "first")), "_anchor1"),
      paste0(seq_along(anchors(sub, "second")), "_anchor2")
    )
    br <- get_between_regions(sub)
    r1 <- analyze_regions(ar, cp, tp, cb, tb, "At Anchors", lt)
    r2 <- analyze_regions(br, cp, tp, cb, tb, "Between Anchors", lt)
    results_list[[paste0(protein, "_", lt)]] <- bind_rows(r1, r2) |>
      mutate(protein = protein)
  }
}

density_data <- bind_rows(results_list)

# Plot helpers ------------------------------------------------------------

create_ma_plot <- function(res_df, protein_name) {
  n_up <- sum(res_df$isDE == "Increased", na.rm = TRUE)
  n_dn <- sum(res_df$isDE == "Decreased", na.rm = TRUE)
  
  ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = isDE)) +
    geom_point(alpha = 0.7, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "grey40", linewidth = 0.3) +
    scale_color_manual(
      values = c(
        "Increased" = "#F8766D",
        "Decreased" = "#619CFF",
        "Not significant" = "grey80"
      )
    ) +
    ylim(c(-4, 4)) +
    scale_x_log10(breaks = c(1, 50, 500)) +
    labs(
      y = "log2FC(sorbitol/control)",
      x = "mean of normalized counts"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.text = element_text(size = 7.5),
      axis.title = element_text(size = 8.5),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      aspect.ratio = 1
    ) +
    annotate("text",
             label = "Increased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1,
             y = 3.6,
             color = "#F8766D",
             fontface = "bold",
             size = 7 / .pt) +
    annotate("text",
             label = paste0("n = ", n_up),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1,
             y = 3.1,
             color = "#F8766D",
             size = 6 / .pt) +
    annotate("text",
             label = "Decreased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1,
             y = -3.4,
             color = "#619CFF",
             fontface = "bold",
             size = 7 / .pt) +
    annotate("text",
             label = paste0("n = ", n_dn),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.1,
             y = -3.9,
             color = "#619CFF",
             size = 6 / .pt)
}

create_ma_density_plot <- function(res_df) {
  ggplot(res_df, aes(y = log2FoldChange)) +
    geom_density(color = "#9370DB", fill = "#9370DB", alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "grey40", linewidth = 0.3) +
    ylim(c(-4, 4)) +
    xlim(c(0, 1)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt"),
      panel.spacing = unit(0, "pt")
    )
}

create_bar_plot <- function(d) {
  ggplot(d, aes(x = Category, y = Percentage, fill = Condition)) +
    geom_hline(yintercept = seq(0, 100, 25), 
               color = "gray90", linetype = "dashed") +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(0.7),
                  width = 0.25,
                  linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_dodge(0.7),
              vjust = -0.5,
              size = 6 / .pt) +
    scale_fill_manual(values = c(control = "#619CFF", sorbitol = "#F8766D")) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      aspect.ratio = 1
    ) +
    labs(y = "% of anchors bound") +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.05)),
      breaks = seq(0, 100, 25)
    )
}

create_density_plot <- function(df) {
  df <- df |>
    mutate(region_type = factor(region_type, 
                                levels = c("At Anchors", "Between Anchors")))
  
  ggplot(df, aes(x = log2FC, fill = region_type, color = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "gray30", linewidth = 0.3) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    scale_fill_manual(
      values = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"),
      name = ""
    ) +
    scale_color_manual(
      values = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"),
      name = ""
    ) +
    labs(
      title = "Gained Loops",
      x = "log2FC(sorbitol/control)",
      y = "Density"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      aspect.ratio = 1
    ) +
    coord_cartesian(xlim = c(-3, 3))
}

# Panel positioning -------------------------------------------------------

panel_width <- 1.8
panel_height <- 1.8
panel_spacing <- 0.2

ma_panel_width <- panel_width * 1.07
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

# Panels A-C --------------------------------------------------------------

## Panel A - H3K27ac MA plot
plotText(label = "A",
         x = x_col1 - 0.15,
         y = y_start_row1 - 0.2,
         fontsize = 12,
         fontface = "bold")

plotGG(create_ma_plot(h3k27ac_ma_data, "H3K27ac"),
       x = x_col1 + ma_dx,
       y = y_start_row1 + ma_dy,
       width = ma_panel_width,
       height = ma_panel_height)

plotGG(create_ma_density_plot(h3k27ac_ma_data),
       x = (x_col1 + ma_dx) + ma_panel_width + density_spacing,
       y = (y_start_row1 + ma_dy) + density_y_offset,
       width = density_panel_width,
       height = density_panel_height)

## YAP1 MA plot
plotGG(create_ma_plot(yap1_ma_data, "YAP1"),
       x = x_col1 + ma_dx,
       y = y_start_row2 + ma_dy,
       width = ma_panel_width,
       height = ma_panel_height)

plotGG(create_ma_density_plot(yap1_ma_data),
       x = (x_col1 + ma_dx) + ma_panel_width + density_spacing,
       y = (y_start_row2 + ma_dy) + density_y_offset,
       width = density_panel_width,
       height = density_panel_height)

## Bar plots
plotText(label = "B",
         x = x_col2 - 0.15,
         y = y_start_row1 - 0.2,
         fontsize = 12,
         fontface = "bold")

h3k27ac_bar_plot <- create_bar_plot(bar_plot_data |> filter(Target == "H3K27ac"))
yap1_bar_plot <- create_bar_plot(bar_plot_data |> filter(Target == "YAP1"))

plotGG(h3k27ac_bar_plot,
       x = x_col2,
       y = y_start_row1,
       width = panel_width,
       height = panel_height)

plotGG(yap1_bar_plot,
       x = x_col2,
       y = y_start_row2,
       width = panel_width,
       height = panel_height)

## Panel C - Density plots
plotText(label = "C",
         x = x_col3 - 0.15,
         y = y_start_row1 - 0.2,
         fontsize = 12,
         fontface = "bold")

plotGG(create_density_plot(density_data |> 
                             filter(protein == "H3K27ac", loop_category == "Gained")),
       x = x_col3,
       y = y_start_row1,
       width = panel_width,
       height = panel_height)

plotGG(create_density_plot(density_data |> 
                             filter(protein == "YAP1", loop_category == "Gained")),
       x = x_col3,
       y = y_start_row2,
       width = panel_width,
       height = panel_height)

## Row labels/segments
row1_y <- y_start_row1 + panel_height + 0.05
row2_y <- y_start_row2 + panel_height + 0.05

plotSegments(x0 = x_col1,
             y0 = row1_y,
             x1 = x_col3 + panel_width,
             y1 = row1_y,
             default.units = "inches",
             linecolor = gray_color,
             lwd = 1)

plotText("H3K27ac",
         x = (x_col1 + (x_col3 + panel_width)) / 2,
         y = row1_y + 0.1,
         fontcolor = gray_color,
         fontsize = 10)

plotSegments(x0 = x_col1,
             y0 = row2_y,
             x1 = x_col3 + panel_width,
             y1 = row2_y,
             default.units = "inches",
             linecolor = gray_color,
             lwd = 1)

plotText("YAP1",
         x = (x_col1 + (x_col3 + panel_width)) / 2,
         y = row2_y + 0.1,
         fontcolor = gray_color,
         fontsize = 10)

dev.off()

cat("\nFigure created successfully!\n")
cat("Saved to: figures/FigureS4.pdf\n")
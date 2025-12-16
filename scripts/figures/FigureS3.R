## Create multipanel figure: CUT&Tag anchor vs between density analysis
## Three-panel figure showing CTCF and RAD21 binding patterns

# Load required libraries -------------------------------------------------
library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(bamsignals)
library(GenomicRanges)
library(plotgardener)

# Helper Functions --------------------------------------------------------

## Extract loop anchors
get_loop_anchors <- function(loops) {
  first_anchors <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  
  all_anchors <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_id <- c(
    paste0(seq_along(first_anchors), "_anchor1"),
    paste0(seq_along(second_anchors), "_anchor2")
  )
  
  return(all_anchors)
}

## Extract regions BETWEEN loop anchors
get_between_regions <- function(loops) {
  first_anchors <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  
  starts <- end(first_anchors) + 1
  ends <- start(second_anchors) - 1
  
  valid_idx <- ends >= starts
  
  if (sum(valid_idx) == 0) {
    return(GRanges())
  }
  
  between_regions <- GRanges(
    seqnames = seqnames(first_anchors)[valid_idx],
    ranges = IRanges(
      start = starts[valid_idx],
      end = ends[valid_idx]
    ),
    loop_id = paste0(which(valid_idx), "_between")
  )
  
  return(between_regions)
}

## Find overlapping peaks and calculate log2FC
analyze_regions <- function(regions, peaks_control, peaks_treat, 
                            bam_control, bam_treat, region_type, loop_category) {
  
  overlaps <- findOverlaps(peaks_control, regions)
  
  if (length(overlaps) == 0) {
    return(data.frame())
  }
  
  overlapping_peaks <- peaks_control[queryHits(overlaps)]
  mcols(overlapping_peaks)$region_id <- mcols(regions)$loop_id[subjectHits(overlaps)]
  
  peak_ids <- paste0(
    as.character(seqnames(overlapping_peaks)), ":",
    start(overlapping_peaks), "-",
    end(overlapping_peaks)
  )
  
  control_signal <- bamCount(bam_control, overlapping_peaks, paired.end = "midpoint")
  treat_signal <- bamCount(bam_treat, overlapping_peaks, paired.end = "midpoint")
  
  log2FC <- case_when(
    control_signal == 0 & treat_signal == 0 ~ 0,
    control_signal == 0 ~ log2(treat_signal + 1),
    treat_signal == 0 ~ -log2(control_signal + 1),
    TRUE ~ log2((treat_signal + 1) / (control_signal + 1))
  )
  
  result_df <- data.frame(
    loop_category = loop_category,
    region_type = region_type,
    region_id = mcols(overlapping_peaks)$region_id,
    peak_id = peak_ids,
    log2FC = log2FC,
    control_signal = control_signal,
    treat_signal = treat_signal,
    peak_score = mcols(overlapping_peaks)$signalValue,
    stringsAsFactors = FALSE
  )
  
  result_df <- result_df |>
    group_by(peak_id) |>
    slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    ungroup()
  
  return(result_df)
}

# Data Loading ------------------------------------------------------------

## Load loops
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

## Categorize loops
gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]
staticLoops <- loops[loops$padj > 0.1]

cat("Loop categories:\n")
cat("  Gained:", length(gainedLoops), "\n")
cat("  Lost:", length(lostLoops), "\n")
cat("  Static:", length(staticLoops), "\n")

## Load CUT&Tag peaks
target <- c("CTCF", "RAD21")
condition <- c("control", "sorbitol")

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(
      "data/processed/cutntag/output/peaks/",
      full.names = TRUE,
      pattern = pattern
    )
    
    file_path <- file_path[grepl("\\.narrowPeak$", file_path)]
    
    if (length(file_path) == 1) {
      peak_name <- paste0(t, "_", cond)
      peak_list[[peak_name]] <- read_narrowpeaks(file_path)
      cat("Loaded peaks:", peak_name, "\n")
    }
  }
}

## Load BAM files
bam_files <- character()
for (t in target) {
  for (cond in condition) {
    pattern <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(
      "data/processed/cutntag/output/mergeAlign/",
      full.names = TRUE,
      pattern = pattern
    )
    
    file_path <- file_path[grepl("\\.bam$", file_path)]
    
    if (length(file_path) == 1) {
      bam_name <- paste0(t, "_", cond)
      bam_files[bam_name] <- file_path
      cat("Loaded BAM:", bam_name, "\n")
    }
  }
}

# Analysis ----------------------------------------------------------------

# Analyze Gained, Static, and Lost loops
loop_categories <- list(
  Gained = gainedLoops,
  Static = staticLoops,
  Lost = lostLoops
)

results_list <- list()

for (protein in target) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing:", protein, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  control_peaks <- peak_list[[paste0(protein, "_control")]]
  sorbitol_peaks <- peak_list[[paste0(protein, "_sorbitol")]]
  control_bam <- bam_files[paste0(protein, "_control")]
  sorbitol_bam <- bam_files[paste0(protein, "_sorbitol")]
  
  for (category in names(loop_categories)) {
    cat("\nProcessing", category, "loops...\n")
    
    loops_subset <- loop_categories[[category]]
    
    anchors <- get_loop_anchors(loops_subset)
    between <- get_between_regions(loops_subset)
    
    cat("  Number of anchors:", length(anchors), "\n")
    cat("  Number of between regions:", length(between), "\n")
    
    anchor_results <- analyze_regions(
      anchors, control_peaks, sorbitol_peaks,
      control_bam, sorbitol_bam, "At Anchors", category
    )
    
    between_results <- analyze_regions(
      between, control_peaks, sorbitol_peaks,
      control_bam, sorbitol_bam, "Between Anchors", category
    )
    
    combined_results <- rbind(anchor_results, between_results) |>
      mutate(protein = protein)
    
    results_list[[paste0(protein, "_", category)]] <- combined_results
    
    cat("  Anchor peaks:", nrow(anchor_results), "\n")
    cat("  Between peaks:", nrow(between_results), "\n")
  }
}

## Combine all results
all_results <- bind_rows(results_list)

# Create plotting function matching Figure 3 style ------------------------

create_density_plot_with_labels <- function(data, protein_name) {
  # Define colors
  at_anchors_color <- "#5DA5DA"
  between_anchors_color <- "#FAA43A"
  
  # Create base plot with labels embedded directly in the plot
  base_plot <- ggplot(data, aes(x = log2FC, fill = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "gray75", linewidth = 0.3) +
    geom_density(alpha = 0.4, color = NA) +
    facet_wrap(~loop_category, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("At Anchors" = at_anchors_color, 
                 "Between Anchors" = between_anchors_color),
      name = ""
    ) +
    labs(
      x = paste0(protein_name, " log2 (sorbitol/control)"),
      y = "Density"
    ) +
    # Add text labels directly to the plot
    annotate("text", x = 0, y = Inf, label = "At\nAnchors",
             color = at_anchors_color, fontface = "bold", 
             size = 7 / .pt, lineheight = 0.8, vjust = 1.5) +
    annotate("text", x = -3, y = Inf, label = "Between\nAnchors",
             color = between_anchors_color, fontface = "bold", 
             size = 7 / .pt, lineheight = 0.8, vjust = 1.5) +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      strip.background = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.5, "lines"),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
    ) +
    coord_cartesian(xlim = c(-4.5, 4.5), clip = "off") +
    scale_x_continuous(limits = c(-4.5, 4.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  return(base_plot)
}

# Create plots for each protein -------------------------------------------

## CTCF plot - Gained, Static, and Lost
ctcf_data <- all_results |> 
  filter(protein == "CTCF") |>
  mutate(
    loop_category = factor(loop_category, levels = c("Gained", "Static", "Lost")),
    region_type = factor(region_type, levels = c("At Anchors", "Between Anchors"))
  )

ctcf_plot <- create_density_plot_with_labels(ctcf_data, "CTCF")

## RAD21 plot - Gained, Static, and Lost
rad21_data <- all_results |> 
  filter(protein == "RAD21") |>
  mutate(
    loop_category = factor(loop_category, levels = c("Gained", "Static", "Lost")),
    region_type = factor(region_type, levels = c("At Anchors", "Between Anchors"))
  )

rad21_plot <- create_density_plot_with_labels(rad21_data, "RAD21")

# Create multipanel figure with plotgardener ------------------------------

## Create output directory if needed
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## Define page dimensions - increased height for 3 facets
page_width <- 8.5
page_height <- 8

## Initialize PDF device
pdf("figures/FigureS3.pdf", width = page_width, height = page_height)

## Create page
pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## Define panel dimensions - increased height for 3 facets
panel_width <- 4
panel_height <- 3.5
panel_spacing <- 0.25

## Panel A: CTCF
x_col1 <- 0.25
x_col2 <- x_col1 + panel_width + panel_spacing
y_start <- 0.5

plotText(
  label = "A",
  fontsize = 12,
  fontface = "bold",
  x = x_col1 - 0.15,
  y = y_start - 0.2
)

plotGG(
  plot = ctcf_plot,
  x = x_col1,
  y = y_start,
  width = panel_width,
  height = panel_height
)

## Panel B: RAD21
plotText(
  label = "B",
  fontsize = 12,
  fontface = "bold",
  x = x_col2 - 0.15,
  y = y_start - 0.2
)

plotGG(
  plot = rad21_plot,
  x = x_col2,
  y = y_start,
  width = panel_width,
  height = panel_height
)

## Close device
dev.off()

cat("\n===========================================\n")
cat("FIGURE CREATED: figures/FigureS3.pdf\n")
cat("===========================================\n")
cat("Panel A: CTCF binding at anchors vs between (Gained, Static & Lost loops)\n")
cat("Panel B: RAD21 binding at anchors vs between (Gained, Static & Lost loops)\n")
cat("===========================================\n\n")

## Print summary statistics
cat("\nCTCF Summary Statistics:\n")
ctcf_summary <- ctcf_data |>
  group_by(loop_category, region_type) |>
  summarize(
    median = median(log2FC, na.rm = TRUE),
    mean = mean(log2FC, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )
print(ctcf_summary)

cat("\nRAD21 Summary Statistics:\n")
rad21_summary <- rad21_data |>
  group_by(loop_category, region_type) |>
  summarize(
    median = median(log2FC, na.rm = TRUE),
    mean = mean(log2FC, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )
print(rad21_summary)
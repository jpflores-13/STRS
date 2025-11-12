## Analyzing Cut&Tag peaks at chromatin loop anchors
## Revised analysis to plot fold changes for specific protein peaks

# Load required libraries ------------------------------------------------
library(InteractionSet)  
library(mariner)         
library(DESeq2)          
library(ggplot2)         
library(dplyr)           
library(plyranges)       
library(tidyr)           
library(forcats)         
library(bamsignals)      
library(GenomicRanges)   
library(GenomeInfoDb)    
library(rtracklayer)     

# Data --------------------------------------------------------------------
## Read in narrowPeak files
cutntag <- list.files("data/processed/cutntag/output/peaks/",
                      full.names = TRUE,
                      pattern = ".narrowPeak")

## Set up sample names
target <- c("CTCF", "H3K27ac", "RAD21", "YAP1")
condition <- c("control", "sorbitol")
peak_names <- paste0(rep(target, each = 2), "_", condition)

## Read in peaks as GRanges objects
peak_list <- list()
for(i in seq_along(cutntag)) {
  peaks <- read_narrowpeaks(cutntag[i])
  peak_list[[peak_names[i]]] <- peaks
}

## Set up paths to BAM files for signal extraction
bam_files <- list.files("data/processed/cutntag/output/mergeAlign/",
                        full.names = TRUE,
                        pattern = ".bam$")
names(bam_files) <- peak_names

## Read in differential loops from HiC analysis
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

## Separate loops into categories based on differential analysis
gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]

lostLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]

staticLoops <- loops[loops$padj > 0.05]

cat("\nLoop counts:\n")
cat("Gained loops:", length(gainedLoops), "\n")
cat("Lost loops:", length(lostLoops), "\n")
cat("Static loops:", length(staticLoops), "\n")

# Helper functions --------------------------------------------------------

## Function to get both anchors from loops as a single GRanges object
get_all_anchors <- function(loops) {
  first_anchors <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")
  
  all_anchors <- c(first_anchors, second_anchors)
  
  mcols(all_anchors)$loop_id <- c(
    paste0(as.character(seqnames(first_anchors)), ":", start(first_anchors), "-", end(first_anchors)),
    paste0(as.character(seqnames(second_anchors)), ":", start(second_anchors), "-", end(second_anchors))
  )
  
  return(all_anchors)
}

## Function to find peaks that overlap with anchors
find_overlapping_peaks <- function(anchors, peaks) {
  overlaps <- findOverlaps(peaks, anchors)
  
  if (length(overlaps) == 0) {
    return(GRanges())
  }
  
  overlapping_peaks <- peaks[queryHits(overlaps)]
  mcols(overlapping_peaks)$anchor_id <- mcols(anchors)$loop_id[subjectHits(overlaps)]
  
  if (is.null(names(overlapping_peaks)) || all(names(overlapping_peaks) == "")) {
    names(overlapping_peaks) <- paste0(
      as.character(seqnames(overlapping_peaks)), ":",
      start(overlapping_peaks), "-",
      end(overlapping_peaks)
    )
  }
  
  return(overlapping_peaks)
}

## Function to calculate signal in peaks from BAM files
calculate_peak_signal <- function(peaks, bam_file) {
  counts <- bamCount(
    bam_file,
    peaks,
    paired.end = "midpoint"
  )
  
  return(counts)
}

## Function to calculate log2FC for peaks between conditions
calculate_peak_lfc <- function(peaks, control_bam, treatment_bam) {
  control_signal <- calculate_peak_signal(peaks, control_bam)
  treatment_signal <- calculate_peak_signal(peaks, treatment_bam)
  
  log2FC <- case_when(
    control_signal == 0 & treatment_signal == 0 ~ 0,
    control_signal == 0 ~ log2(treatment_signal + 1),
    treatment_signal == 0 ~ -log2(control_signal + 1),
    TRUE ~ log2((treatment_signal + 1)/(control_signal + 1))
  )
  
  return(log2FC)
}

# Analysis ----------------------------------------------------------------

## Get all anchors from each loop category
gained_anchors <- get_all_anchors(gainedLoops)
lost_anchors <- get_all_anchors(lostLoops)
static_anchors <- get_all_anchors(staticLoops)

mcols(gained_anchors)$loop_type <- "Gained"
mcols(lost_anchors)$loop_type <- "Lost"
mcols(static_anchors)$loop_type <- "Static"

all_anchors <- c(gained_anchors, lost_anchors, static_anchors)

cat("\nAnchor counts:\n")
cat("Gained loop anchors:", length(gained_anchors), "\n")
cat("Lost loop anchors:", length(lost_anchors), "\n")
cat("Static loop anchors:", length(static_anchors), "\n")

# Create a data frame to store results
results <- data.frame()

# Process each protein target
for (t in target) {
  if (!paste0(t, "_control") %in% names(peak_list)) {
    next
  }
  
  control_peaks <- peak_list[[paste0(t, "_control")]]
  
  if (length(control_peaks) == 0) {
    next
  }
  
  if (!paste0(t, "_control") %in% names(bam_files) || !paste0(t, "_sorbitol") %in% names(bam_files)) {
    next
  }
  
  control_bam <- bam_files[paste0(t, "_control")]
  treatment_bam <- bam_files[paste0(t, "_sorbitol")]
  
  overlapping_peaks <- find_overlapping_peaks(all_anchors, control_peaks)
  
  if (length(overlapping_peaks) == 0) {
    next
  }
  
  peak_lfc <- calculate_peak_lfc(overlapping_peaks, control_bam, treatment_bam)
  
  if (is.null(names(overlapping_peaks)) || all(names(overlapping_peaks) == "")) {
    peak_ids <- paste0(
      as.character(seqnames(overlapping_peaks)), ":",
      start(overlapping_peaks), "-",
      end(overlapping_peaks)
    )
  } else {
    peak_ids <- names(overlapping_peaks)
  }
  
  target_df <- data.frame(
    target = rep(t, length(overlapping_peaks)),
    peak_id = peak_ids,
    anchor_id = mcols(overlapping_peaks)$anchor_id,
    loop_type = mcols(all_anchors)[match(mcols(overlapping_peaks)$anchor_id, mcols(all_anchors)$loop_id), "loop_type"],
    log2FC = peak_lfc,
    peak_score = mcols(overlapping_peaks)$signalValue
  )
  
  target_df <- target_df |>
    group_by(peak_id) |>
    slice_max(order_by = peak_score, n = 1, with_ties = FALSE) |>
    ungroup()
  
  results <- bind_rows(results, target_df)
}

# Create facet labels
facet_labels <- c(
  "CTCF" = "CTCF",
  "H3K27ac" = "H3K27ac",
  "RAD21" = "RAD21",
  "YAP1" = "YAP1"
)

# Print summary of results
cat("\nSummary of overlapping peaks by protein:\n")
summary_by_target <- results |>
  group_by(target) |>
  summarize(
    total_peaks = dplyr::n(),
    gained_peaks = sum(loop_type == "Gained"),
    lost_peaks = sum(loop_type == "Lost"),
    static_peaks = sum(loop_type == "Static"),
    .groups = "drop"
  )
print(summary_by_target)

# Create individual plots for each protein --------------------------------
pdf("plots/cutntag_peak_HiCanchors_individual_plots.pdf", width = 8, height = 8)

protein_plots <- list()

create_protein_plot <- function(data, protein_name) {
  category_counts <- data |>
    group_by(loop_type) |>
    summarize(count = dplyr::n(), .groups = "drop")
  
  p <- ggplot(data, aes(x = fct_rev(loop_type), y = log2FC, fill = loop_type)) +
    geom_hline(yintercept = 0,
               color = "gray",
               linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(
      values = c(
        "Gained" = "#E64B35",
        "Lost" = "#4DBBD5",
        "Static" = "#999999"
      )
    ) +
    scale_x_discrete() +
    labs(
      x = "",
      y = "log2FoldChange(Treated/Control)",
      title = protein_name,
      subtitle = ""
    ) +
    coord_cartesian(ylim = c(-5, 5)) +
    theme(aspect.ratio = 0.8)
  
  return(p)
}

for (t in unique(results$target)) {
  protein_data <- results |> filter(target == t)
  
  if (nrow(protein_data) == 0) {
    next
  }
  
  protein_plot <- create_protein_plot(protein_data, facet_labels[t])
  protein_plots[[t]] <- protein_plot
  print(protein_plot)
}

if (length(protein_plots) > 0) {
  n_plots <- length(protein_plots)
  n_cols <- 2
  n_rows <- 2
  
  combined_plot <- gridExtra::grid.arrange(
    grobs = protein_plots,
    ncol = n_cols,
    nrow = n_rows,
    widths = rep(1, n_cols),
    heights = rep(1, n_rows)
  )
  
  print(combined_plot)
}

dev.off()

summary_stats <- results |>
  group_by(target, loop_type) |>
  summarize(
    mean_lfc = mean(log2FC, na.rm = TRUE),
    median_lfc = median(log2FC, na.rm = TRUE),
    n_peaks = dplyr::n(),
    sd_lfc = sd(log2FC, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nDetailed summary statistics:\n")
print(summary_stats)

cat("\nAnalysis complete! Plots saved to: plots/cutntag_peak_HiCanchors_individual_plots.pdf\n")

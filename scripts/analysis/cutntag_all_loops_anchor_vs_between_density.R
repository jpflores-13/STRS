## Density plot analysis: Anchors vs Between regions across ALL loop categories

library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(bamsignals)
library(GenomicRanges)

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

loop_categories <- list(
  Gained = gainedLoops,
  Lost = lostLoops,
  Static = staticLoops
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

## Save results
saveRDS(all_results, "data/processed/cutntag/all_loops_anchor_vs_between_density_analysis.rds")
write.csv(all_results, "data/processed/cutntag/all_loops_anchor_vs_between_density_analysis.csv", 
          row.names = FALSE)

# Visualization -----------------------------------------------------------

pdf("plots/cutntag_all_loops_anchor_vs_between_density.pdf", width = 10, height = 8)

for (protein in target) {
  protein_data <- all_results |> 
    filter(protein == !!protein) |>
    mutate(
      loop_category = factor(loop_category, levels = c("Static", "Lost", "Gained")),
      region_type = factor(region_type, levels = c("At Anchors", "Between Anchors"))
    )
  
  # Calculate summary statistics for reporting
  summary_stats <- protein_data |>
    group_by(loop_category, region_type) |>
    summarize(
      median = median(log2FC, na.rm = TRUE),
      mean = mean(log2FC, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  cat("\n", protein, "summary statistics:\n")
  print(summary_stats)
  
  # Perform statistical tests
  for (category in c("Static", "Lost", "Gained")) {
    anchor_vals <- protein_data |> 
      filter(loop_category == category, region_type == "At Anchors") |> 
      pull(log2FC)
    between_vals <- protein_data |> 
      filter(loop_category == category, region_type == "Between Anchors") |> 
      pull(log2FC)
    
    if (length(anchor_vals) > 0 && length(between_vals) > 0) {
      test_result <- wilcox.test(anchor_vals, between_vals)
      cat(category, "loops - Wilcoxon test p-value:", test_result$p.value, "\n")
    }
  }
  
  # Create density plots
  p <- ggplot(protein_data, aes(x = log2FC, fill = region_type, color = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    geom_density(alpha = 0.4, linewidth = 1.2) +
    facet_wrap(~loop_category, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"),
      name = ""
    ) +
    scale_color_manual(
      values = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"),
      name = ""
    ) +
    labs(
      title = paste(protein, "Binding Distribution"),
      x = "log2FoldChange(Treated/Untreated)",
      y = "Density"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(1, "lines"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      aspect.ratio = 1
    ) +
    coord_cartesian(xlim = c(-3, 3))
  
  print(p)
}

dev.off()

cat("\n\nAnalysis complete!\n")
cat("Plots saved to: plots/cutntag_all_loops_anchor_vs_between_density.pdf\n")

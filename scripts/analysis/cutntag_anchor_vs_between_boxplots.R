# ##############################################################################
# filename:    cutntag_anchor_vs_between_boxplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Boxplots comparing CTCF and RAD21 log2FC at gained loop anchors
#              vs between-anchor regions; Wilcoxon tests for significance
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(mariner)
library(ggplot2)
library(dplyr)
library(plyranges)
library(bamsignals)
library(GenomicRanges)

# Parameters ----
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
cutntag_peaks_dir <- "data/processed/cutntag/output/peaks/"
bam_files_dir     <- "data/processed/cutntag/output/mergeAlign/"
output_rds        <- "data/processed/cutntag/anchor_vs_between_analysis.rds"
output_csv        <- "data/processed/cutntag/anchor_vs_between_analysis.csv"
output_pdf        <- "plots/cutntag_anchor_vs_between_boxplots.pdf"
target            <- c("CTCF", "RAD21")
condition         <- c("control", "sorbitol")

# Helper functions ----

get_loop_anchors <- function(loops) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")

  all_anchors <- c(first_anchors, second_anchors)
  mcols(all_anchors)$loop_id <- c(
    paste0(seq_along(first_anchors),  "_anchor1"),
    paste0(seq_along(second_anchors), "_anchor2")
  )
  all_anchors
}

get_between_regions <- function(loops) {
  first_anchors  <- anchors(loops, "first")
  second_anchors <- anchors(loops, "second")

  between_regions <- GRanges(
    seqnames = seqnames(first_anchors),
    ranges   = IRanges(start = end(first_anchors) + 1,
                       end   = start(second_anchors) - 1),
    loop_id  = paste0(seq_along(first_anchors), "_between")
  )

  between_regions[width(between_regions) > 0]
}

analyze_regions <- function(regions, peaks_control, peaks_treat,
                            bam_control, bam_treat, region_type) {
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
    control_signal == 0 ~ log2(treat_signal + 1),
    treat_signal   == 0 ~ -log2(control_signal + 1),
    TRUE ~ log2((treat_signal + 1) / (control_signal + 1))
  )

  data.frame(
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

# Data import ----
loops <- readRDS(diff_loops_rds) |> interactions()

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
cat("Number of gained loops:", length(gainedLoops), "\n")

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern   <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = pattern)
    file_path <- file_path[grepl("\\.narrowPeak$", file_path)]
    if (length(file_path) == 1) {
      peak_list[[paste0(t, "_", cond)]] <- read_narrowpeaks(file_path)
      cat("Loaded peaks:", t, cond, "\n")
    } else if (length(file_path) == 0) {
      warning("No peak file found for ", t, "_", cond)
    } else {
      warning("Multiple peak files found for ", t, "_", cond)
    }
  }
}

bam_files <- character()
for (t in target) {
  for (cond in condition) {
    pattern   <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(bam_files_dir, full.names = TRUE, pattern = pattern)
    file_path <- file_path[grepl("\\.bam$", file_path)]
    if (length(file_path) == 1) {
      bam_files[paste0(t, "_", cond)] <- file_path
      cat("Loaded BAM:", t, cond, "\n")
    } else if (length(file_path) == 0) {
      warning("No BAM file found for ", t, "_", cond)
    } else {
      warning("Multiple BAM files found for ", t, "_", cond)
    }
  }
}

# Analysis ----
gained_anchors <- get_loop_anchors(gainedLoops)
gained_between <- get_between_regions(gainedLoops)

cat("\nGained loop regions:\n")
cat("  Anchors:", length(gained_anchors), "\n")
cat("  Between regions:", length(gained_between), "\n")
cat("  Mean between region size:", mean(width(gained_between)), "bp\n")

results_list <- list()

for (protein in target) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing:", protein, "\n")
  cat(rep("=", 60), "\n", sep = "")

  control_peaks  <- peak_list[[paste0(protein, "_control")]]
  sorbitol_peaks <- peak_list[[paste0(protein, "_sorbitol")]]
  control_bam    <- bam_files[paste0(protein, "_control")]
  sorbitol_bam   <- bam_files[paste0(protein, "_sorbitol")]

  anchor_results  <- analyze_regions(gained_anchors, control_peaks, sorbitol_peaks,
                                     control_bam, sorbitol_bam, "Loop Anchor")
  between_results <- analyze_regions(gained_between, control_peaks, sorbitol_peaks,
                                     control_bam, sorbitol_bam, "Between Anchors")

  results_list[[protein]] <- rbind(anchor_results, between_results) |>
    mutate(protein = protein)

  cat("  Anchor peaks (unique):", nrow(anchor_results), "\n")
  cat("  Between peaks (unique):", nrow(between_results), "\n")
}

all_results <- bind_rows(results_list)

# Save outputs ----
saveRDS(all_results, output_rds)
write.csv(all_results, output_csv, row.names = FALSE)

pdf(output_pdf, width = 8, height = 6)

for (protein in target) {
  protein_data <- all_results |> filter(protein == !!protein)

  stats <- protein_data |>
    group_by(region_type) |>
    summarize(median = median(log2FC, na.rm = TRUE),
              mean   = mean(log2FC, na.rm = TRUE),
              n      = n(), .groups = "drop")

  cat("\n", protein, "statistics:\n")
  print(stats)

  anchor_vals  <- protein_data |> filter(region_type == "Loop Anchor")     |> pull(log2FC)
  between_vals <- protein_data |> filter(region_type == "Between Anchors") |> pull(log2FC)

  if (length(anchor_vals) > 0 && length(between_vals) > 0) {
    cat("Wilcoxon test p-value:", wilcox.test(anchor_vals, between_vals)$p.value, "\n")
  }

  p <- ggplot(protein_data, aes(x = region_type, y = log2FC, fill = region_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = c("Loop Anchor" = "#E64B35", "Between Anchors" = "#4DBBD5")) +
    labs(title    = paste(protein, "Binding: Loop Anchors vs Between Regions"),
         subtitle = paste("Gained loops (n =", length(gainedLoops), ")"),
         x = "", y = "log2FoldChange(Treated/Control)") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x   = element_text(size = 11)) +
    coord_cartesian(ylim = c(-3, 3))

  print(p)
}

dev.off()

cat("\n\nAnalysis complete!\n")
cat("Results:", output_csv, "\n")
cat("Plots:", output_pdf, "\n")

sessionInfo()

# ##############################################################################
# filename:    cutntag_gained_loops_anchor_vs_between_density.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Density plot analysis of CTCF and RAD21 log2FC at gained loop
#              anchors vs between-anchor regions; Wilcoxon test for each protein
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
output_rds        <- "data/processed/cutntag/gained_loops_anchor_vs_between_density_analysis.rds"
output_csv        <- "data/processed/cutntag/gained_loops_anchor_vs_between_density_analysis.csv"
output_pdf        <- "plots/cutntag_gained_loops_anchor_vs_between_density.pdf"
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

  starts    <- end(first_anchors) + 1
  ends      <- start(second_anchors) - 1
  valid_idx <- ends >= starts

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
    control_signal == 0 ~ log2(treat_signal + 1),
    treat_signal   == 0 ~ -log2(control_signal + 1),
    TRUE ~ log2((treat_signal + 1) / (control_signal + 1))
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

# Data import ----
loops <- readRDS(diff_loops_rds) |> interactions()

gainedLoops <- loops[loops$padj < 0.1 & loops$log2FoldChange > 0]
lostLoops   <- loops[loops$padj < 0.1 & loops$log2FoldChange < 0]
staticLoops <- loops[loops$padj > 0.1]

cat("Loop categories:\n")
cat("  Gained:", length(gainedLoops), "\n")
cat("  Lost:",   length(lostLoops),   "\n")
cat("  Static:", length(staticLoops), "\n")

peak_list <- list()
for (t in target) {
  for (cond in condition) {
    pattern   <- paste0(t, "_", ifelse(cond == "control", "cont", cond))
    file_path <- list.files(cutntag_peaks_dir, full.names = TRUE, pattern = pattern)
    file_path <- file_path[grepl("\\.narrowPeak$", file_path)]
    if (length(file_path) == 1) {
      peak_list[[paste0(t, "_", cond)]] <- read_narrowpeaks(file_path)
      cat("Loaded peaks:", t, cond, "\n")
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
    }
  }
}

# Analysis ----
results_list <- list()

for (protein in target) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Analyzing:", protein, "\n")
  cat(rep("=", 60), "\n", sep = "")

  control_peaks  <- peak_list[[paste0(protein, "_control")]]
  sorbitol_peaks <- peak_list[[paste0(protein, "_sorbitol")]]
  control_bam    <- bam_files[paste0(protein, "_control")]
  sorbitol_bam   <- bam_files[paste0(protein, "_sorbitol")]

  cat("\nProcessing Gained loops...\n")
  anchors_gr <- get_loop_anchors(gainedLoops)
  between_gr <- get_between_regions(gainedLoops)

  cat("  Anchors:", length(anchors_gr), "\n")
  cat("  Between regions:", length(between_gr), "\n")

  anchor_results  <- analyze_regions(anchors_gr, control_peaks, sorbitol_peaks,
                                     control_bam, sorbitol_bam, "At Anchors", "Gained")
  between_results <- analyze_regions(between_gr, control_peaks, sorbitol_peaks,
                                     control_bam, sorbitol_bam, "Between Anchors", "Gained")

  results_list[[paste0(protein, "_Gained")]] <- rbind(anchor_results, between_results) |>
    mutate(protein = protein)

  cat("  Anchor peaks:", nrow(anchor_results), "\n")
  cat("  Between peaks:", nrow(between_results), "\n")
}

all_results <- bind_rows(results_list)

# Save outputs ----
saveRDS(all_results, output_rds)
write.csv(all_results, output_csv, row.names = FALSE)

pdf(output_pdf, width = 8, height = 5)

for (protein in target) {
  protein_data <- all_results |>
    filter(protein == !!protein) |>
    mutate(region_type = factor(region_type, levels = c("At Anchors", "Between Anchors")))

  summary_stats <- protein_data |>
    group_by(region_type) |>
    summarize(median = median(log2FC, na.rm = TRUE),
              mean   = mean(log2FC, na.rm = TRUE),
              n      = n(), .groups = "drop")

  cat("\n", protein, "summary statistics:\n")
  print(summary_stats)

  anchor_vals  <- protein_data |> filter(region_type == "At Anchors")     |> pull(log2FC)
  between_vals <- protein_data |> filter(region_type == "Between Anchors") |> pull(log2FC)

  if (length(anchor_vals) > 0 && length(between_vals) > 0) {
    cat("Gained loops - Wilcoxon p:", wilcox.test(anchor_vals, between_vals)$p.value, "\n")
  }

  p <- ggplot(protein_data, aes(x = log2FC, fill = region_type, color = region_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    geom_density(alpha = 0.4, linewidth = 1.2) +
    scale_fill_manual(values  = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"), name = "") +
    scale_color_manual(values = c("At Anchors" = "#5DA5DA", "Between Anchors" = "#FAA43A"), name = "") +
    labs(title = paste(protein, "Binding Distribution - Gained Loops"),
         x = "log2FoldChange(Treated/Untreated)", y = "Density") +
    theme_classic() +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "bottom",
          axis.line  = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          aspect.ratio = 1) +
    coord_cartesian(xlim = c(-3, 3))

  print(p)
}

dev.off()

cat("\n\nAnalysis complete!\n")
cat("Results:", output_csv, "\n")
cat("Plots:", output_pdf, "\n")

sessionInfo()

# ##############################################################################
# filename:    ctcf_cooccupancy_pairedBarplots.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Paired barplots showing co-occupancy of retained vs lost CTCF
#              peaks with CTCF, RAD21, H3K27ac, and YAP1 peaks; chi-squared
#              tests for independence at each factor
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters ----
peaks_dir  <- "data/processed/cutntag/output/peaks"
plots_dir  <- "plots"
deseq_dir  <- "data/processed/cutntag/deseq2"
condition  <- "cont"
proteins   <- c("CTCF", "RAD21", "H3K27ac", "YAP1")

# Data import ----
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## Helper to read narrowPeak files
read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name",
                                    "score", "strand", "signalValue",
                                    "pValue", "qValue", "peak"))

  GRanges(seqnames = peaks$seqnames,
          ranges = IRanges(start = peaks$start + 1, end = peaks$end),
          strand = "*",
          score = peaks$score,
          signalValue = peaks$signalValue,
          pValue = peaks$pValue,
          qValue = peaks$qValue,
          peak = peaks$peak)
}

cat("Loading CTCF control peaks...\n")

ctcf_control_file <- list.files(
  path = peaks_dir,
  pattern = paste0(".*CTCF_", condition, ".*\\.narrowPeak$"),
  full.names = TRUE
)

if (length(ctcf_control_file) == 0) stop("No CTCF control narrowPeak file found!")

ctcf_control <- read_narrowpeaks(ctcf_control_file[1])
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("  CTCF control peaks:", length(ctcf_control), "\n\n")

ctcf_deseq_file <- file.path(deseq_dir, "diff_CTCF_counts.rds")
if (!file.exists(ctcf_deseq_file)) stop("CTCF DESeq2 file not found: ", ctcf_deseq_file)

ctcf_deseq <- readRDS(ctcf_deseq_file)
cat("  DESeq2 regions loaded:", length(ctcf_deseq), "\n\n")

# Analysis ----

## Classify CTCF peaks as Retained vs Lost ----
overlaps <- findOverlaps(ctcf_control, ctcf_deseq)
matched_peaks <- ctcf_control[queryHits(overlaps)]

mcols(matched_peaks)$log2FoldChange <-
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(overlaps)]
mcols(matched_peaks)$padj <-
  mcols(ctcf_deseq)$padj[subjectHits(overlaps)]

cat("  Matched peaks:", length(matched_peaks), "\n\n")

mcols(matched_peaks)$peak_status <- case_when(
  !is.na(matched_peaks$padj) &
    matched_peaks$padj < 0.1 &
    matched_peaks$log2FoldChange > 0 ~ "Gained",
  !is.na(matched_peaks$padj) &
    matched_peaks$padj < 0.1 &
    matched_peaks$log2FoldChange < 0 ~ "Lost",
  TRUE ~ "Static"
)

mcols(matched_peaks)$retention_category <- ifelse(
  matched_peaks$peak_status %in% c("Gained", "Static"), "Retained", "Lost"
)

cat("Peak classification:\n")
cat("  Retained (Gained + Static):", sum(matched_peaks$retention_category == "Retained"), "\n")
cat("  Lost:", sum(matched_peaks$retention_category == "Lost"), "\n\n")

retained_peaks <- matched_peaks[matched_peaks$retention_category == "Retained"]
lost_peaks     <- matched_peaks[matched_peaks$retention_category == "Lost"]

## Load other factor peaks ----
cat("Loading other factor peaks (", condition, " condition)...\n", sep = "")

narrowpeak_files <- list.files(
  path = peaks_dir,
  pattern = paste0(".*_", condition, "_.*\\.narrowPeak$"),
  full.names = TRUE
)

peaks_list <- lapply(proteins, function(prot) {
  prot_file <- grep(prot, narrowpeak_files, value = TRUE)

  if (length(prot_file) == 0) stop("No narrowPeak file found for protein: ", prot)
  if (length(prot_file) > 1) {
    warning("Multiple narrowPeak files found for ", prot, ". Using first one.")
    prot_file <- prot_file[1]
  }

  peaks <- read_narrowpeaks(prot_file)
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  cat("  - ", prot, ": ", length(peaks), " peaks\n", sep = "")
  return(peaks)
})
names(peaks_list) <- proteins

## Calculate overlaps ----
count_overlaps <- function(query_peaks, subject_peaks) {
  overlaps <- findOverlaps(query_peaks, subject_peaks, ignore.strand = TRUE)
  length(unique(queryHits(overlaps)))
}

retained_counts <- sapply(proteins, function(prot) {
  n_overlap <- count_overlaps(retained_peaks, peaks_list[[prot]])
  pct <- round(n_overlap / length(retained_peaks) * 100, 1)
  cat("  - ", prot, ": ", n_overlap, " / ", length(retained_peaks),
      " (", pct, "%)\n", sep = "")
  return(n_overlap)
})

lost_counts <- sapply(proteins, function(prot) {
  n_overlap <- count_overlaps(lost_peaks, peaks_list[[prot]])
  pct <- round(n_overlap / length(lost_peaks) * 100, 1)
  cat("  - ", prot, ": ", n_overlap, " / ", length(lost_peaks),
      " (", pct, "%)\n", sep = "")
  return(n_overlap)
})

## Peaks co-bound by all 4 factors ----
find_all_cobound <- function(query_peaks, peaks_list) {
  cobound_indices <- unique(queryHits(
    findOverlaps(query_peaks, peaks_list[[1]], ignore.strand = TRUE)
  ))

  for (i in 2:length(peaks_list)) {
    current_overlaps <- unique(queryHits(
      findOverlaps(query_peaks, peaks_list[[i]], ignore.strand = TRUE)
    ))
    cobound_indices <- intersect(cobound_indices, current_overlaps)
  }

  return(length(cobound_indices))
}

retained_all4 <- find_all_cobound(retained_peaks, peaks_list)
lost_all4     <- find_all_cobound(lost_peaks,     peaks_list)

## Prepare plot data ----
retained_counts_extended <- c(retained_counts, "All 4" = retained_all4)
lost_counts_extended     <- c(lost_counts,     "All 4" = lost_all4)

plot_data <- data.frame(
  Factor   = rep(c(proteins, "All 4"), 2),
  Count    = c(retained_counts_extended, lost_counts_extended),
  Category = rep(c("Retained", "Lost"), each = length(retained_counts_extended)),
  Total    = rep(c(length(retained_peaks), length(lost_peaks)),
                 each = length(retained_counts_extended))
) |>
  mutate(
    Percentage = Count / Total * 100,
    Factor     = factor(Factor, levels = c(proteins, "All 4"))
  )

# Visualization ----

colors <- c("Retained" = "#7FCDBB", "Lost" = "#E34A33")

p <- ggplot(plot_data, aes(x = Factor, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors) +
  labs(title = "Co-occupancy of Retained vs Lost CTCF Peaks",
       subtitle = sprintf("Retained: n=%d | Lost: n=%d",
                          length(retained_peaks), length(lost_peaks)),
       x = "Factor", y = "Number of CTCF peaks overlapping factor",
       fill = "CTCF Peak Status") +
  theme_bw() +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 11, hjust = 0.5),
    axis.title         = element_text(size = 12),
    axis.text.x        = element_text(size = 11, angle = 0),
    legend.position    = "top",
    legend.title       = element_text(size = 11, face = "bold"),
    panel.grid.minor   = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

p_pct <- ggplot(plot_data, aes(x = Factor, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors) +
  labs(title = "Co-occupancy of Retained vs Lost CTCF Peaks (Percentage)",
       subtitle = sprintf("Retained: n=%d | Lost: n=%d",
                          length(retained_peaks), length(lost_peaks)),
       x = "Factor", y = "Percentage of CTCF peaks overlapping factor",
       fill = "CTCF Peak Status") +
  theme_bw() +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position   = "top",
    panel.grid.minor  = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 100))

# Save outputs ----
pdf_file     <- file.path(plots_dir, "ctcf_cobinding_retained_vs_lost.pdf")
pdf_file_pct <- file.path(plots_dir, "ctcf_cobinding_retained_vs_lost_percentage.pdf")

ggsave(pdf_file,     plot = p,     width = 10, height = 6, device = "pdf")
ggsave(pdf_file_pct, plot = p_pct, width = 10, height = 6, device = "pdf")

cat("  Saved:", pdf_file, "\n")
cat("  Saved:", pdf_file_pct, "\n")

sessionInfo()

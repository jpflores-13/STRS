## Paired Barplots: Co-occupancy Analysis of Retained vs Lost CTCF Peaks
## Shows how many retained or lost CTCF peaks overlap with various factors

# Load libraries ----------------------------------------------------------

library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(tidyr)

# Configuration -----------------------------------------------------------

peaks_dir <- "data/processed/cutntag/output/peaks"
plots_dir <- "plots"
deseq_dir <- "data/processed/cutntag/deseq2"

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

proteins <- c("CTCF", "RAD21", "H3K27ac", "YAP1")
condition <- "cont"  # Use control condition peaks

cat("\n", rep("=", 70), "\n", sep = "")
cat("CO-BINDING ANALYSIS: RETAINED VS LOST CTCF PEAKS\n")
cat(rep("=", 70), "\n\n", sep = "")

# Helper function to read narrowPeak files --------------------------------

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

# Load and classify CTCF peaks --------------------------------------------

cat("Loading CTCF control peaks...\n")

## Load CTCF control peak calls
ctcf_control_file <- list.files(
  path = peaks_dir,
  pattern = paste0(".*CTCF_", condition, ".*\\.narrowPeak$"),
  full.names = TRUE
)

if (length(ctcf_control_file) == 0) {
  stop("No CTCF control narrowPeak file found!")
}

ctcf_control <- read_narrowpeaks(ctcf_control_file[1])
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")

cat("  CTCF control peaks:", length(ctcf_control), "\n\n")

## Load CTCF DESeq2 results
cat("Loading CTCF DESeq2 results...\n")

ctcf_deseq_file <- file.path(deseq_dir, "diff_CTCF_counts.rds")

if (!file.exists(ctcf_deseq_file)) {
  stop("CTCF DESeq2 file not found: ", ctcf_deseq_file)
}

ctcf_deseq <- readRDS(ctcf_deseq_file)
cat("  DESeq2 regions loaded:", length(ctcf_deseq), "\n\n")

## Match control peaks to DESeq2 results
cat("Matching control peaks to DESeq2 results...\n")

overlaps <- findOverlaps(ctcf_control, ctcf_deseq)
matched_peaks <- ctcf_control[queryHits(overlaps)]

## Add DESeq2 statistics
mcols(matched_peaks)$log2FoldChange <- 
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(overlaps)]
mcols(matched_peaks)$padj <- 
  mcols(ctcf_deseq)$padj[subjectHits(overlaps)]

cat("  Matched peaks:", length(matched_peaks), "\n\n")

## Classify peaks as Retained vs Lost
## Following the same logic as the binding strength script:
## Retained = Gained (padj<0.1, log2FC>0) OR Static (padj>=0.1)
## Lost = Lost (padj<0.1, log2FC<0)

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
  matched_peaks$peak_status %in% c("Gained", "Static"),
  "Retained",
  "Lost"
)

cat("Peak classification:\n")
cat("  Retained (Gained + Static):", 
    sum(matched_peaks$retention_category == "Retained"), "\n")
cat("  Lost:", sum(matched_peaks$retention_category == "Lost"), "\n")
cat("    - Gained:", sum(matched_peaks$peak_status == "Gained"), "\n")
cat("    - Static:", sum(matched_peaks$peak_status == "Static"), "\n")
cat("    - Lost:", sum(matched_peaks$peak_status == "Lost"), "\n\n")

## Separate into Retained and Lost peak sets
retained_peaks <- matched_peaks[matched_peaks$retention_category == "Retained"]
lost_peaks <- matched_peaks[matched_peaks$retention_category == "Lost"]

# Load other factor peaks -------------------------------------------------

cat("Loading other factor peaks (", condition, " condition)...\n", sep = "")

## Find narrowPeak files
narrowpeak_files <- list.files(
  path = peaks_dir,
  pattern = paste0(".*_", condition, "_.*\\.narrowPeak$"),
  full.names = TRUE
)

## Load peaks for each protein
peaks_list <- lapply(proteins, function(prot) {
  
  prot_file <- grep(prot, narrowpeak_files, value = TRUE)
  
  if (length(prot_file) == 0) {
    stop("No narrowPeak file found for protein: ", prot)
  }
  
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
cat("\n")

# Calculate overlaps for each category ------------------------------------

cat("Calculating overlaps with other factors...\n\n")

## Function to count overlaps
count_overlaps <- function(query_peaks, subject_peaks) {
  overlaps <- findOverlaps(query_peaks, subject_peaks, ignore.strand = TRUE)
  length(unique(queryHits(overlaps)))
}

## Calculate for Retained peaks
cat("Retained CTCF peaks:\n")
retained_counts <- sapply(proteins, function(prot) {
  n_overlap <- count_overlaps(retained_peaks, peaks_list[[prot]])
  pct <- round(n_overlap / length(retained_peaks) * 100, 1)
  cat("  - ", prot, ": ", n_overlap, " / ", length(retained_peaks), 
      " (", pct, "%)\n", sep = "")
  return(n_overlap)
})

cat("\nLost CTCF peaks:\n")
lost_counts <- sapply(proteins, function(prot) {
  n_overlap <- count_overlaps(lost_peaks, peaks_list[[prot]])
  pct <- round(n_overlap / length(lost_peaks) * 100, 1)
  cat("  - ", prot, ": ", n_overlap, " / ", length(lost_peaks), 
      " (", pct, "%)\n", sep = "")
  return(n_overlap)
})

cat("\n")

# Count peaks with all 4 factors ------------------------------------------

cat("Calculating peaks cobound by all 4 factors...\n")

## Function to find peaks overlapping all factors
find_all_cobound <- function(query_peaks, peaks_list) {
  
  ## Start with indices that overlap first factor
  cobound_indices <- unique(queryHits(
    findOverlaps(query_peaks, peaks_list[[1]], ignore.strand = TRUE)
  ))
  
  ## Iteratively intersect with other factors
  for (i in 2:length(peaks_list)) {
    current_overlaps <- unique(queryHits(
      findOverlaps(query_peaks, peaks_list[[i]], ignore.strand = TRUE)
    ))
    cobound_indices <- intersect(cobound_indices, current_overlaps)
  }
  
  return(length(cobound_indices))
}

retained_all4 <- find_all_cobound(retained_peaks, peaks_list)
lost_all4 <- find_all_cobound(lost_peaks, peaks_list)

cat("  Retained peaks with all 4 factors:", retained_all4, 
    "(", round(retained_all4 / length(retained_peaks) * 100, 1), "%)\n")
cat("  Lost peaks with all 4 factors:", lost_all4,
    "(", round(lost_all4 / length(lost_peaks) * 100, 1), "%)\n\n")

# Prepare data for plotting -----------------------------------------------

## Add "All 4" category to the counts
retained_counts_extended <- c(retained_counts, "All 4" = retained_all4)
lost_counts_extended <- c(lost_counts, "All 4" = lost_all4)

## Create data frame
plot_data <- data.frame(
  Factor = rep(c(proteins, "All 4"), 2),
  Count = c(retained_counts_extended, lost_counts_extended),
  Category = rep(c("Retained", "Lost"), each = length(retained_counts_extended)),
  Total = rep(c(length(retained_peaks), length(lost_peaks)), 
              each = length(retained_counts_extended))
) |>
  mutate(
    Percentage = Count / Total * 100,
    Factor = factor(Factor, levels = c(proteins, "All 4"))
  )

# Create paired barplot ---------------------------------------------------

cat("Creating paired barplot...\n")

## Define colors
colors <- c("Retained" = "#7FCDBB", "Lost" = "#E34A33")

## Create plot
p <- ggplot(plot_data, aes(x = Factor, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors) +
  labs(
    title = "Co-occupancy of Retained vs Lost CTCF Peaks",
    subtitle = sprintf("Retained: n=%d | Lost: n=%d", 
                       length(retained_peaks), length(lost_peaks)),
    x = "Factor",
    y = "Number of CTCF peaks overlapping factor",
    fill = "CTCF Peak Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

## Save plot
pdf_file <- file.path(plots_dir, "ctcf_cobinding_retained_vs_lost.pdf")
ggsave(pdf_file, plot = p, width = 10, height = 6, device = "pdf")

cat("  Saved:", pdf_file, "\n\n")

# Alternative: Percentage plot --------------------------------------------

cat("Creating percentage-based paired barplot...\n")

p_pct <- ggplot(plot_data, aes(x = Factor, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = colors) +
  labs(
    title = "Co-occupancy of Retained vs Lost CTCF Peaks (Percentage)",
    subtitle = sprintf("Retained: n=%d | Lost: n=%d", 
                       length(retained_peaks), length(lost_peaks)),
    x = "Factor",
    y = "Percentage of CTCF peaks overlapping factor",
    fill = "CTCF Peak Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 0),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, 100))

## Save percentage plot
pdf_file_pct <- file.path(plots_dir, "ctcf_cobinding_retained_vs_lost_percentage.pdf")
ggsave(pdf_file_pct, plot = p_pct, width = 10, height = 6, device = "pdf")

cat("  Saved:", pdf_file_pct, "\n\n")

# Summary statistics table ------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("SUMMARY TABLE\n")
cat(rep("=", 70), "\n\n", sep = "")

summary_table <- plot_data |>
  select(Factor, Category, Count, Percentage) |>
  pivot_wider(
    names_from = Category,
    values_from = c(Count, Percentage),
    names_glue = "{Category}_{.value}"
  ) |>
  mutate(
    Difference = Retained_Percentage - Lost_Percentage,
    Ratio = Retained_Percentage / Lost_Percentage
  )

print(summary_table)

cat("\n")

# Statistical tests -------------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("STATISTICAL TESTS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Chi-squared tests for independence:\n")
cat("(Testing if co-occupancy rate differs between Retained vs Lost)\n\n")

for (factor in c(proteins, "All 4")) {
  
  ## Get counts
  retained_with <- plot_data$Count[plot_data$Factor == factor & 
                                     plot_data$Category == "Retained"]
  retained_without <- length(retained_peaks) - retained_with
  
  lost_with <- plot_data$Count[plot_data$Factor == factor & 
                                 plot_data$Category == "Lost"]
  lost_without <- length(lost_peaks) - lost_with
  
  ## Create contingency table
  cont_table <- matrix(
    c(retained_with, retained_without, lost_with, lost_without),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Retained", "Lost"),
      c("With factor", "Without factor")
    )
  )
  
  ## Perform chi-squared test
  test_result <- chisq.test(cont_table)
  
  cat(factor, ":\n")
  cat("  Retained: ", retained_with, "/", length(retained_peaks), 
      " (", round(retained_with/length(retained_peaks)*100, 1), "%)\n", sep = "")
  cat("  Lost:     ", lost_with, "/", length(lost_peaks),
      " (", round(lost_with/length(lost_peaks)*100, 1), "%)\n", sep = "")
  cat("  Chi-squared p-value:", sprintf("%.2e", test_result$p.value), "\n")
  cat("  Effect size (OR):", 
      sprintf("%.2f", (retained_with * lost_without) / (retained_without * lost_with)), 
      "\n\n")
}

# Final summary -----------------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Output files:\n")
cat("  - ", pdf_file, "\n", sep = "")
cat("  - ", pdf_file_pct, "\n\n", sep = "")

cat("Summary:\n")
cat("  Total CTCF control peaks analyzed:", length(matched_peaks), "\n")
cat("  - Retained:", length(retained_peaks), "\n")
cat("  - Lost:", length(lost_peaks), "\n\n")

cat("Key findings:\n")
for (factor in c(proteins[-1], "All 4")) {  # Skip CTCF itself
  retained_pct <- plot_data$Percentage[plot_data$Factor == factor & 
                                         plot_data$Category == "Retained"]
  lost_pct <- plot_data$Percentage[plot_data$Factor == factor & 
                                     plot_data$Category == "Lost"]
  diff <- retained_pct - lost_pct
  
  cat("  - ", factor, ": ", 
      ifelse(diff > 0, "Higher", "Lower"), 
      " in Retained (", sprintf("%.1f%%", abs(diff)), " difference)\n", 
      sep = "")
}

cat("\nSession info:\n")
print(sessionInfo())

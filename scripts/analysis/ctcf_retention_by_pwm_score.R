## PWM Score Analysis: CTCF Motif Strength at Retained vs Lost Peaks
## Calculates CTCF PWM scores at 50bp windows around peak summits
## Compares motif strength between retained and lost CTCF peaks

# Load libraries ----------------------------------------------------------

library(GenomicRanges)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
## Note: JASPAR2020 has more stable API than JASPAR2024
## Install with: BiocManager::install("JASPAR2020")
library(ggplot2)
library(dplyr)

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
          peak = peaks$peak)  # This is the summit offset from start
}

# Load control peak calls -------------------------------------------------

cat("Loading CONTROL peak calls...\n")

## CTCF control peaks
ctcf_control <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

# Load DESeq2 results -----------------------------------------------------

cat("\nLoading DESeq2 differential analysis results...\n")

ctcf_deseq <- readRDS("data/processed/cutntag/deseq2/diff_CTCF_counts.rds")

# Classify CTCF peaks as Retained vs Lost --------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Classifying CTCF peaks\n")
cat(rep("=", 60), "\n", sep = "")

## Find overlaps between CTCF control peaks and CTCF DESeq2 results
ctcf_overlaps <- findOverlaps(ctcf_control, ctcf_deseq)

cat("CTCF control peaks matched to DESeq2 results:", length(ctcf_overlaps), "\n")

## Extract matched control peaks
ctcf_matched <- ctcf_control[queryHits(ctcf_overlaps)]

## Add DESeq2 classification
mcols(ctcf_matched)$log2FoldChange_CTCF <- 
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(ctcf_overlaps)]
mcols(ctcf_matched)$padj_CTCF <- 
  mcols(ctcf_deseq)$padj[subjectHits(ctcf_overlaps)]

## Classify peaks
mcols(ctcf_matched)$peak_status <- case_when(
  !is.na(ctcf_matched$padj_CTCF) & 
    ctcf_matched$padj_CTCF < 0.1 & 
    ctcf_matched$log2FoldChange_CTCF > 0 ~ "Gained",
  !is.na(ctcf_matched$padj_CTCF) & 
    ctcf_matched$padj_CTCF < 0.1 & 
    ctcf_matched$log2FoldChange_CTCF < 0 ~ "Lost",
  !is.na(ctcf_matched$padj_CTCF) ~ "Static",
  TRUE ~ "Filtered"
)

## Retention category
mcols(ctcf_matched)$retention_category <- ifelse(
  ctcf_matched$peak_status %in% c("Gained", "Static"),
  "Retained",
  ifelse(ctcf_matched$peak_status == "Lost", "Lost", "Filtered")
)

## Print summary
cat("\nCTCF peak classification:\n")
cat("  Gained:", sum(ctcf_matched$peak_status == "Gained"), "\n")
cat("  Lost:", sum(ctcf_matched$peak_status == "Lost"), "\n")
cat("  Static:", sum(ctcf_matched$peak_status == "Static"), "\n")
cat("  Retained (Gained+Static):", 
    sum(ctcf_matched$retention_category == "Retained"), "\n")
cat("  Filtered (no DESeq2 data):", 
    sum(ctcf_matched$peak_status == "Filtered"), "\n")

# Create 50bp windows around peak summits ---------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Creating 50bp windows around peak summits\n")
cat(rep("=", 60), "\n", sep = "")

## EXPLANATION:
## The 'peak' column from narrowPeak contains the 0-based offset from chromStart
## to the summit (point of highest enrichment).
## 
## Example: If chromStart=1000 (in 0-based coords) and peak=50,
## then summit is at position 1050 (0-based) or 1051 (1-based).
##
## Since our GRanges uses 1-based coords (start already +1 from chromStart),
## the summit position in 1-based coords is: start + peak
##
## To create a 50bp window centered on summit:
## - summit_center = start + peak (this is the 1-based summit position)
## - window_start = summit_center - 24 (25bp upstream, 1-based)
## - window_end = summit_center + 25 (25bp downstream, 1-based)
## This gives us exactly 50bp centered on the summit

## Calculate summit positions (1-based coordinates)
summit_positions <- start(ctcf_matched) + ctcf_matched$peak

## Create 50bp windows centered on summits
## We want 25bp on each side of the summit
window_size <- 50
flank_size <- (window_size / 2) - 1  # 24bp on each side, plus summit = 50bp

ctcf_summits <- GRanges(
  seqnames = seqnames(ctcf_matched),
  ranges = IRanges(
    start = summit_positions - flank_size,
    end = summit_positions + flank_size + 1
  ),
  strand = strand(ctcf_matched)
)

## Transfer metadata
mcols(ctcf_summits) <- mcols(ctcf_matched)

## Filter to only Retained vs Lost peaks
ctcf_summits_filtered <- ctcf_summits[
  ctcf_summits$retention_category %in% c("Retained", "Lost")
]

cat("Total summit windows created:", length(ctcf_summits), "\n")
cat("Summit windows for analysis (Retained + Lost):", 
    length(ctcf_summits_filtered), "\n")
cat("  Retained:", sum(ctcf_summits_filtered$retention_category == "Retained"), "\n")
cat("  Lost:", sum(ctcf_summits_filtered$retention_category == "Lost"), "\n")

# Get DNA sequences for summit windows ------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Extracting DNA sequences from genome\n")
cat(rep("=", 60), "\n", sep = "")

## EXPLANATION:
## We use BSgenome to extract the actual DNA sequences at our 50bp windows.
## BSgenome.Hsapiens.UCSC.hg38 provides access to the human reference genome.
## getSeq() extracts sequences at the coordinates specified by our GRanges object.

## Extract sequences
summit_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ctcf_summits_filtered)

cat("Sequences extracted:", length(summit_sequences), "\n")
cat("Sequence width (should all be 50bp):", 
    unique(width(summit_sequences)), "\n")

# Get CTCF PWM from JASPAR ------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Loading CTCF Position Weight Matrix (PWM)\n")
cat(rep("=", 60), "\n", sep = "")

## EXPLANATION:
## JASPAR is a database of transcription factor binding profiles.
## Each profile is represented as a Position Weight Matrix (PWM) that
## describes the probability of each nucleotide at each position in the
## binding site. We retrieve the CTCF PWM to score our sequences.
##
## The PWM scoring process:
## 1. For each position in the sequence, calculate a score based on how well
##    it matches the expected nucleotide frequencies in the PWM
## 2. Sum these scores across the entire motif length
## 3. Convert to a relative score (0-1 or percentage)

## Get CTCF PWM from JASPAR
## MA0139.1 is the canonical CTCF matrix
## JASPAR2024 seems to have API issues, trying alternative approaches

## Approach 1: Try JASPAR2020 which has a more stable API
tryCatch({
  if (requireNamespace("JASPAR2020", quietly = TRUE)) {
    library(JASPAR2020)
    opts <- list()
    opts[["species"]] <- 9606  # Human taxon ID
    opts[["name"]] <- "CTCF"
    pfm_list <- getMatrixSet(JASPAR2020, opts)
    ctcf_pfm <- pfm_list[[1]]
    cat("\nUsing JASPAR2020 database\n")
  } else {
    stop("JASPAR2020 not available")
  }
}, error = function(e) {
  ## Approach 2: Try direct object access
  tryCatch({
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]] <- "CTCF"
    pfm_list <- getMatrixSet(JASPAR2024, opts)
    ctcf_pfm <- pfm_list[[1]]
    cat("\nUsing JASPAR2024 database\n")
  }, error = function(e2) {
    stop("Could not access JASPAR database. Please install JASPAR2020:\n  BiocManager::install('JASPAR2020')")
  })
})

cat("Using CTCF matrix:", ID(ctcf_pfm), "-", name(ctcf_pfm), "\n")
cat("Matrix width:", ncol(Matrix(ctcf_pfm)), "bp\n")

## Convert PFM (Position Frequency Matrix) to PWM (Position Weight Matrix)
## PWM is log-odds scored version of PFM
ctcf_pwm <- toPWM(ctcf_pfm, type = "log2probratio", pseudocounts = 0.8)

cat("PWM created successfully\n")

# Score sequences with CTCF PWM -------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("Scoring sequences with CTCF PWM\n")
cat(rep("=", 60), "\n", sep = "")

## EXPLANATION:
## searchSeq() scans each sequence for matches to the CTCF PWM.
## It returns:
## - Best match position in the sequence
## - Match score (how well the sequence matches the PWM)
## - Match strand (+ or -)
## - Relative score (score normalized to 0-1 range)
##
## The relative score represents how well the sequence matches the CTCF
## motif compared to the best possible match (1.0) and worst possible
## match (0.0). Higher scores indicate stronger CTCF binding motifs.

## Score sequences
## min.score="80%" means we report matches with relative score >= 80%
## We use a permissive threshold to ensure we capture all peaks
ctcf_matches <- searchSeq(
  ctcf_pwm, 
  summit_sequences, 
  seqname = as.character(seqnames(ctcf_summits_filtered)),
  min.score = "60%",  # Permissive threshold to capture all peaks
  strand = "*"         # Search both strands
)

cat("Total motif matches found:", length(ctcf_matches), "\n")

## Convert to data frame for easier manipulation
matches_df <- writeGFF3(ctcf_matches) |>
  as.data.frame() |>
  dplyr::select(seqname, start, end, score, strand) |>
  mutate(peak_index = as.integer(seqname))  # seqname contains the peak index

cat("Matches data frame created\n")

## EXPLANATION of score extraction:
## For each peak, we may have multiple PWM matches (both strands, different positions).
## We take the MAXIMUM score for each peak, as this represents the best CTCF motif
## match within that 50bp window. This is the most biologically relevant score
## since CTCF will bind to the strongest available motif.

## Extract best score per peak (maximum across all matches)
best_scores <- matches_df |>
  group_by(peak_index) |>
  summarise(
    max_pwm_score = max(score),
    .groups = "drop"
  )

cat("Best scores calculated for", nrow(best_scores), "peaks\n")

## Add PWM scores to our GRanges object
mcols(ctcf_summits_filtered)$pwm_score <- NA_real_

## Match scores to peaks
mcols(ctcf_summits_filtered)$pwm_score[best_scores$peak_index] <- 
  best_scores$max_pwm_score

cat("PWM scores assigned to peaks\n")
cat("Peaks with PWM scores:", 
    sum(!is.na(ctcf_summits_filtered$pwm_score)), "\n")
cat("Peaks without PWM scores (no match >= 60%):", 
    sum(is.na(ctcf_summits_filtered$pwm_score)), "\n")

# Handle peaks without scores ---------------------------------------------

## EXPLANATION:
## Some peaks may not have any PWM matches above our threshold (60%).
## These could be:
## 1. False positive peaks (not real CTCF sites)
## 2. Peaks with very weak/degraded CTCF motifs
## 3. CTCF binding through indirect mechanisms
##
## For visualization, we have options:
## A. Exclude them (shows only confident CTCF motifs)
## B. Assign them a score of 0 (conservative, assumes no motif)
## C. Score them with a lower threshold
##
## Let's try option C first - rescore with no threshold to get all peaks

cat("\nRescoring peaks without matches using no threshold...\n")

## Get indices of peaks without scores
peaks_without_scores <- which(is.na(ctcf_summits_filtered$pwm_score))
cat("Rescoring", length(peaks_without_scores), "peaks\n")

if (length(peaks_without_scores) > 0) {
  ## Extract sequences for these peaks
  missing_sequences <- summit_sequences[peaks_without_scores]
  
  ## Search with no minimum score threshold
  missing_matches <- searchSeq(
    ctcf_pwm, 
    missing_sequences, 
    min.score = "0%",  # No threshold - get best match regardless
    strand = "*"
  )
  
  ## Extract best scores
  missing_matches_df <- writeGFF3(missing_matches) |>
    as.data.frame() |>
    dplyr::select(seqname, score) |>
    mutate(original_index = peaks_without_scores[as.integer(seqname)]) |>
    group_by(original_index) |>
    summarise(max_pwm_score = max(score), .groups = "drop")
  
  ## Assign these scores
  mcols(ctcf_summits_filtered)$pwm_score[missing_matches_df$original_index] <- 
    missing_matches_df$max_pwm_score
  
  cat("Rescoring complete\n")
}

cat("\nFinal PWM score statistics:\n")
cat("Peaks with PWM scores:", 
    sum(!is.na(ctcf_summits_filtered$pwm_score)), "\n")
cat("Score range:", 
    sprintf("[%.2f, %.2f]", 
            min(ctcf_summits_filtered$pwm_score, na.rm = TRUE),
            max(ctcf_summits_filtered$pwm_score, na.rm = TRUE)), "\n")

# Prepare data for plotting -----------------------------------------------

## Create data frame
plot_df <- as.data.frame(mcols(ctcf_summits_filtered)) |>
  dplyr::select(retention_category, pwm_score) |>
  filter(!is.na(pwm_score))  # Remove any remaining NAs

cat("\nFinal dataset for plotting:\n")
cat("  Retained CTCF peaks:", 
    sum(plot_df$retention_category == "Retained"), "\n")
cat("  Lost CTCF peaks:", 
    sum(plot_df$retention_category == "Lost"), "\n")

# Statistics --------------------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("STATISTICS\n")
cat(rep("=", 60), "\n", sep = "")

## Helper function to format p-values with significance stars
format_pvalue <- function(p) {
  stars <- case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  if (p < 0.001) {
    label <- sprintf("p < 0.001 %s", stars)
  } else if (p < 0.01) {
    label <- sprintf("p = %.3f %s", p, stars)
  } else {
    label <- sprintf("p = %.3f %s", p, stars)
  }
  
  return(list(label = label, stars = stars))
}

## Statistical test
wilcox_result <- wilcox.test(pwm_score ~ retention_category, data = plot_df)
pval_info <- format_pvalue(wilcox_result$p.value)

cat("\nCTCF PWM scores at peak summits:\n")
cat(sprintf("  Retained CTCF: n=%d, median PWM score=%.2f\n",
            sum(plot_df$retention_category == "Retained"),
            median(plot_df$pwm_score[plot_df$retention_category == "Retained"])))
cat(sprintf("  Lost CTCF: n=%d, median PWM score=%.2f\n",
            sum(plot_df$retention_category == "Lost"),
            median(plot_df$pwm_score[plot_df$retention_category == "Lost"])))
cat(sprintf("  Median difference: %.2f\n",
            median(plot_df$pwm_score[plot_df$retention_category == "Retained"]) -
              median(plot_df$pwm_score[plot_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e %s\n", 
            wilcox_result$p.value, pval_info$stars))

## Additional statistics
cat("\nMean PWM scores:\n")
cat(sprintf("  Retained: %.2f\n",
            mean(plot_df$pwm_score[plot_df$retention_category == "Retained"])))
cat(sprintf("  Lost: %.2f\n",
            mean(plot_df$pwm_score[plot_df$retention_category == "Lost"])))

## Effect size (Cohen's d)
retained_scores <- plot_df$pwm_score[plot_df$retention_category == "Retained"]
lost_scores <- plot_df$pwm_score[plot_df$retention_category == "Lost"]
cohens_d <- (mean(retained_scores) - mean(lost_scores)) / 
  sqrt((var(retained_scores) + var(lost_scores)) / 2)
cat(sprintf("\nEffect size (Cohen's d): %.3f\n", cohens_d))

# Create density plot -----------------------------------------------------

p_density <- ggplot(plot_df, 
                    aes(x = pwm_score, 
                        fill = retention_category,
                        color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  scale_color_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "CTCF Motif Strength at Retained vs Lost Peaks",
    subtitle = sprintf("%s | 50bp windows at peak summits", 
                       pval_info$label),
    x = "CTCF PWM Score",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "top"
  )

# Create violin plot ------------------------------------------------------

## Calculate y-position for significance bracket
y_max <- max(plot_df$pwm_score)
y_range <- diff(range(plot_df$pwm_score))
bracket_y <- y_max + 0.05 * y_range

p_violin <- ggplot(plot_df, 
                   aes(x = retention_category, 
                       y = pwm_score,
                       fill = retention_category)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
  ## Add significance bracket
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
           linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.01 * y_range,
           linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.01 * y_range,
           linewidth = 0.5) +
  annotate("text", x = 1.5, y = bracket_y + 0.02 * y_range, 
           label = pval_info$stars, size = 6) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "CTCF Motif Strength at Retained vs Lost Peaks",
    subtitle = sprintf("%s | n(Retained)=%d, n(Lost)=%d",
                       pval_info$label,
                       sum(plot_df$retention_category == "Retained"),
                       sum(plot_df$retention_category == "Lost")),
    x = "CTCF Peak Status",
    y = "CTCF PWM Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    legend.position = "none"
  )

# Create boxplot ----------------------------------------------------------

p_boxplot <- ggplot(plot_df, 
                    aes(x = retention_category, 
                        y = pwm_score,
                        fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  ## Add significance bracket (reusing same y positions)
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
           linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.01 * y_range,
           linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.01 * y_range,
           linewidth = 0.5) +
  annotate("text", x = 1.5, y = bracket_y + 0.02 * y_range, 
           label = pval_info$stars, size = 6) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    title = "CTCF Motif Strength at Retained vs Lost Peaks",
    # subtitle = sprintf("%s | n(Retained)=%d, n(Lost)=%d",
    #                    pval_info$label,
    #                    sum(plot_df$retention_category == "Retained"),
    #                    sum(plot_df$retention_category == "Lost")),
    x = "CTCF Peak Status",
    y = "CTCF PWM Score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    legend.position = "none"
  )

# Save plots --------------------------------------------------------------

dir.create("plots", showWarnings = FALSE)

pdf("plots/ctcf_retention_by_pwm_score_density.pdf", width = 8, height = 5)
print(p_density)
dev.off()

pdf("plots/ctcf_retention_by_pwm_score_violin.pdf", width = 6, height = 6)
print(p_violin)
dev.off()

pdf("plots/ctcf_retention_by_pwm_score_boxplot.pdf", width = 6, height = 6)
print(p_boxplot)
dev.off()

# Save data ---------------------------------------------------------------

## Save GRanges object with PWM scores
saveRDS(ctcf_summits_filtered, 
        "data/processed/cutntag/ctcf_summits_with_pwm_scores.rds")

## Save plot data
write.table(plot_df,
            "data/processed/cutntag/ctcf_pwm_scores_for_plotting.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("DONE!\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nPlots saved to:\n")
cat("  plots/ctcf_retention_by_pwm_score_density.pdf\n")
cat("  plots/ctcf_retention_by_pwm_score_violin.pdf\n")
cat("  plots/ctcf_retention_by_pwm_score_boxplot.pdf\n")
cat("\nData saved to:\n")
cat("  data/processed/cutntag/ctcf_summits_with_pwm_scores.rds\n")
cat("  data/processed/cutntag/ctcf_pwm_scores_for_plotting.txt\n")
cat("\nThis analysis shows:\n")
cat("  CTCF motif strength (PWM score) at 50bp windows around peak summits\n")
cat("  Comparison between retained and lost CTCF peaks after sorbitol treatment\n")
cat("  Higher PWM scores indicate stronger CTCF binding motifs\n")

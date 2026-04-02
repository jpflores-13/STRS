# ##############################################################################
# filename:    ctcf_retention_by_pwm_score.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Compares CTCF PWM scores at 50bp windows around peak summits
#              between retained vs lost CTCF peaks; density, violin, and
#              boxplot outputs; data exported as RDS and TSV
# ##############################################################################

# Libraries ----
library(GenomicRanges)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(ggplot2)
library(dplyr)

# Parameters ----
ctcf_narrowpeak <- "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
ctcf_deseq_rds  <- "data/processed/cutntag/deseq2/diff_CTCF_counts.rds"
plots_dir       <- "plots"
pwm_scores_rds  <- "data/processed/cutntag/ctcf_summits_with_pwm_scores.rds"
pwm_txt         <- "data/processed/cutntag/ctcf_pwm_scores_for_plotting.txt"

# Helper functions ----

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

# Data import ----
cat("Loading CONTROL peak calls...\n")

ctcf_control <- read_narrowpeaks(ctcf_narrowpeak)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

cat("\nLoading DESeq2 differential analysis results...\n")
ctcf_deseq <- readRDS(ctcf_deseq_rds)

# Analysis ----

## Classify CTCF peaks as Retained vs Lost ----
ctcf_overlaps <- findOverlaps(ctcf_control, ctcf_deseq)
cat("CTCF control peaks matched to DESeq2 results:", length(ctcf_overlaps), "\n")

ctcf_matched <- ctcf_control[queryHits(ctcf_overlaps)]

mcols(ctcf_matched)$log2FoldChange_CTCF <-
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(ctcf_overlaps)]
mcols(ctcf_matched)$padj_CTCF <-
  mcols(ctcf_deseq)$padj[subjectHits(ctcf_overlaps)]

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

mcols(ctcf_matched)$retention_category <- ifelse(
  ctcf_matched$peak_status %in% c("Gained", "Static"),
  "Retained",
  ifelse(ctcf_matched$peak_status == "Lost", "Lost", "Filtered")
)

cat("\nCTCF peak classification:\n")
cat("  Gained:", sum(ctcf_matched$peak_status == "Gained"), "\n")
cat("  Lost:", sum(ctcf_matched$peak_status == "Lost"), "\n")
cat("  Static:", sum(ctcf_matched$peak_status == "Static"), "\n")
cat("  Retained (Gained+Static):",
    sum(ctcf_matched$retention_category == "Retained"), "\n")
cat("  Filtered (no DESeq2 data):",
    sum(ctcf_matched$peak_status == "Filtered"), "\n")

## Create 50bp windows around peak summits ----
# The 'peak' column contains 0-based offset from chromStart to the summit.
# Since our GRanges uses 1-based coords (start already +1), summit = start + peak.
# A 50bp window: start = summit - 24, end = summit + 25.
summit_positions <- start(ctcf_matched) + ctcf_matched$peak
window_size <- 50
flank_size  <- (window_size / 2) - 1

ctcf_summits <- GRanges(
  seqnames = seqnames(ctcf_matched),
  ranges   = IRanges(
    start = summit_positions - flank_size,
    end   = summit_positions + flank_size + 1
  ),
  strand = strand(ctcf_matched)
)
mcols(ctcf_summits) <- mcols(ctcf_matched)

ctcf_summits_filtered <- ctcf_summits[
  ctcf_summits$retention_category %in% c("Retained", "Lost")
]

cat("Total summit windows created:", length(ctcf_summits), "\n")
cat("Summit windows for analysis (Retained + Lost):",
    length(ctcf_summits_filtered), "\n")

## Extract genome sequences ----
summit_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ctcf_summits_filtered)
cat("Sequences extracted:", length(summit_sequences), "\n")

## Load CTCF PWM from JASPAR ----
tryCatch({
  if (requireNamespace("JASPAR2020", quietly = TRUE)) {
    library(JASPAR2020)
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]]    <- "CTCF"
    pfm_list  <- getMatrixSet(JASPAR2020, opts)
    ctcf_pfm  <- pfm_list[[1]]
    cat("\nUsing JASPAR2020 database\n")
  } else {
    stop("JASPAR2020 not available")
  }
}, error = function(e) {
  tryCatch({
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]]    <- "CTCF"
    pfm_list <- getMatrixSet(JASPAR2024, opts)
    ctcf_pfm <<- pfm_list[[1]]
    cat("\nUsing JASPAR2024 database\n")
  }, error = function(e2) {
    stop("Could not access JASPAR database. Install JASPAR2020:\n  BiocManager::install('JASPAR2020')")
  })
})

cat("Using CTCF matrix:", ID(ctcf_pfm), "-", name(ctcf_pfm), "\n")
ctcf_pwm <- toPWM(ctcf_pfm, type = "log2probratio", pseudocounts = 0.8)

## Score sequences with CTCF PWM ----
ctcf_matches <- searchSeq(
  ctcf_pwm,
  summit_sequences,
  seqname   = as.character(seqnames(ctcf_summits_filtered)),
  min.score = "60%",
  strand    = "*"
)

cat("Total motif matches found:", length(ctcf_matches), "\n")

matches_df <- writeGFF3(ctcf_matches) |>
  as.data.frame() |>
  dplyr::select(seqname, start, end, score, strand) |>
  mutate(peak_index = as.integer(seqname))

best_scores <- matches_df |>
  group_by(peak_index) |>
  summarise(max_pwm_score = max(score), .groups = "drop")

mcols(ctcf_summits_filtered)$pwm_score <- NA_real_
mcols(ctcf_summits_filtered)$pwm_score[best_scores$peak_index] <-
  best_scores$max_pwm_score

## Rescore peaks without matches using no threshold ----
peaks_without_scores <- which(is.na(ctcf_summits_filtered$pwm_score))
cat("Rescoring", length(peaks_without_scores), "peaks with no threshold...\n")

if (length(peaks_without_scores) > 0) {
  missing_sequences <- summit_sequences[peaks_without_scores]

  missing_matches <- searchSeq(ctcf_pwm, missing_sequences,
                               min.score = "0%", strand = "*")

  missing_matches_df <- writeGFF3(missing_matches) |>
    as.data.frame() |>
    dplyr::select(seqname, score) |>
    mutate(original_index = peaks_without_scores[as.integer(seqname)]) |>
    group_by(original_index) |>
    summarise(max_pwm_score = max(score), .groups = "drop")

  mcols(ctcf_summits_filtered)$pwm_score[missing_matches_df$original_index] <-
    missing_matches_df$max_pwm_score
}

cat("Peaks with PWM scores:",
    sum(!is.na(ctcf_summits_filtered$pwm_score)), "\n")

## Prepare plot data ----
plot_df <- as.data.frame(mcols(ctcf_summits_filtered)) |>
  dplyr::select(retention_category, pwm_score) |>
  filter(!is.na(pwm_score))

cat("\nFinal dataset: Retained =",
    sum(plot_df$retention_category == "Retained"),
    ", Lost =", sum(plot_df$retention_category == "Lost"), "\n")

## Statistics ----
format_pvalue <- function(p) {
  stars <- case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
  label <- if (p < 0.001) sprintf("p < 0.001 %s", stars) else sprintf("p = %.3f %s", p, stars)
  list(label = label, stars = stars)
}

wilcox_result <- wilcox.test(pwm_score ~ retention_category, data = plot_df)
pval_info     <- format_pvalue(wilcox_result$p.value)

cat("\nCTCF PWM scores at peak summits:\n")
cat(sprintf("  Retained: n=%d, median=%.2f\n",
            sum(plot_df$retention_category == "Retained"),
            median(plot_df$pwm_score[plot_df$retention_category == "Retained"])))
cat(sprintf("  Lost: n=%d, median=%.2f\n",
            sum(plot_df$retention_category == "Lost"),
            median(plot_df$pwm_score[plot_df$retention_category == "Lost"])))
cat(sprintf("  Wilcoxon p-value: %.2e %s\n", wilcox_result$p.value, pval_info$stars))

retained_scores <- plot_df$pwm_score[plot_df$retention_category == "Retained"]
lost_scores     <- plot_df$pwm_score[plot_df$retention_category == "Lost"]
cohens_d <- (mean(retained_scores) - mean(lost_scores)) /
  sqrt((var(retained_scores) + var(lost_scores)) / 2)
cat(sprintf("  Effect size (Cohen's d): %.3f\n", cohens_d))

# Visualization ----

## Density plot ----
p_density <- ggplot(plot_df,
                    aes(x = pwm_score,
                        fill  = retention_category,
                        color = retention_category)) +
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values  = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  scale_color_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                     name = "CTCF Peak Status") +
  labs(title    = "CTCF Motif Strength at Retained vs Lost Peaks",
       subtitle = sprintf("%s | 50bp windows at peak summits", pval_info$label),
       x = "CTCF PWM Score", y = "Density") +
  theme_bw() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "top")

## Violin plot ----
y_max     <- max(plot_df$pwm_score)
y_range   <- diff(range(plot_df$pwm_score))
bracket_y <- y_max + 0.05 * y_range

p_violin <- ggplot(plot_df,
                   aes(x = retention_category,
                       y = pwm_score,
                       fill = retention_category)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.01 * y_range, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.01 * y_range, linewidth = 0.5) +
  annotate("text", x = 1.5, y = bracket_y + 0.02 * y_range, label = pval_info$stars, size = 6) +
  scale_fill_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  labs(title    = "CTCF Motif Strength at Retained vs Lost Peaks",
       subtitle = sprintf("%s | n(Retained)=%d, n(Lost)=%d", pval_info$label,
                          sum(plot_df$retention_category == "Retained"),
                          sum(plot_df$retention_category == "Lost")),
       x = "CTCF Peak Status", y = "CTCF PWM Score") +
  theme_bw() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        legend.position = "none")

## Boxplot ----
p_boxplot <- ggplot(plot_df,
                    aes(x = retention_category,
                        y = pwm_score,
                        fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.01 * y_range, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.01 * y_range, linewidth = 0.5) +
  annotate("text", x = 1.5, y = bracket_y + 0.02 * y_range, label = pval_info$stars, size = 6) +
  scale_fill_manual(values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
                    name = "CTCF Peak Status") +
  labs(title = "CTCF Motif Strength at Retained vs Lost Peaks",
       x = "CTCF Peak Status", y = "CTCF PWM Score") +
  theme_bw() +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 9),
        legend.position = "none")

# Save outputs ----
dir.create(plots_dir, showWarnings = FALSE)

pdf(file.path(plots_dir, "ctcf_retention_by_pwm_score_density.pdf"), width = 8, height = 5)
print(p_density)
dev.off()

pdf(file.path(plots_dir, "ctcf_retention_by_pwm_score_violin.pdf"), width = 6, height = 6)
print(p_violin)
dev.off()

pdf(file.path(plots_dir, "ctcf_retention_by_pwm_score_boxplot.pdf"), width = 6, height = 6)
print(p_boxplot)
dev.off()

saveRDS(ctcf_summits_filtered, pwm_scores_rds)
write.table(plot_df, pwm_txt, quote = FALSE, sep = "\t", row.names = FALSE)

sessionInfo()

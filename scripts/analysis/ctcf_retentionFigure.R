## Figure 4

# Load required libraries -------------------------------------------------

library(plotgardener)
library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)
library(GenomeInfoDb)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(ggsignif)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Helper functions --------------------------------------------------------

## Read narrowPeak files
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

## Get library size from BAM
bam_libsize <- function(bam) {
  sum(idxstatsBam(bam)$mapped, na.rm = TRUE)
}

## Analyze peaks by CPM bins
analyze_by_cpm_bins <- function(peaks_gr,
                                bam_control,
                                bam_treat,
                                protein_name = "CTCF",
                                n_bins = 20,
                                pseudocount_cpm = 1e-6,
                                min_peaks_per_bin = 20,
                                paired_end_mode = "midpoint") {
  
  message("Processing ", protein_name, ": ", length(peaks_gr), " peaks")
  
  ctl_counts <- bamCount(bam_control, peaks_gr, paired.end = paired_end_mode)
  trt_counts <- bamCount(bam_treat, peaks_gr, paired.end = paired_end_mode)
  
  ctl_cpm <- ctl_counts / (bam_libsize(bam_control) / 1e6)
  trt_cpm <- trt_counts / (bam_libsize(bam_treat) / 1e6)
  
  log2fc <- log2((trt_cpm + pseudocount_cpm) / (ctl_cpm + pseudocount_cpm))
  
  ctl_cpm_pos <- pmax(ctl_cpm, pseudocount_cpm)
  qbreaks <- unique(quantile(ctl_cpm, probs = seq(0, 1, length.out = n_bins + 1),
                             na.rm = TRUE, type = 7))
  bin_factor <- cut(ctl_cpm, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)
  
  df <- tibble(
    protein = protein_name,
    bin = bin_factor,
    ctl_cpm = ctl_cpm,
    trt_cpm = trt_cpm,
    log2fc = log2fc,
    log10_ctl_cpm = log10(ctl_cpm_pos)
  ) |> 
    filter(!is.na(bin))
  
  bin_summ <- df |>
    group_by(bin) |>
    summarise(
      n = dplyr::n(),
      median_log10_ctl_cpm = median(log10_ctl_cpm),
      median_log2fc = median(log2fc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(as.integer(bin))
  
  keep_bins <- bin_summ$bin[bin_summ$n >= min_peaks_per_bin]
  df <- df |> filter(bin %in% keep_bins)
  bin_summ <- bin_summ |> filter(bin %in% keep_bins)
  
  bin_summ$protein <- protein_name
  
  list(data = df, bins = bin_summ)
}

# =========================================================================
# PANEL A: log2FC vs CPM bins
# =========================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PANEL A: log2FC vs CPM bins (combined CTCF and RAD21)\n")
cat(rep("=", 70), "\n\n")

## Load peaks
ctcf_control_peaks <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
rad21_control_peaks <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
)

## BAM files
ctcf_bam_ctl <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_bam_trt <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"
rad21_bam_ctl <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_nodups_sorted.bam"
rad21_bam_trt <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_nodups_sorted.bam"

## Analyze with 100 bins
ctcf_res <- analyze_by_cpm_bins(ctcf_control_peaks, ctcf_bam_ctl, ctcf_bam_trt,
                                protein_name = "CTCF", n_bins = 100)
rad21_res <- analyze_by_cpm_bins(rad21_control_peaks, rad21_bam_ctl, rad21_bam_trt,
                                 protein_name = "RAD21", n_bins = 100)

## Combine and calculate percentiles
combined_bins <- bind_rows(ctcf_res$bins, rad21_res$bins) |>
  group_by(protein) |>
  mutate(percentile = (as.integer(bin) - 0.5) / max(as.integer(bin)) * 100) |>
  ungroup()

## Create plot with smaller text for compact display
protein_colors <- c("CTCF" = "#008B8B", "RAD21" = "#FF6B6B")

p_panelA <- ggplot(combined_bins, aes(x = percentile, y = median_log2fc, 
                                      color = protein, group = protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = protein_colors) +
  # moved labels up & left a bit
  annotate("text", x = 95, y = -1.0, label = "CTCF",
           color = protein_colors["CTCF"], fontface = "bold", size = 3,
           hjust = 0, vjust = 0.5) +
  annotate("text", x = 96, y = 0.35, label = "RAD21",
           color = protein_colors["RAD21"], fontface = "bold", size = 3,
           hjust = 0, vjust = 0.5) +
  labs(
    x = "CPM percentile",
    y = "log2FC(sorbitol/control)"
  ) +
  coord_cartesian(xlim = c(58, 100), ylim = c(-2.3, 0.5), clip = "off") +
  scale_x_continuous(limits = c(58, 100), expand = c(0, 0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(2, 20, 2, 2, "pt")
  )

cat("✓ Panel A plot created\n")

# =========================================================================
# PANEL B: CTCF retention by PWM score (boxplot) - MODIFIED
# =========================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PANEL B: CTCF retention by PWM score (MODIFIED: outliers removed, y-axis adjusted)\n")
cat(rep("=", 70), "\n\n")

## Load CTCF control peaks
ctcf_control <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
ctcf_control <- keepStandardChromosomes(ctcf_control, pruning.mode = "coarse")
cat("CTCF control peaks:", length(ctcf_control), "\n")

## Load DESeq2 results
ctcf_deseq <- readRDS("data/processed/cutntag/deseq2/diff_CTCF_counts.rds")

## Match and classify peaks
ctcf_overlaps <- findOverlaps(ctcf_control, ctcf_deseq)
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

cat("Classification:\n")
cat("  Retained:", sum(ctcf_matched$retention_category == "Retained"), "\n")
cat("  Lost:", sum(ctcf_matched$retention_category == "Lost"), "\n")

## Create summit windows
summit_positions <- start(ctcf_matched) + ctcf_matched$peak
window_size <- 50
flank_size <- (window_size / 2) - 1

ctcf_summits <- GRanges(
  seqnames = seqnames(ctcf_matched),
  ranges = IRanges(
    start = summit_positions - flank_size,
    end = summit_positions + flank_size + 1
  ),
  strand = strand(ctcf_matched)
)

mcols(ctcf_summits) <- mcols(ctcf_matched)

ctcf_summits_filtered <- ctcf_summits[
  ctcf_summits$retention_category %in% c("Retained", "Lost")
]

## Extract sequences
summit_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ctcf_summits_filtered)

## Get CTCF PWM
tryCatch({
  if (requireNamespace("JASPAR2020", quietly = TRUE)) {
    library(JASPAR2020)
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]] <- "CTCF"
    pfm_list <- getMatrixSet(JASPAR2020, opts)
    ctcf_pfm <- pfm_list[[1]]
    cat("Using JASPAR2020 database\n")
  } else {
    stop("JASPAR2020 not available")
  }
}, error = function(e) {
  stop("Could not load CTCF PWM from JASPAR. Please install JASPAR2020.")
})

## Calculate PWM scores (using searchSeq approach like original)
ctcf_pwm <- toPWM(ctcf_pfm)
pwm_scores <- sapply(summit_sequences, function(seq) {
  tryCatch({
    site_set <- searchSeq(ctcf_pwm, seq, min.score = "0%", strand = "*")
    if (length(site_set) > 0) {
      max(relScore(site_set))
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)
})

mcols(ctcf_summits_filtered)$pwm_score <- pwm_scores

## Create plotting data
plot_df_pwm <- as.data.frame(mcols(ctcf_summits_filtered)) |>
  dplyr::select(retention_category, pwm_score) |>
  filter(!is.na(pwm_score))

cat("Final data: Retained =", sum(plot_df_pwm$retention_category == "Retained"),
    ", Lost =", sum(plot_df_pwm$retention_category == "Lost"), "\n")

## MODIFICATION: Remove outliers using IQR method
remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  x[x >= (qnt[1] - H) & x <= (qnt[2] + H)]
}

plot_df_pwm_filtered <- plot_df_pwm |>
  group_by(retention_category) |>
  mutate(pwm_score_filtered = ifelse(
    pwm_score %in% remove_outliers(pwm_score),
    pwm_score,
    NA
  )) |>
  ungroup() |>
  filter(!is.na(pwm_score_filtered))

cat("After outlier removal: Retained =", 
    sum(plot_df_pwm_filtered$retention_category == "Retained"),
    ", Lost =", sum(plot_df_pwm_filtered$retention_category == "Lost"), "\n")

## MODIFICATION: Calculate new y-axis limits based on filtered data
y_range <- range(plot_df_pwm_filtered$pwm_score_filtered, na.rm = TRUE)
y_buffer <- diff(y_range) * 0.1
y_limits <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

cat("Y-axis range after outlier removal: [", 
    round(y_limits[1], 2), ",", round(y_limits[2], 2), "]\n")

## Statistical test
wilcox_result_pwm <- wilcox.test(
  pwm_score_filtered ~ retention_category, 
  data = plot_df_pwm_filtered
)
sig_stars_pwm <- case_when(
  wilcox_result_pwm$p.value < 0.001 ~ "***",
  wilcox_result_pwm$p.value < 0.01 ~ "**",
  wilcox_result_pwm$p.value < 0.05 ~ "*",
  TRUE ~ "ns"
)

## Create boxplot with adjusted y-axis and outliers removed
y_max_pwm <- max(plot_df_pwm_filtered$pwm_score_filtered)
y_range_pwm <- diff(range(plot_df_pwm_filtered$pwm_score_filtered))
bracket_y_pwm <- y_max_pwm + 0.05 * y_range_pwm

p_panelB <- ggplot(plot_df_pwm_filtered, 
                   aes(x = retention_category, 
                       y = pwm_score_filtered,
                       fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  annotate("segment", x = 1, xend = 2, y = bracket_y_pwm, yend = bracket_y_pwm,
           linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1, y = bracket_y_pwm, 
           yend = bracket_y_pwm - 0.01 * y_range_pwm, linewidth = 0.3) +
  annotate("segment", x = 2, xend = 2, y = bracket_y_pwm, 
           yend = bracket_y_pwm - 0.01 * y_range_pwm, linewidth = 0.3) +
  annotate("text", x = 1.5, y = bracket_y_pwm + 0.02 * y_range_pwm, 
           label = sig_stars_pwm, size = 4) +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    x = "CTCF Peaks",
    y = "CTCF PWM Score"
  ) +
  coord_cartesian(ylim = y_limits) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank()
  )

cat("✓ Panel B plot created (outliers removed, y-axis adjusted)\n")

# =========================================================================
# PANEL C: CTCF retention by RAD21 FC (boxplot) - MODIFIED
# =========================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PANEL C: CTCF retention by RAD21 FC (MODIFIED: boxplot instead of violin, outliers removed, y-axis adjusted)\n")
cat(rep("=", 70), "\n\n")

## Load RAD21 DESeq2 results (reuse ctcf_matched from Panel B)
rad21_deseq <- readRDS("data/processed/cutntag/deseq2/diff_RAD21_counts.rds")

## Match CTCF peaks to RAD21
ctcf_rad21_overlaps <- findOverlaps(ctcf_matched, rad21_deseq)

rad21_fc_by_ctcf <- mcols(rad21_deseq)$log2FoldChange[subjectHits(ctcf_rad21_overlaps)] |>
  split(queryHits(ctcf_rad21_overlaps)) |>
  sapply(mean, na.rm = TRUE)

ctcf_with_rad21_fc <- ctcf_matched[as.integer(names(rad21_fc_by_ctcf))]
mcols(ctcf_with_rad21_fc)$rad21_log2fc <- rad21_fc_by_ctcf

## Prepare plotting data
rad21_data <- as.data.frame(mcols(ctcf_with_rad21_fc)) |>
  filter(retention_category %in% c("Retained", "Lost"),
         !is.na(rad21_log2fc),
         is.finite(rad21_log2fc))

cat("RAD21 data: Retained =", sum(rad21_data$retention_category == "Retained"),
    ", Lost =", sum(rad21_data$retention_category == "Lost"), "\n")

## MODIFICATION: Remove outliers using IQR method
remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  x[x >= (qnt[1] - H) & x <= (qnt[2] + H)]
}

rad21_data_filtered <- rad21_data |>
  group_by(retention_category) |>
  mutate(rad21_log2fc_filtered = ifelse(
    rad21_log2fc %in% remove_outliers(rad21_log2fc),
    rad21_log2fc,
    NA
  )) |>
  ungroup() |>
  filter(!is.na(rad21_log2fc_filtered))

cat("After outlier removal: Retained =", 
    sum(rad21_data_filtered$retention_category == "Retained"),
    ", Lost =", sum(rad21_data_filtered$retention_category == "Lost"), "\n")

## MODIFICATION: Calculate new y-axis limits based on filtered data
y_range_c <- range(rad21_data_filtered$rad21_log2fc_filtered, na.rm = TRUE)
y_buffer_c <- diff(y_range_c) * 0.15
y_limits_c <- c(y_range_c[1] - y_buffer_c, y_range_c[2] + y_buffer_c)

cat("Y-axis range after outlier removal: [", 
    round(y_limits_c[1], 2), ",", round(y_limits_c[2], 2), "]\n")

## MODIFICATION: Create boxplot instead of violin plot with manual stats positioning
## Calculate statistics for manual annotation
wilcox_test_c <- wilcox.test(
  rad21_log2fc_filtered ~ retention_category,
  data = rad21_data_filtered
)

sig_label_c <- case_when(
  wilcox_test_c$p.value < 0.001 ~ "***",
  wilcox_test_c$p.value < 0.01 ~ "**",
  wilcox_test_c$p.value < 0.05 ~ "*",
  TRUE ~ "ns"
)

cat("Panel C Wilcoxon test p-value:", format.pval(wilcox_test_c$p.value), "\n")
cat("Significance label:", sig_label_c, "\n")

p_panelC <- ggplot(rad21_data_filtered, 
                   aes(x = retention_category, 
                       y = rad21_log2fc_filtered, 
                       fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.3) +
  # Manual statistics annotation positioned at y = 0.95-1.05
  annotate("segment", x = 1, xend = 2, y = 0.95, yend = 0.95,
           linewidth = 0.3, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 0.95, yend = 0.88,
           linewidth = 0.3, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 0.95, yend = 0.88,
           linewidth = 0.3, color = "black") +
  annotate("text", x = 1.5, y = 1.05, label = sig_label_c,
           size = 3, color = "black") +
  scale_fill_manual(
    values = c("Retained" = "#7FCDBB", "Lost" = "#E34A33"),
    name = "CTCF Peak Status"
  ) +
  labs(
    x = "CTCF Peaks",
    y = "RAD21 log2FC(sorb/ctrl)"
  ) +
  scale_y_continuous(limits = c(-2, 1.2)) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank()
  )

cat("✓ Panel C plot created (boxplot, outliers removed, y-axis adjusted)\n")

# =========================================================================
# PANEL D: CTCF promoter density
# =========================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("PANEL D: CTCF promoter enrichment analysis\n")
cat(rep("=", 70), "\n\n")

## Load CTCF control peaks (fresh for this panel)
ctcf_control_d <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
ctcf_control_d <- keepStandardChromosomes(ctcf_control_d, pruning.mode = "coarse")
seqlevelsStyle(ctcf_control_d) <- "UCSC"

## Load promoters
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
promoter_regions <- promoters(genes, upstream = 2000, downstream = 2000)

## Classify by promoter overlap
promoter_overlaps <- countOverlaps(ctcf_control_d, promoter_regions) > 0
mcols(ctcf_control_d)$overlaps_promoter <- promoter_overlaps
mcols(ctcf_control_d)$promoter_category <- ifelse(
  promoter_overlaps, "Promoter", "Non-promoter"
)

cat("Promoter peaks:", sum(promoter_overlaps), 
    "(", round(100 * mean(promoter_overlaps), 1), "%)\n")

## BAM files
ctcf_control_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_sorbitol_bam <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"

## Calculate log2FC
control_counts <- bamCount(ctcf_control_bam, ctcf_control_d, 
                           paired.end = "midpoint")
sorbitol_counts <- bamCount(ctcf_sorbitol_bam, ctcf_control_d, 
                            paired.end = "midpoint")

mcols(ctcf_control_d)$control_counts <- control_counts
mcols(ctcf_control_d)$sorbitol_counts <- sorbitol_counts
mcols(ctcf_control_d)$log2fc <- log2(
  (sorbitol_counts + 1) / (control_counts + 1)
)

## Filter by minimum signal
min_counts <- 10
ctcf_filtered_d <- ctcf_control_d[control_counts >= min_counts]

cat("Peaks after filtering (>=", min_counts, "counts):", length(ctcf_filtered_d), "\n")

## Prepare plotting data
peak_data <- as.data.frame(mcols(ctcf_filtered_d))

## Statistical test
wilcox_result_prom <- wilcox.test(
  log2fc ~ promoter_category, 
  data = peak_data,
  alternative = "two.sided"
)

sig_label_prom <- case_when(
  wilcox_result_prom$p.value < 0.0001 ~ "****",
  wilcox_result_prom$p.value < 0.001 ~ "***",
  wilcox_result_prom$p.value < 0.01 ~ "**",
  wilcox_result_prom$p.value < 0.05 ~ "*",
  TRUE ~ "ns"
)

## Create density plot with smaller features
p_panelD <- ggplot(peak_data, 
                   aes(x = log2fc, 
                       fill = promoter_category,
                       color = promoter_category)) +
  geom_density(alpha = 0.5, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
    name = ""
  ) +
  scale_color_manual(
    values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
    name = ""
  ) +
  labs(
    title = "Retained CTCF Peak Overlap",
    x = "log2FC(sorbitol/control)",
    y = ""
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "top",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, -5, 0),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2)
  )

cat("✓ Panel D plot created\n")

# =========================================================================
# Create combined multi-panel figure using plotgardener
# =========================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("Creating combined multi-panel figure (8.5\" width)\n")
cat(rep("=", 70), "\n\n")

## Define page dimensions to fit on letter-size paper
page_width <- 8.5     # Standard US letter width
page_height <- 3      # Adjust height to maintain readability

## Create PDF output
pdf("plots/ctcf_retentionFigure.pdf",
    width = page_width,
    height = page_height,
    useDingbats = FALSE)

## Initialize page
pageCreate(width = page_width, 
           height = page_height,
           showGuides = FALSE)

## Panel dimensions (smaller to fit 4 panels in 8.5")
panel_width <- 1.8
panel_height <- 1.8
spacing <- 0.2        # Smaller gaps between panels
margin_left <- 0.3
label_offset_x <- 0.15
label_offset_y <- 0.15

## Calculate panel x-positions
panel_x_positions <- c(
  margin_left,                                    # Panel A
  margin_left + panel_width + spacing,            # Panel B
  margin_left + 2 * (panel_width + spacing),      # Panel C
  margin_left + 3 * (panel_width + spacing)       # Panel D
)

panel_y <- 0.6  # Y position for all panels (from top)

## Plot Panel A
plotGG(
  plot = p_panelA,
  x = panel_x_positions[1],
  y = panel_y,
  width = panel_width,
  height = panel_height,
  just = c("left", "top")
)

plotText(
  label = "A",
  x = panel_x_positions[1] - label_offset_x,
  y = panel_y - label_offset_y,
  fontface = "bold",
  fontsize = 12
)

## Plot Panel B
plotGG(
  plot = p_panelB,
  x = panel_x_positions[2],
  y = panel_y,
  width = panel_width,
  height = panel_height,
  just = c("left", "top")
)

plotText(
  label = "B",
  x = panel_x_positions[2] - label_offset_x,
  y = panel_y - label_offset_y,
  fontface = "bold",
  fontsize = 12
)

## Plot Panel C
plotGG(
  plot = p_panelC,
  x = panel_x_positions[3],
  y = panel_y,
  width = panel_width,
  height = panel_height,
  just = c("left", "top")
)

plotText(
  label = "C",
  x = panel_x_positions[3] - label_offset_x,
  y = panel_y - label_offset_y,
  fontface = "bold",
  fontsize = 12
)

## Plot Panel D
plotGG(
  plot = p_panelD,
  x = panel_x_positions[4],
  y = panel_y,
  width = panel_width,
  height = panel_height,
  just = c("left", "top")
)

plotText(
  label = "D",
  x = panel_x_positions[4] - label_offset_x,
  y = panel_y - label_offset_y,
  fontface = "bold",
  fontsize = 12
)

## Close the PDF
dev.off()
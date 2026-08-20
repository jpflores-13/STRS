# FigureS8 ----------------------------------------------------------------
# Author:      JP Flores
# Date:        2026-04-23
# Project:     STRS
# Description: Four-panel CTCF retention supplementary figure. Panel A: median
#              log2FC vs. CPM percentile bin for CTCF and RAD21 CUT&Tag peaks.
#              Panel B: CTCF PWM score boxplot comparing retained vs. lost peaks.
#              Panel C: RAD21 log2FC boxplot at CTCF peak sites, stratified by
#              CTCF retention. Panel D: Density plot of CTCF CUT&Tag log2FC at
#              promoter vs. non-promoter peaks. All panels composed with
#              plotgardener.
# Input:       - CUT&Tag peaks:   data/processed/cutntag/output/peaks/*.narrowPeak
#              - DESeq2 results:  data/processed/cutntag/deseq2/diff_CTCF_counts.rds
#                                 data/processed/cutntag/deseq2/diff_RAD21_counts.rds
#              - Merged BAMs:     data/processed/cutntag/output/mergeAlign/*.bam
# Output:      - figures/FigureS8.pdf
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------

project_dir      <- "/work/users/j/p/jpflores/projects/STRS"
peaks_dir        <- file.path(project_dir, "data/processed/cutntag/output/peaks")
deseq_dir        <- file.path(project_dir, "data/processed/cutntag/deseq2")
bam_dir          <- file.path(project_dir, "data/processed/cutntag/output/mergeAlign")
output_pdf       <- file.path(project_dir, "figures/FigureS8.pdf")

## Peak files
ctcf_narrowpeak  <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak")
rad21_narrowpeak <- file.path(peaks_dir, "STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak")

## BAM files
ctcf_bam_ctl  <- file.path(bam_dir, "STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam")
ctcf_bam_trt  <- file.path(bam_dir, "STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam")
rad21_bam_ctl <- file.path(bam_dir, "STRS_HEK293_eGFP-YAP_RAD21_cont_0h_nodups_sorted.bam")
rad21_bam_trt <- file.path(bam_dir, "STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_nodups_sorted.bam")

## DESeq2 objects
ctcf_deseq_rds  <- file.path(deseq_dir, "diff_CTCF_counts.rds")
rad21_deseq_rds <- file.path(deseq_dir, "diff_RAD21_counts.rds")

## Analysis thresholds
padj_cutoff       <- 0.1
n_cpm_bins        <- 100
min_peaks_per_bin <- 20
pseudocount_cpm   <- 1e-6
min_bam_counts    <- 10    ## minimum control BAM counts to retain a peak (Panel D)
promoter_up       <- 2000  ## bp upstream for promoter definition
promoter_down     <- 2000  ## bp downstream for promoter definition
pwm_window        <- 50    ## bp window centered on CTCF summit for PWM scoring

## Page / panel geometry (inches)
page_width   <- 8.5
page_height  <- 3.0
panel_w      <- 1.8
panel_h      <- 1.8
col_gap      <- 0.2
top_margin   <- 0.6
left_margin  <- 0.3
label_off_x  <- 0.15
label_off_y  <- 0.15

## Derived x origins per panel
x_panels <- left_margin + (0:3) * (panel_w + col_gap)
y_panel  <- top_margin


# Libraries ---------------------------------------------------------------

library(bamsignals)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(JASPAR2020)
library(plotgardener)
library(plyranges)
library(Rsamtools)
library(TFBSTools)
library(tibble)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# Utility functions -------------------------------------------------------

read_narrowpeaks <- function(file) {
  peaks <- read.table(file, sep = "\t", header = FALSE,
                      col.names = c("seqnames", "start", "end", "name",
                                    "score", "strand", "signalValue",
                                    "pValue", "qValue", "peak"))
  GRanges(seqnames = peaks$seqnames,
          ranges   = IRanges(start = peaks$start + 1, end = peaks$end),
          strand   = "*",
          score = peaks$score, signalValue = peaks$signalValue,
          pValue = peaks$pValue, qValue = peaks$qValue, peak = peaks$peak)
}

bam_libsize <- function(bam) {
  sum(idxstatsBam(bam)$mapped, na.rm = TRUE)
}

remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  H   <- 1.5 * IQR(x, na.rm = TRUE)
  x[x >= (qnt[1] - H) & x <= (qnt[2] + H)]
}

sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

analyze_by_cpm_bins <- function(peaks_gr, bam_ctl, bam_trt,
                                protein_name, n_bins = n_cpm_bins,
                                min_n = min_peaks_per_bin) {
  message("Processing ", protein_name, ": ", length(peaks_gr), " peaks")
  ctl_counts <- bamCount(bam_ctl, peaks_gr, paired.end = "midpoint")
  trt_counts <- bamCount(bam_trt, peaks_gr, paired.end = "midpoint")
  ctl_cpm    <- ctl_counts / (bam_libsize(bam_ctl) / 1e6)
  trt_cpm    <- trt_counts / (bam_libsize(bam_trt) / 1e6)
  log2fc     <- log2((trt_cpm + pseudocount_cpm) / (ctl_cpm + pseudocount_cpm))
  qbreaks    <- unique(quantile(ctl_cpm, probs = seq(0, 1, length.out = n_bins + 1),
                                na.rm = TRUE, type = 7))
  bin_factor <- cut(ctl_cpm, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)
  tibble(protein       = protein_name,
         bin           = bin_factor,
         ctl_cpm       = ctl_cpm,
         log2fc        = log2fc,
         log10_ctl_cpm = log10(pmax(ctl_cpm, pseudocount_cpm))) |>
    dplyr::filter(!is.na(bin)) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(n                    = dplyr::n(),
                     median_log10_ctl_cpm = median(log10_ctl_cpm),
                     median_log2fc        = median(log2fc, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(n >= min_n) |>
    dplyr::mutate(percentile = (as.integer(bin) - 0.5) / max(as.integer(bin)) * 100,
                  protein    = protein_name)
}


# Load data ---------------------------------------------------------------

## CUT&Tag peaks
ctcf_ctrl_peaks  <- read_narrowpeaks(ctcf_narrowpeak)
rad21_ctrl_peaks <- read_narrowpeaks(rad21_narrowpeak)

## DESeq2 results
ctcf_deseq  <- readRDS(ctcf_deseq_rds)
rad21_deseq <- readRDS(rad21_deseq_rds)

## TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# Wrangle data ------------------------------------------------------------

## ---- Panel A: CPM-binned log2FC ----

combined_bins <- dplyr::bind_rows(
  analyze_by_cpm_bins(ctcf_ctrl_peaks,  ctcf_bam_ctl,  ctcf_bam_trt,  "CTCF"),
  analyze_by_cpm_bins(rad21_ctrl_peaks, rad21_bam_ctl, rad21_bam_trt, "RAD21")
)

## ---- Shared: annotate CTCF peaks with DESeq2 retention status ----
## Used by Panels B and C.

ctcf_ctrl_std <- keepStandardChromosomes(ctcf_ctrl_peaks, pruning.mode = "coarse")
ctcf_overlaps <- findOverlaps(ctcf_ctrl_std, ctcf_deseq)
ctcf_matched  <- ctcf_ctrl_std[queryHits(ctcf_overlaps)]

mcols(ctcf_matched)$log2FoldChange_CTCF <-
  mcols(ctcf_deseq)$log2FoldChange[subjectHits(ctcf_overlaps)]
mcols(ctcf_matched)$padj_CTCF <-
  mcols(ctcf_deseq)$padj[subjectHits(ctcf_overlaps)]

mcols(ctcf_matched)$retention_category <- dplyr::case_when(
  !is.na(ctcf_matched$padj_CTCF) &
    ctcf_matched$padj_CTCF < padj_cutoff &
    ctcf_matched$log2FoldChange_CTCF < 0 ~ "Lost",
  !is.na(ctcf_matched$padj_CTCF)         ~ "Retained",
  TRUE                                    ~ "Filtered"
)

## ---- Panel B: CTCF PWM score ----

flank_size   <- as.integer((pwm_window / 2) - 1)
summit_pos   <- start(ctcf_matched) + ctcf_matched$peak
ctcf_summits <- GRanges(
  seqnames = seqnames(ctcf_matched),
  ranges   = IRanges(summit_pos - flank_size, summit_pos + flank_size + 1L),
  strand   = strand(ctcf_matched)
)
mcols(ctcf_summits) <- mcols(ctcf_matched)
ctcf_summits_filt   <- ctcf_summits[
  ctcf_summits$retention_category %in% c("Retained", "Lost")
]

summit_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ctcf_summits_filt)
ctcf_pfm    <- getMatrixSet(JASPAR2020, list(species = 9606, name = "CTCF"))[[1]]
ctcf_pwm    <- toPWM(ctcf_pfm)

pwm_scores <- sapply(summit_seqs, function(seq) {
  sites <- searchSeq(ctcf_pwm, seq, min.score = "0%", strand = "*")
  if (length(sites) > 0) max(relScore(sites)) else NA_real_
})
mcols(ctcf_summits_filt)$pwm_score <- pwm_scores

pwm_df <- as.data.frame(mcols(ctcf_summits_filt)) |>
  dplyr::select(retention_category, pwm_score) |>
  dplyr::filter(!is.na(pwm_score)) |>
  dplyr::group_by(retention_category) |>
  dplyr::mutate(pwm_score = ifelse(pwm_score %in% remove_outliers(pwm_score),
                                   pwm_score, NA_real_)) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(pwm_score))

wilcox_pwm    <- wilcox.test(pwm_score ~ retention_category, data = pwm_df)
sig_pwm       <- sig_stars(wilcox_pwm$p.value)
y_max_pwm     <- max(pwm_df$pwm_score)
y_rng_pwm     <- diff(range(pwm_df$pwm_score))
bracket_y_pwm <- y_max_pwm + 0.05 * y_rng_pwm

## ---- Panel C: RAD21 log2FC at CTCF sites ----

ctcf_rad21_hits  <- findOverlaps(ctcf_matched, rad21_deseq)
rad21_fc_by_ctcf <- mcols(rad21_deseq)$log2FoldChange[subjectHits(ctcf_rad21_hits)] |>
  split(queryHits(ctcf_rad21_hits)) |>
  sapply(mean, na.rm = TRUE)

ctcf_with_rad21 <- ctcf_matched[as.integer(names(rad21_fc_by_ctcf))]
mcols(ctcf_with_rad21)$rad21_log2fc <- rad21_fc_by_ctcf

rad21_df <- as.data.frame(mcols(ctcf_with_rad21)) |>
  dplyr::filter(retention_category %in% c("Retained", "Lost"),
                !is.na(rad21_log2fc), is.finite(rad21_log2fc)) |>
  dplyr::group_by(retention_category) |>
  dplyr::mutate(rad21_log2fc = ifelse(rad21_log2fc %in% remove_outliers(rad21_log2fc),
                                      rad21_log2fc, NA_real_)) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(rad21_log2fc))

wilcox_rad21 <- wilcox.test(rad21_log2fc ~ retention_category, data = rad21_df)
sig_rad21    <- sig_stars(wilcox_rad21$p.value)

## ---- Panel D: CTCF promoter density ----

ctcf_prom <- keepStandardChromosomes(read_narrowpeaks(ctcf_narrowpeak),
                                     pruning.mode = "coarse")
seqlevelsStyle(ctcf_prom) <- "UCSC"

promoter_regions <- promoters(genes(txdb),
                              upstream   = promoter_up,
                              downstream = promoter_down)
mcols(ctcf_prom)$promoter_category <- ifelse(
  countOverlaps(ctcf_prom, promoter_regions) > 0, "Promoter", "Non-promoter"
)

ctl_cnt <- bamCount(ctcf_bam_ctl, ctcf_prom, paired.end = "midpoint")
trt_cnt <- bamCount(ctcf_bam_trt, ctcf_prom, paired.end = "midpoint")
mcols(ctcf_prom)$log2fc <- log2((trt_cnt + 1) / (ctl_cnt + 1))

prom_df <- as.data.frame(mcols(ctcf_prom[ctl_cnt >= min_bam_counts])) |>
  dplyr::select(promoter_category, log2fc)


# Visualization -----------------------------------------------------------

## Shared colors
retention_fills <- c("Retained" = "#7FCDBB", "Lost" = "#E34A33")
protein_colors  <- c("CTCF" = "#008B8B", "RAD21" = "#FF6B6B")

## Shared theme (matches S10)
theme_s10 <- theme_classic(base_size = 8) +
  theme(plot.title   = element_text(size = 8, face = "bold", hjust = 0.5,
                                    margin = margin(b = 3)),
        axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6),
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold"))

## ---- Panel A: CPM percentile line plot ----

p_a <- ggplot(combined_bins, aes(x = percentile, y = median_log2fc,
                                 color = protein, group = protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(values = protein_colors) +
  annotate("text", x = 95, y = -1.0, label = "CTCF",
           color = protein_colors["CTCF"],  fontface = "bold", size = 2.5, hjust = 0) +
  annotate("text", x = 96, y =  0.35, label = "RAD21",
           color = protein_colors["RAD21"], fontface = "bold", size = 2.5, hjust = 0) +
  labs(x = "CPM percentile", y = "log2FC(sorbitol/control)") +
  coord_cartesian(xlim = c(58, 100), ylim = c(-2.3, 0.5), clip = "off") +
  scale_x_continuous(limits = c(58, 100), expand = c(0, 0)) +
  theme_s10 +
  theme(legend.position = "none",
        plot.margin     = margin(2, 20, 2, 2, "pt"))

## ---- Panel B: CTCF PWM score boxplot ----

p_b <- ggplot(pwm_df, aes(x = retention_category, y = pwm_score,
                          fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  annotate("segment",
           x = 1, xend = 2, y = bracket_y_pwm, yend = bracket_y_pwm,
           linewidth = 0.3) +
  annotate("segment",
           x = 1, xend = 1,
           y = bracket_y_pwm, yend = bracket_y_pwm - 0.01 * y_rng_pwm,
           linewidth = 0.3) +
  annotate("segment",
           x = 2, xend = 2,
           y = bracket_y_pwm, yend = bracket_y_pwm - 0.01 * y_rng_pwm,
           linewidth = 0.3) +
  annotate("text", x = 1.5, y = bracket_y_pwm + 0.02 * y_rng_pwm,
           label = sig_pwm, size = 3) +
  scale_fill_manual(values = retention_fills) +
  labs(x = "CTCF Peaks", y = "CTCF PWM Score") +
  theme_s10 +
  theme(legend.position = "none",
        panel.grid      = element_blank())

## ---- Panel C: RAD21 log2FC by CTCF retention boxplot ----

p_c <- ggplot(rad21_df, aes(x = retention_category, y = rad21_log2fc,
                            fill = retention_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.3) +
  annotate("segment", x = 1, xend = 2, y = 0.95, yend = 0.95,
           linewidth = 0.3, color = "black") +
  annotate("segment", x = 1, xend = 1, y = 0.95, yend = 0.88,
           linewidth = 0.3, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 0.95, yend = 0.88,
           linewidth = 0.3, color = "black") +
  annotate("text", x = 1.5, y = 1.05, label = sig_rad21, size = 3) +
  scale_fill_manual(values = retention_fills) +
  labs(x = "CTCF Peaks", y = "RAD21 log2FC(sorb/ctrl)") +
  scale_y_continuous(limits = c(-2, 1.2)) +
  theme_s10 +
  theme(legend.position = "none",
        panel.grid      = element_blank())

## ---- Panel D: CTCF promoter density plot ----

p_d <- ggplot(prom_df, aes(x = log2fc, fill = promoter_category,
                           color = promoter_category)) +
  geom_density(alpha = 0.5, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  scale_fill_manual(values  = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
                    name = "") +
  scale_color_manual(values = c("Promoter" = "#F8766D", "Non-promoter" = "#619CFF"),
                     name = "") +
  labs(title = "Retained CTCF Peak Overlap",
       x = "log2FC(sorbitol/control)", y = "") +
  theme_bw(base_size = 8) +
  theme(legend.position   = "top",
        legend.key.size   = unit(0.3, "cm"),
        legend.text       = element_text(size = 7),
        legend.margin     = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, -5, 0),
        axis.line         = element_line(linewidth = 0.3),
        axis.ticks        = element_line(linewidth = 0.3),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(linewidth = 0.2))


# Save outputs ------------------------------------------------------------

pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

## ---- Panel labels ----
for (i in seq_along(x_panels)) {
  plotText(
    label         = LETTERS[i],
    fontsize      = 12,
    fontface      = "bold",
    x             = x_panels[i] - label_off_x,
    y             = y_panel     - label_off_y,
    just          = c("left", "top"),
    default.units = "inches"
  )
}

## ---- Panels ----
panels <- list(p_a, p_b, p_c, p_d)

for (i in seq_along(panels)) {
  plotGG(
    plot          = panels[[i]],
    x             = x_panels[i],
    y             = y_panel,
    width         = panel_w,
    height        = panel_h,
    just          = c("left", "top"),
    default.units = "inches"
  )
}

dev.off()

message("Saved: ", output_pdf)


# Session info ------------------------------------------------------------

sessionInfo()
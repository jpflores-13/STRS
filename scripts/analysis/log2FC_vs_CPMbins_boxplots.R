## Combined plot of log2FC vs CPM bins for CTCF and RAD21
## Shows median trend lines for both proteins on the same plot

library(GenomicRanges)
library(plyranges)
library(ggplot2)
library(dplyr)
library(bamsignals)
library(GenomeInfoDb)
library(Rsamtools)

# read narrowPeak ---------------------------------------------------------
read_narrowpeaks <- function(file) {
  gr <- plyranges::read_narrowpeaks(file)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

# library size from BAM (mapped reads) ------------------------------------
bam_libsize <- function(bam) {
  sum(idxstatsBam(bam)$mapped, na.rm = TRUE)
}

# quantify, bin by CPM (quantiles), compute log2FC, and build plot --------
analyze_by_cpm_bins <- function(peaks_gr,
                                bam_control,
                                bam_treat,
                                protein_name = "CTCF",
                                n_bins = 20,
                                pseudocount_cpm = 1e-6,
                                min_peaks_per_bin = 20,
                                paired_end_mode = "midpoint") {
  
  message("Peaks input: ", length(peaks_gr))
  
  # Count reads at control peaks in each condition
  ctl_counts <- bamCount(bam_control, peaks_gr, paired.end = paired_end_mode)
  trt_counts <- bamCount(bam_treat,   peaks_gr, paired.end = paired_end_mode)
  
  # CPM normalize
  ctl_cpm <- ctl_counts / (bam_libsize(bam_control) / 1e6)
  trt_cpm <- trt_counts / (bam_libsize(bam_treat)  / 1e6)
  
  # log2FC
  log2fc <- log2((trt_cpm + pseudocount_cpm) / (ctl_cpm + pseudocount_cpm))
  
  # Quantile bins on CONTROL CPM
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
  ) |> filter(!is.na(bin))
  
  bin_summ <- df |>
    group_by(bin) |>
    summarise(n = n(),
              median_log10_ctl_cpm = median(log10_ctl_cpm),
              median_log2fc = median(log2fc, na.rm = TRUE),
              .groups = "drop") |>
    arrange(as.integer(bin))
  
  # drop very small bins if needed
  keep_bins <- bin_summ$bin[bin_summ$n >= min_peaks_per_bin]
  df <- df |> filter(bin %in% keep_bins)
  bin_summ <- bin_summ |> filter(bin %in% keep_bins)
  
  # Add protein name to bin_summ for combined plotting
  bin_summ$protein <- protein_name
  
  list(data = df, bins = bin_summ)
}

# Inputs ------------------------------------------------------------------
ctcf_control_peaks <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
)
rad21_control_peaks <- read_narrowpeaks(
  "data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak"
)

ctcf_bam_ctl <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_nodups_sorted.bam"
ctcf_bam_trt <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_nodups_sorted.bam"
rad21_bam_ctl <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_cont_0h_nodups_sorted.bam"
rad21_bam_trt <- "data/processed/cutntag/output/mergeAlign/STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_nodups_sorted.bam"

# Run with 100 bins -------------------------------------------------------
ctcf_res  <- analyze_by_cpm_bins(ctcf_control_peaks, ctcf_bam_ctl, ctcf_bam_trt,
                                 protein_name = "CTCF", n_bins = 100)
rad21_res <- analyze_by_cpm_bins(rad21_control_peaks, rad21_bam_ctl, rad21_bam_trt,
                                 protein_name = "RAD21", n_bins = 100)

# Combine bin summaries for plotting --------------------------------------
combined_bins <- bind_rows(ctcf_res$bins, rad21_res$bins)

# Calculate percentile position for x-axis
combined_bins <- combined_bins |>
  group_by(protein) |>
  mutate(percentile = (as.integer(bin) - 0.5) / max(as.integer(bin)) * 100) |>
  ungroup()

# Combined plot -----------------------------------------------------------
# Define colors
protein_colors <- c("CTCF" = "#008B8B",   # Dark teal
                    "RAD21" = "#FF6B6B")   # Coral

# Get position for CTCF label (at end of line)
ctcf_end <- combined_bins |> filter(protein == "CTCF", percentile == max(percentile))

p_combined <- ggplot(combined_bins, aes(x = percentile, y = median_log2fc, 
                                        color = protein, group = protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.6) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = protein_colors) +
  # Add direct labels on the plot
  annotate("text", x = 98, y = -1.2, 
           label = "CTCF", color = protein_colors["CTCF"], 
           fontface = "bold", size = 4.5, hjust = 0, vjust = 0.5) +
  annotate("text", x = 100, y = 0.25, 
           label = "RAD21", color = protein_colors["RAD21"], 
           fontface = "bold", size = 4.5, hjust = 0.5, vjust = 0.5) +
  labs(
    x = "CPM percentile",
    y = "log2FC (sorbitol / control)"
  ) +
  coord_cartesian(xlim = c(58, 100), ylim = c(-2.3, 0.5), clip = "off") +
  scale_x_continuous(limits = c(58, 100), expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    plot.margin = margin(5.5, 40, 5.5, 5.5, "pt")
  )

# Save combined plot ------------------------------------------------------
dir.create("plots", showWarnings = FALSE)

pdf("plots/log2FC_vs_CPMbins_combined_100bins.pdf", width = 9, height = 6, useDingbats = FALSE)
print(p_combined)
dev.off()

# ENCODE ChIP-seq Overlap Scatter --------------------------------------------
# Author:      JP Flores
# Date:        2026-04-16
# Project:     STRS
# Description: For each ENCODE HEK293/HEK293T ChIP-seq experiment (TF and
#              histone marks), computes the percentage of peaks from six
#              query peak sets (CTCF retained/lost, RAD21 retained/lost,
#              ATAC gained/lost loop anchors) that overlap the ENCODE peaks.
#              Multiple experiments per target are collapsed by averaging.
#              Results are displayed as three QQ-style scatter plots
#              (one per comparison pair), colored by factor class.
# Input:       data/external/encode_chipseq/metadata_filtered.tsv
#              data/external/encode_chipseq/beds/*.bed.gz
#              data/processed/cutntag/output/peaks/*.narrowPeak
#              data/processed/cutntag/deseq2/diff_CTCF_counts.rds (unused directly)
#              data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds
# Output:      plots/encode_overlap_scatter.pdf
# -----------------------------------------------------------------------------


# Parameters ------------------------------------------------------------------

project_dir     <- "/work/users/j/p/jpflores/projects/STRS"
encode_dir      <- file.path(project_dir, "data/external/encode_chipseq")
beds_dir        <- file.path(encode_dir, "beds")
metadata_file   <- file.path(encode_dir, "metadata_filtered.tsv")
peaks_dir       <- file.path(project_dir, "data/processed/cutntag/output/peaks")
diff_loops_rds  <- file.path(project_dir, "data/processed/hic/diffLoops",
                             "diffLoops_eGFP-YAP_noDroso_10kb.rds")
output_pdf      <- file.path(project_dir, "plots/encode_overlap_scatter.pdf")

## Thresholds
padj_threshold  <- 0.1      # for gained/lost loop calls
peak_buffer     <- 1000     # bp buffer used in retained/lost peak classification
min_pct         <- 1        # minimum % overlap in at least one condition to label
top_n_label     <- 15       # number of top targets to label per panel

## Page layout
page_width      <- 14
page_height     <- 5.5
panel_width     <- 4.0
panel_height    <- 4.5


# Libraries -------------------------------------------------------------------

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(InteractionSet)
library(mariner)
library(plotgardener)
library(rtracklayer)
library(tidyr)


# Load data -------------------------------------------------------------------

## ENCODE metadata (filtered to IDR thresholded + replicated peaks, GRCh38)
metadata <- read.delim(metadata_file, header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)

## Differential loops
loops <- readRDS(diff_loops_rds) |>
  interactions() |>
  as.data.frame() |>
  as_ginteractions()

## CUT&Tag narrowPeak files — CTCF and RAD21 control + sorbitol
ctcf_control_peaks  <- file.path(peaks_dir,
                                 "STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak")   |> import()
ctcf_sorb_peaks     <- file.path(peaks_dir,
                                 "STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak") |> import()
rad21_control_peaks <- file.path(peaks_dir,
                                 "STRS_HEK293_eGFP-YAP_RAD21_cont_0h_peaks.narrowPeak")  |> import()
rad21_sorb_peaks    <- file.path(peaks_dir,
                                 "STRS_HEK293_eGFP-YAP_RAD21_sorbitol_1h_peaks.narrowPeak") |> import()


# Derive query peak sets ------------------------------------------------------

## Identify retained vs. lost peaks via countOverlaps with a bp buffer.
## Retained = control peaks overlapping >= 1 sorbitol peak (within buffer)
## Lost     = control peaks with no overlap in sorbitol
identify_peak_changes <- function(control_peaks, sorbitol_peaks,
                                  buffer = peak_buffer) {
  control_overlaps  <- countOverlaps(control_peaks,  sorbitol_peaks + buffer) > 0
  list(
    retained = control_peaks[control_overlaps],
    lost     = control_peaks[!control_overlaps]
  )
}

ctcf_changes  <- identify_peak_changes(ctcf_control_peaks,  ctcf_sorb_peaks)
rad21_changes <- identify_peak_changes(rad21_control_peaks, rad21_sorb_peaks)

## Extract unique loop anchors for gained and lost loops
extract_anchors <- function(gi) {
  unique(c(anchors(gi, "first"), anchors(gi, "second")))
}

gained_loops <- loops[loops$padj < padj_threshold & loops$log2FoldChange > 0]
lost_loops   <- loops[loops$padj < padj_threshold & loops$log2FoldChange < 0]

atac_gained_anchors <- extract_anchors(gained_loops)
atac_lost_anchors   <- extract_anchors(lost_loops)

## Named list of the 6 query peak sets
query_sets <- list(
  ctcf_retained  = ctcf_changes$retained,
  ctcf_lost      = ctcf_changes$lost,
  rad21_retained = rad21_changes$retained,
  rad21_lost     = rad21_changes$lost,
  atac_gained    = atac_gained_anchors,
  atac_lost      = atac_lost_anchors
)

## Report sizes
message("Query peak set sizes:")
for (nm in names(query_sets)) {
  message("  ", nm, ": ", length(query_sets[[nm]]))
}


# Compute overlaps ------------------------------------------------------------

## For each ENCODE BED file, compute % of each query set overlapping it.
## Option A: what fraction of YOUR peaks are hit by the ENCODE peaks?
##   sum(countOverlaps(your_peaks, encode_peaks) > 0) / length(your_peaks) * 100

bed_files <- list.files(beds_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)
message("Found ", length(bed_files), " ENCODE BED files")

compute_pct_overlap <- function(bed_file, query_sets) {
  encode_peaks <- tryCatch(
    import(bed_file, format = "narrowPeak"),
    error = function(e) {
      warning("Failed to import: ", bed_file, " — skipping")
      return(NULL)
    }
  )
  if (is.null(encode_peaks) || length(encode_peaks) == 0) return(NULL)
  
  accession <- gsub("\\.bed\\.gz$", "", basename(bed_file))
  
  pcts <- vapply(query_sets, function(qset) {
    sum(countOverlaps(qset, encode_peaks) > 0) / length(qset) * 100
  }, numeric(1))
  
  data.frame(
    file_accession = accession,
    as.list(pcts),
    check.names   = FALSE
  )
}

overlap_list <- lapply(bed_files, compute_pct_overlap,
                       query_sets = query_sets)

## Drop any NULLs (failed imports) and bind rows
overlap_df <- overlap_list[!sapply(overlap_list, is.null)] |>
  dplyr::bind_rows()

message("Overlap computed for ", nrow(overlap_df), " files")


# Join metadata and collapse by target ----------------------------------------

## Standardise column name (metadata uses backtick-unfriendly names)
meta_slim <- metadata |>
  dplyr::select(
    file_accession   = `File accession`,
    target_raw       = `Experiment target`,
    biosample        = `Biosample term name`,
    assay            = `Assay`
  ) |>
  ## Strip the "-human" suffix ENCODE appends to every target name
  dplyr::mutate(target = gsub("-human$", "", target_raw))

overlap_annotated <- overlap_df |>
  dplyr::left_join(meta_slim, by = "file_accession")

## Average across multiple experiments for the same target
overlap_collapsed <- overlap_annotated |>
  dplyr::group_by(target, assay) |>
  dplyr::summarise(
    dplyr::across(c(ctcf_retained, ctcf_lost,
                    rad21_retained, rad21_lost,
                    atac_gained, atac_lost),
                  mean, na.rm = TRUE),
    n_experiments = dplyr::n(),
    .groups = "drop"
  )

message("Unique targets after collapsing: ", nrow(overlap_collapsed))


# Assign factor classes -------------------------------------------------------

assign_factor_class <- function(target, assay) {
  dplyr::case_when(
    assay == "Histone ChIP-seq" &
      grepl("H3K27ac|H3K4me1|H3K4me3|H3K36me3|H3K9ac", target,
            ignore.case = TRUE)                                  ~ "Active mark",
    assay == "Histone ChIP-seq" &
      grepl("H3K27me3|H3K9me3|H3K9me2", target,
            ignore.case = TRUE)                                  ~ "Repressive mark",
    grepl("^CTCF$|^RAD21$|^SMC|^STAG", target,
          ignore.case = TRUE)                                    ~ "CTCF/Cohesin",
    grepl("^SP[0-9]|^KLF", target,
          ignore.case = TRUE)                                    ~ "SP/KLF",
    TRUE                                                         ~ "TF"
  )
}

overlap_collapsed <- overlap_collapsed |>
  dplyr::mutate(
    factor_class = assign_factor_class(target, assay),
    factor_class = factor(factor_class,
                          levels = c("Active mark", "Repressive mark",
                                     "CTCF/Cohesin", "SP/KLF", "TF"))
  )


# Build scatter panels --------------------------------------------------------

## Color palette — 5 factor classes, Okabe-Ito compatible
class_colors <- c(
  "Active mark"     = "#E69F00",
  "Repressive mark" = "#56B4E9",
  "CTCF/Cohesin"   = "#009E73",
  "SP/KLF"         = "#0072B2",
  "TF"             = "#CCCCCC"
)

## Point sizes — highlight biologically interesting classes
class_sizes <- c(
  "Active mark"     = 2.2,
  "Repressive mark" = 2.2,
  "CTCF/Cohesin"   = 2.2,
  "SP/KLF"         = 2.2,
  "TF"             = 1.4
)

## Alpha values
class_alpha <- c(
  "Active mark"     = 0.9,
  "Repressive mark" = 0.9,
  "CTCF/Cohesin"   = 0.9,
  "SP/KLF"         = 0.9,
  "TF"             = 0.3
)

make_scatter_panel <- function(df, x_col, y_col, x_label, y_label, title) {
  
  df <- df |>
    dplyr::mutate(
      x        = .data[[x_col]],
      y        = .data[[y_col]],
      diag_dist = y - x,
      abs_dist  = abs(diag_dist)
    )
  
  ax_max <- max(c(df$x, df$y), na.rm = TRUE) * 1.05
  
  ## Label top N by absolute diagonal distance, excluding the grey TF cloud
  ## unless they happen to be among the most extreme
  top_motifs <- df |>
    dplyr::filter(pmax(x, y) >= min_pct) |>
    dplyr::slice_max(order_by = abs_dist, n = top_n_label, with_ties = FALSE)
  
  df_bg <- df |> dplyr::filter(factor_class == "TF")
  df_fg <- df |> dplyr::filter(factor_class != "TF")
  
  ggplot(df, aes(x = x, y = y)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey60", linetype = "dashed", linewidth = 0.4) +
    ## Grey background cloud first
    geom_point(data  = df_bg,
               aes(color = factor_class, size = factor_class,
                   alpha = factor_class)) +
    ## Colored foreground points on top
    geom_point(data  = df_fg,
               aes(color = factor_class, size = factor_class,
                   alpha = factor_class)) +
    ggrepel::geom_text_repel(
      data          = top_motifs,
      aes(label = target, color = factor_class),
      size          = 2.5,
      fontface      = "bold",
      max.overlaps  = 30,
      segment.size  = 0.3,
      segment.color = "grey50",
      show.legend   = FALSE
    ) +
    scale_color_manual(values = class_colors, name = "Factor class") +
    scale_size_manual(values  = class_sizes,  name = "Factor class") +
    scale_alpha_manual(values = class_alpha,  name = "Factor class") +
    coord_fixed(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(
      title = title,
      x     = paste0("% overlap  [", x_label, "]"),
      y     = paste0("% overlap  [", y_label, "]")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(size = 10, face = "bold",
                                     hjust = 0.5, margin = margin(b = 4)),
      axis.title      = element_text(size = 9),
      axis.text       = element_text(size = 8),
      legend.position = "right",
      legend.title    = element_text(size = 8, face = "bold"),
      legend.text     = element_text(size = 7)
    )
}

p_ctcf  <- make_scatter_panel(
  df      = overlap_collapsed,
  x_col   = "ctcf_lost",
  y_col   = "ctcf_retained",
  x_label = "CTCF lost peaks",
  y_label = "CTCF retained peaks",
  title   = "CTCF: Retained vs. Lost peaks"
)

p_rad21 <- make_scatter_panel(
  df      = overlap_collapsed,
  x_col   = "rad21_lost",
  y_col   = "rad21_retained",
  x_label = "RAD21 lost peaks",
  y_label = "RAD21 retained peaks",
  title   = "RAD21: Retained vs. Lost peaks"
)

p_atac  <- make_scatter_panel(
  df      = overlap_collapsed,
  x_col   = "atac_lost",
  y_col   = "atac_gained",
  x_label = "ATAC lost anchors",
  y_label = "ATAC gained anchors",
  title   = "ATAC: Gained vs. Lost loop anchors"
)


# Save outputs ----------------------------------------------------------------

## Save overlap table for downstream inspection
write.table(
  overlap_collapsed,
  file      = file.path(encode_dir, "overlap_results_collapsed.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

## Save figure via plotgardener
pdf(output_pdf, width = page_width, height = page_height)

pageCreate(width = page_width, height = page_height, showGuides = FALSE)

x_positions <- c(0.5,
                 0.5 + panel_width + 0.6,
                 0.5 + (panel_width + 0.6) * 2)

panel_labels <- c("A", "B", "C")
panels       <- list(p_ctcf, p_rad21, p_atac)

for (i in seq_along(panels)) {
  plotText(
    label    = panel_labels[i],
    x        = x_positions[i] - 0.1,
    y        = 0.35,
    fontsize = 14,
    fontface = "bold"
  )
  plotGG(
    plot   = panels[[i]],
    x      = x_positions[i],
    y      = 0.5,
    width  = panel_width,
    height = panel_height,
    just   = c("left", "top")
  )
}

dev.off()

message("Figure saved to: ", output_pdf)


# Session info ----------------------------------------------------------------

sessionInfo()
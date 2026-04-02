# ##############################################################################
# filename:    hic_timecourse_apa_enrichment.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Calculates loop enrichment scores from pre-computed APA matrices
#              using center pixel vs top-left + bottom-right corner background;
#              outputs per-timepoint enrichment tables for gained and lost loops
# ##############################################################################

# Libraries ----
library(mariner)
library(tidyverse)
library(SummarizedExperiment)

# Parameters ----
apa_dir         <- "data/processed/hic/normalizedAPA"
output_dir      <- "data/processed/hic/apa_enrichment_matrices"
output_all_rds  <- "data/processed/hic/apa_loop_enrichment_all_original_bg.rds"
output_summ_rds <- "data/processed/hic/apa_loop_enrichment_summary_original_bg.rds"
timepoints      <- c("0h", "10m", "30m", "1h", "3h", "6h", "12h", "24h")
apa_buffer      <- 10

# Data import ----
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time_order      <- setNames(seq_along(timepoints), timepoints)
timepoint_hours <- c(0, 10/60, 30/60, 1, 3, 6, 12, 24)
names(timepoint_hours) <- timepoints

# Analysis ----
calculate_enrichment <- function(loop_type = c("gained", "lost"),
                                 buffer = apa_buffer,
                                 prefix = "") {
  loop_type <- match.arg(loop_type)

  pattern   <- paste0("tc_", loop_type, "Loops.*normalized\\.rds")
  apa_files <- list.files(apa_dir, full.names = TRUE, pattern = pattern)

  if (length(apa_files) == 0) stop("No APA matrix files found for ", loop_type, " loops")

  apa_files <- apa_files[order(sapply(apa_files, function(x) {
    tp <- str_extract(x, "0h|10m|30m|1h|3h|6h|12h|24h")
    time_order[tp]
  }))]

  message("Found ", length(apa_files), " APA matrix files for ", loop_type, " loops")

  apa_matrices    <- vector("list", length(apa_files))
  timepoint_names <- character(length(apa_files))

  for (i in seq_along(apa_files)) {
    tp               <- str_extract(apa_files[i], "0h|10m|30m|1h|3h|6h|12h|24h")
    timepoint_names[i] <- tp
    message("  Loading timepoint: ", tp)
    apa_matrices[[i]] <- readRDS(apa_files[i])
  }

  saveRDS(apa_matrices, file.path(output_dir, paste0(prefix, "_apa_matrices.rds")))

  enrichment_scores <- numeric(length(apa_matrices))

  for (i in seq_along(apa_matrices)) {
    mat          <- apa_matrices[[i]]
    center_idx   <- buffer + 1
    fg_val       <- mat[center_idx, center_idx]
    mat_size     <- nrow(mat)
    br_start     <- mat_size - 4
    bg_vals      <- c(as.vector(mat[1:5, 1:5]),
                      as.vector(mat[br_start:mat_size, br_start:mat_size]))
    bg_vals      <- bg_vals[!is.na(bg_vals)]

    enrichment_scores[i] <- if (length(bg_vals) > 0) fg_val / median(bg_vals) else NA
  }

  tibble(
    timepoint = factor(timepoint_names, levels = timepoints),
    hours     = timepoint_hours[timepoint_names],
    file      = paste0(loop_type, "_", timepoint_names, "_apa_matrix"),
    score     = enrichment_scores
  )
}

gained_enrichment           <- calculate_enrichment("gained", prefix = "gained_apa_bg")
gained_enrichment$loop_type <- "Gained"

lost_enrichment           <- calculate_enrichment("lost",   prefix = "lost_apa_bg")
lost_enrichment$loop_type <- "Lost"

all_enrichment <- bind_rows(gained_enrichment, lost_enrichment)

enrichment_summary <- all_enrichment |>
  group_by(loop_type, timepoint, hours) |>
  summarize(
    mean_score   = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score     = sd(score, na.rm = TRUE),
    se_score     = sd_score / sqrt(sum(!is.na(score))),
    lower_ci     = mean_score - 1.96 * se_score,
    upper_ci     = mean_score + 1.96 * se_score,
    n            = sum(!is.na(score)),
    .groups      = "drop"
  )

# Save outputs ----
saveRDS(all_enrichment,    output_all_rds)
saveRDS(enrichment_summary, output_summ_rds)

message("Saved: ", output_all_rds)
message("Saved: ", output_summ_rds)

sessionInfo()

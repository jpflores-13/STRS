# ##############################################################################
# filename:    strengthened_weakened_loop_analysis.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: Overlap between differential loop calls and pre-existing untreated
#              loops; quantifies de novo vs strengthened gained loops and anchor
#              sharing; barplots of loop-level and anchor-level overlap
# ##############################################################################

# Libraries ----
library(InteractionSet)
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(ggplot2)
library(glue)
library(data.table)
library(mariner)
library(purrr)
library(tibble)
library(readr)

# Parameters ----
loops_dir         <- "data/processed/hic/loops"
diff_loops_rds    <- "data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds"
output_overlap    <- "plots/strengthened_weakened_loop_analysis.pdf"
output_anchor     <- "plots/gained_anchor_sharing.pdf"
padj_threshold    <- 0.1
lfc_threshold     <- 0
anchor_width      <- 10000L

# Data import ----
loop_files <- list.files(loops_dir, pattern = "control_sip-loops_5kbLoops.txt$",
                          full.names = TRUE)
stopifnot(length(loop_files) > 0)

untreated_loops <- map_dfr(loop_files, fread, .id = "dataset") |>
  dplyr::rename(chr1    = chromosome1, start1 = x1, end1 = x2,
                chr2    = chromosome2, start2 = y1, end2 = y2) |>
  dplyr::distinct(chr1, start1, end1, chr2, start2, end2, .keep_all = TRUE) |>
  dplyr::relocate(dataset, .after = dplyr::last_col()) |>
  as_ginteractions()

diff_loopCounts <- readRDS(diff_loops_rds)
deseq_results   <- rowData(diff_loopCounts)

# Analysis ----
loop_categories <- deseq_results |>
  as.data.frame() |>
  mutate(category = case_when(
    log2FoldChange >  lfc_threshold & padj < padj_threshold ~ "Gained",
    log2FoldChange < -lfc_threshold & padj < padj_threshold ~ "Lost",
    TRUE ~ "Static"
  ))

rowData(diff_loopCounts)$category <- loop_categories$category

all_gi    <- interactions(diff_loopCounts)
gained_gi <- all_gi[loop_categories$category == "Gained"]
lost_gi   <- all_gi[loop_categories$category == "Lost"]
static_gi <- all_gi[loop_categories$category == "Static"]

calculate_overlap_percentage <- function(query_gi, subject_gi, category_name) {
  if (length(query_gi) == 0L) {
    return(data.frame(category = category_name, total_loops = 0L,
                      overlapping_loops = 0L, overlap_percentage = 0))
  }
  overlaps <- findOverlaps(query_gi, subject_gi, type = "any")
  n_query  <- length(query_gi)
  n_over   <- length(unique(queryHits(overlaps)))
  data.frame(category = category_name, total_loops = n_query,
             overlapping_loops = n_over, overlap_percentage = n_over / n_query * 100)
}

overlap_results <- list(
  calculate_overlap_percentage(gained_gi, untreated_loops, "Gained"),
  calculate_overlap_percentage(lost_gi,   untreated_loops, "Lost"),
  calculate_overlap_percentage(static_gi, untreated_loops, "Static")
) |> bind_rows()

gained_overlap_pct <- overlap_results$overlap_percentage[overlap_results$category == "Gained"]
de_novo_pct        <- 100 - gained_overlap_pct

extract_anchors_std <- function(gi, width = anchor_width) {
  if (length(gi) == 0L) return(GRanges())
  resize(c(anchors(gi, "first"), anchors(gi, "second")), width = width, fix = "center")
}

gAnch_std <- extract_anchors_std(gained_gi)
uAnch_std <- extract_anchors_std(untreated_loops)

if (length(gAnch_std) > 0L && length(uAnch_std) > 0L) {
  h <- findOverlaps(gAnch_std, uAnch_std, ignore.strand = TRUE)
  anchor_shared_pct <- length(unique(queryHits(h))) / length(gAnch_std) * 100
} else {
  h <- NULL
  anchor_shared_pct <- 0
}

sentence_numbers <- tibble(
  metric = c("gained_overlap_pct", "gained_de_novo_pct", "anchor_shared_pct"),
  value  = c(gained_overlap_pct, de_novo_pct, anchor_shared_pct)
)

cat(glue("Among sorbitol-induced loops, {round(de_novo_pct, 1)}% are de novo, ",
         "{round(gained_overlap_pct, 1)}% overlap pre-existing loops. ",
         "{round(anchor_shared_pct, 1)}% of sorbitol-induced loop anchors ",
         "are shared with untreated anchors.\n\n"))
print(sentence_numbers)

# Visualization ----
category_colors <- c("Gained" = "#F8766D", "Lost" = "#619CFF", "Static" = "#999999")

overlap_plot <- overlap_results |>
  mutate(category = factor(category, levels = c("Gained", "Lost", "Static"))) |>
  ggplot(aes(x = category, y = overlap_percentage, fill = category)) +
  geom_col(width = 0.7, linewidth = 0.5) +
  geom_text(aes(label = glue("{round(overlap_percentage, 1)}%")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = category_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(y = "% Overlap with Pre-Existing Loops", x = "", title = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 9))

anchor_total      <- length(gAnch_std)
anchor_shared_n   <- if (!is.null(h)) length(unique(queryHits(h))) else 0L
anchor_not_shared <- anchor_total - anchor_shared_n

anchor_share_df <- tibble(
  status  = c("Shared with untreated", "Not shared"),
  count   = c(anchor_shared_n, anchor_not_shared)
) |>
  mutate(percent = count / anchor_total * 100,
         status  = factor(status, levels = c("Shared with untreated", "Not shared")))

anchor_plot <- ggplot(anchor_share_df, aes(x = status, y = percent, fill = status)) +
  geom_col(width = 0.7, linewidth = 0.5) +
  geom_text(aes(label = glue("{round(percent, 1)}% (n={count})")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Shared with untreated" = "#6AA84F", "Not shared" = "#B7B7B7")) +
  ylim(0, 100) +
  labs(title = "Sorbitol-Induced Loop Anchors",
       y = "% of gained loop anchors", x = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"))

# Save outputs ----
ggsave(output_overlap, overlap_plot, width = 8, height = 6)
ggsave(output_anchor,  anchor_plot,  width = 6, height = 5)

sessionInfo()

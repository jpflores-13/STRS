## Overlap Analysis: Untreated vs Differential Loop Calls
## Compare original loop calls with differential analysis categories
## Includes: (1) de novo vs strengthened among gained loops
##           (2) anchor-level sharing percentages + bar plot

# Load required libraries ------------------------------------------------
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

# Load and combine control loop data -------------------------------------

loop_files <- list.files(
  path = "data/processed/hic/loops",
  pattern = "control_sip-loops_5kbLoops.txt$",
  full.names = TRUE
)

stopifnot(length(loop_files) > 0)

untreated_loops <- map_dfr(loop_files, fread, .id = "dataset") |>
  rename(chr1 = chromosome1,
         start1 = x1,
         end1 = x2,
         chr2 = chromosome2,
         start2 = y1,
         end2 = y2) |>
  distinct(chr1, start1, end1, chr2, start2, end2, .keep_all = TRUE) |>
  relocate(dataset, .after = dplyr::last_col()) |>
  as_ginteractions()

# Load differential loop results -----------------------------------------
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds")
deseq_results <- rowData(diff_loopCounts)

# Categorize loops -------------------------------------------------------
padj_threshold <- 0.1
lfc_threshold  <- 0

loop_categories <- deseq_results |>
  as.data.frame() |>
  mutate(
    category = case_when(
      log2FoldChange >  lfc_threshold & padj < padj_threshold ~ "Gained",
      log2FoldChange < -lfc_threshold & padj < padj_threshold ~ "Lost",
      TRUE ~ "Static"
    )
  )

rowData(diff_loopCounts)$category <- loop_categories$category

all_gi    <- interactions(diff_loopCounts)
gained_gi <- all_gi[loop_categories$category == "Gained"]
lost_gi   <- all_gi[loop_categories$category == "Lost"]
static_gi <- all_gi[loop_categories$category == "Static"]

# Loop-level overlap function --------------------------------------------
calculate_overlap_percentage <- function(query_gi, subject_gi, category_name) {
  
  if (length(query_gi) == 0L) {
    return(data.frame(
      category = category_name,
      total_loops = 0L,
      overlapping_loops = 0L,
      overlap_percentage = 0
    ))
  }
  
  overlaps <- findOverlaps(query_gi, subject_gi, type = "any")
  n_query <- length(query_gi)
  n_overlaps <- length(unique(queryHits(overlaps)))
  overlap_percentage <- (n_overlaps / n_query) * 100
  
  data.frame(
    category = category_name,
    total_loops = n_query,
    overlapping_loops = n_overlaps,
    overlap_percentage = overlap_percentage
  )
}

overlap_results <- list(
  calculate_overlap_percentage(gained_gi, untreated_loops, "Gained"),
  calculate_overlap_percentage(lost_gi, untreated_loops, "Lost"),
  calculate_overlap_percentage(static_gi, untreated_loops, "Static")
) |>
  bind_rows()

# Manuscript-ready numbers ------------------------------------------------
gained_overlap_pct <- as.numeric(overlap_results$overlap_percentage[overlap_results$category == "Gained"])
de_novo_pct <- 100 - gained_overlap_pct

anchor_width <- 10000L

extract_anchors_std <- function(gi, width = anchor_width) {
  if (length(gi) == 0L) return(GRanges())
  a1 <- anchors(gi, type = "first")
  a2 <- anchors(gi, type = "second")
  resize(c(a1, a2), width = width, fix = "center")
}

gAnch_std <- extract_anchors_std(gained_gi)
uAnch_std <- extract_anchors_std(untreated_loops)

if (length(gAnch_std) > 0L && length(uAnch_std) > 0L) {
  h <- findOverlaps(gAnch_std, uAnch_std, ignore.strand = TRUE)
  anchor_shared_pct <- length(unique(queryHits(h))) / length(gAnch_std) * 100
} else {
  anchor_shared_pct <- 0
}

sentence_numbers <- tibble(
  metric = c("gained_overlap_pct", "gained_de_novo_pct", "anchor_shared_pct"),
  value  = c(gained_overlap_pct, de_novo_pct, anchor_shared_pct)
)

summary_text <- glue(
  "Among sorbitol-induced loops, {round(de_novo_pct, 1)}% are de novo, ",
  "{round(gained_overlap_pct, 1)}% overlap pre-existing loops. ",
  "{round(anchor_shared_pct, 1)}% of sorbitol-induced loop anchors are shared with untreated anchors."
)

# Loop-level overlap plot -------------------------------------------------
category_colors <- c("Gained"="#F8766D","Lost"="#619CFF","Static"="#999999")

(overlap_plot <- overlap_results |>
    mutate(category = factor(category, levels = c("Gained","Lost","Static"))) |>
    ggplot(aes(x = category, y = overlap_percentage, fill = category)) +
    geom_col(width = 0.7, linewidth = 0.5) +
    geom_text(aes(label = glue("{round(overlap_percentage, 1)}%")),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(y="% Overlap with Pre-Existing Loops", x="", title="") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9)))

ggsave("plots/strengthened_weakened_loop_analysis.pdf", overlap_plot, width=8, height=6)

# Anchor-level sharing bar plot -------------------------------------------
anchor_total <- length(gAnch_std)
anchor_shared_n <- if (exists("h")) length(unique(queryHits(h))) else 0L
anchor_not_shared_n <- anchor_total - anchor_shared_n

anchor_share_df <- tibble(
  status = c("Shared with untreated", "Not shared"),
  count  = c(anchor_shared_n, anchor_not_shared_n)
) |>
  mutate(percent = count / anchor_total * 100)

anchor_colors <- c("Shared with untreated"="#6AA84F","Not shared"="#B7B7B7")

gained_anchor_share_plot <- ggplot(anchor_share_df,
                                   aes(x=status, y=percent, fill=status)) +
  geom_col(width=0.7, linewidth=0.5) +
  geom_text(aes(label=glue("{round(percent,1)}% (n={count})")),
            vjust=-0.3, size=3.5, fontface="bold") +
  scale_fill_manual(values=anchor_colors) +
  ylim(0,100) +
  labs(title="Sorbitol-Induced Loop Anchors",
       y="% of gained loop anchors", x="") +
  theme_classic() +
  theme(legend.position="none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"))

ggsave("plots/gained_anchor_sharing.pdf", gained_anchor_share_plot, width=6, height=5)

cat("\n===========================================\n")
cat("MANUSCRIPT STATEMENT:\n")
cat(summary_text)
cat("\n===========================================\n\n")

print("Summary of values:")
print(sentence_numbers)
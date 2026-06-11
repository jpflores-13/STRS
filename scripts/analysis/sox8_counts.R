# ##############################################################################
# filename:    sox8_counts.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: Quick plot of SOX8 normalized counts across untreated, sorbitol,
#              and sorbitol+auxin conditions using DESeq2::plotCounts on the
#              AuxinRAD21 exonic dds object. Builds dds from scratch from the
#              existing tximport RDS and samplesheet so no prior analysis script
#              needs to have been run.
# input:       data/processed/rna/AuxinRAD21/output/quant/YAPP_HCT116_mAID2-RAD21_1_tximport.rds
#              data/processed/rna/AuxinRAD21/output/YAPP_HCT116_mAID2-RAD21_1_RNApipeSamplesheet.txt
# output:      plots/sox8_counts.pdf
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(data.table)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(stringr)


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"
plots_dir   <- file.path(project_dir, "plots")

tximport_rds     <- file.path(project_dir,
                              "data/processed/rna/AuxinRAD21/output/quant",
                              "YAPP_HCT116_mAID2-RAD21_1_tximport.rds")
samplesheet_file <- file.path(project_dir,
                              "data/processed/rna/AuxinRAD21/output",
                              "YAPP_HCT116_mAID2-RAD21_1_RNApipeSamplesheet.txt")

## SOX8 Ensembl ID (version suffix matched by grep)
sox8_ensembl <- "ENSG00000005513"

## Colors
color_untreated <- "#999999"
color_sorbitol  <- "#F8766D"
color_auxin     <- "#009E73"

## Output
output_pdf <- file.path(plots_dir, "sox8_counts.pdf")


# Load and prepare data -------------------------------------------------------

se <- readRDS(tximport_rds)

## Assign column names alphabetically — same logic as other AuxinRAD21 scripts
samplesheet   <- fread(samplesheet_file) |> as.data.frame()
# sn_ordered    <- sort(samplesheet$sn)
sn_ordered    <- samplesheet$sn
stopifnot(ncol(se$counts) == length(sn_ordered))
colnames(se$counts)    <- sn_ordered
colnames(se$abundance) <- sn_ordered
colnames(se$length)    <- sn_ordered

sample_names  <- colnames(se$counts)
treatment_vec <- dplyr::case_when(
  str_detect(sample_names, "untreated") ~ "untreated",
  str_detect(sample_names, "sorbitol")  ~ "sorbitol",
  str_detect(sample_names, "auxin")     ~ "auxin"
)
rep_vec <- str_match(sample_names, "_(\\d+)_\\d+$")[, 2]

coldata <- data.frame(
  Treatment = factor(treatment_vec,
                     levels = c("untreated", "sorbitol", "auxin")),
  Bio_Rep   = factor(rep_vec),
  row.names = sample_names
)


# Build DESeq2 object ---------------------------------------------------------

dds <- DESeqDataSetFromTximport(se, colData = coldata,
                                design = ~Bio_Rep + Treatment)
dds <- DESeq(dds)


# Plot SOX8 counts ------------------------------------------------------------

sox8_id <- grep(sox8_ensembl, rownames(dds), value = TRUE)

if (length(sox8_id) == 0L) {
  stop("SOX8 (", sox8_ensembl, ") not found in rownames. ",
       "Check that the gene passed filtering.")
}

message("Plotting: ", sox8_id)

sox8_df <- plotCounts(dds,
                      gene       = sox8_id,
                      intgroup   = "Treatment",
                      returnData = TRUE)

(p <- ggplot(sox8_df, aes(x = Treatment, y = count, color = Treatment)) +
  geom_jitter(width = 0.15, size = 3) +
  scale_y_log10() +
  scale_color_manual(values = c(
    "untreated" = color_untreated,
    "sorbitol"  = color_sorbitol,
    "auxin"     = color_auxin
  )) +
  labs(
    title = "SOX8",
    x     = NULL,
    y     = "Normalized counts (log10)"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "italic", size = 11),
    legend.position = "none",
    axis.text       = element_text(size = 9, color = "black"),
    axis.title      = element_text(size = 9),
    axis.line       = element_line(linewidth = 0.3),
    axis.ticks      = element_line(linewidth = 0.3)
  ))


# Save ------------------------------------------------------------------------

ggsave(output_pdf, p, width = 3.5, height = 4, units = "in")
message("Saved: ", output_pdf)


# Session info ----------------------------------------------------------------

sessionInfo()
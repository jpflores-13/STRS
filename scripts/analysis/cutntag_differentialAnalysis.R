## Differential CUT&Tag peak analysis with DESeq2

library(DESeq2)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(apeglm)
library(pheatmap)

# Load peak counts --------------------------------------------------------

peaks <- fread("data/processed/cutntag/output/peaks/STRS_HEK293_eGFP-YAP_1_peakCounts.tsv")

cat("Loaded peak counts:\n")
cat("Number of peaks:", nrow(peaks), "\n")

## Convert to GRanges object
peaks <- GRanges(seqnames = Rle(peaks$chr),
                 ranges = IRanges(start = peaks$start, end = peaks$stop),
                 mcols = peaks[, -c(1:3)])

colnames(mcols(peaks)) <- colnames(mcols(peaks)) |>
  str_remove(pattern = "mcols.")

cat("Converted to GRanges object with", length(peaks), "peaks\n")

# Helper Functions --------------------------------------------------------

## Function to run DESeq2 for a single protein
run_deseq_protein <- function(peaks, protein_name) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing:", protein_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Create matrix for countData -------------------------------------------
  
  ## Extract columns for this protein
  pattern <- paste0(".*", protein_name, ".*")
  m <- mcols(peaks)[, grep(pattern, colnames(mcols(peaks)))] |>
    as.matrix()
  
  cat("Number of samples:", ncol(m), "\n")
  
  # Construct colData/metadata ----------------------------------------------
  
  ## String split the colnames
  colData <- as.data.frame(
    do.call(rbind, strsplit(colnames(m), "_")),
    stringsAsFactors = TRUE
  )
  rownames(colData) <- colnames(m)
  colnames(colData) <- c("Project", "Cell_Type", "Genotype", "Protein",
                         "Treatment", "Time", "Replicate", "Extra")
  
  ## Keep relevant columns and clean up
  colData <- colData |>
    dplyr::select(Project, Cell_Type, Protein, Treatment, Replicate) |>
    mutate(
      Replicate = factor(Replicate),
      Treatment = factor(Treatment, levels = c("cont", "sorbitol"))
    )
  
  ## Make sure sample names match
  stopifnot(all(colnames(m) == rownames(colData)))
  
  cat("\nColumn data:\n")
  print(colData)
  
  # Run DESeq2 --------------------------------------------------------------
  
  dds <- DESeqDataSetFromMatrix(
    countData = m,
    colData = colData,
    design = ~ Replicate + Treatment
  )
  
  ## Disable DESeq's default normalization
  sizeFactors(dds) <- rep(1, ncol(dds))
  
  ## Run DESeq2
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  pdf(paste0("plots/", protein_name, "_dispersion.pdf"), width = 6, height = 5)
  plotDispEsts(dds)
  dev.off()
  
  # QC via data visualization -----------------------------------------------
  
  ## Plot PCA
  pca_plot <- plotPCA(vst(dds), intgroup = "Treatment") +
    theme(aspect.ratio = 1) +
    ggtitle(paste(protein_name, "- PCA"))
  
  ggsave(paste0("plots/", protein_name, "_PCA.pdf"),
         plot = pca_plot,
         width = 6,
         height = 6)
  
  ## Transform counts for hierarchical clustering
  rld <- rlog(dds, blind = TRUE)
  
  ## Extract the rlog matrix
  rld_mat <- assay(rld)
  
  ## Compute pairwise correlation values
  rld_cor <- cor(rld_mat)
  
  ## Plot heatmap
  pdf(paste0("plots/", protein_name, "_correlation_heatmap.pdf"),
      width = 6,
      height = 6)
  pheatmap(rld_cor,
           main = paste(protein_name, "- Sample Correlation"))
  dev.off()
  
  # Get results -------------------------------------------------------------
  
  res <- results(dds)
  cat("\nDESeq2 results:\n")
  cat("Result names:", resultsNames(dds), "\n")
  print(summary(res))
  
  ## Apply LFC shrinkage
  res <- lfcShrink(dds, coef = "Treatment_sorbitol_vs_cont", type = "apeglm")
  
  cat("\nAfter LFC shrinkage:\n")
  print(summary(res))
  
  return(list(
    dds = dds,
    results = res,
    peaks = peaks
  ))
}


# CUSTOM PLOTTING FUNCTIONS -----------------------------------------------

## Function to create MA plot matching your style
create_ma_plot <- function(res, protein_name) {
  
  ## Convert to data frame
  res_df <- as.data.frame(res) |>
    mutate(
      isDE = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Increased",
        log2FoldChange < -1 & padj < 0.05 ~ "Decreased",
        TRUE ~ "Not significant"
      )
    ) |>
    arrange(isDE)
  
  ## Count differential peaks
  n_up <- sum(res_df$isDE == "Increased", na.rm = TRUE)
  n_down <- sum(res_df$isDE == "Decreased", na.rm = TRUE)
  
  ## Create MA plot
  res_gg <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = isDE)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey40") +
    scale_color_manual(values = c(
      "Increased" = "#F8766D",
      "Decreased" = "#619CFF",
      "Not significant" = "grey80"
    )) +
    ylim(c(-4, 4)) +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    labs(
      title = paste(protein_name),
      y = "log2FoldChange(sorbitol/control)",
      x = "mean of normalized counts"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    ) +
    annotate(geom = "text",
             label = "Increased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.5,
             y = 3.5,
             color = "#F8766D",
             fontface = "bold") +
    annotate(geom = "text",
             label = paste0("n = ", n_up),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.5,
             y = 3.1,
             color = "#F8766D") +
    annotate(geom = "text",
             label = "Decreased",
             x = max(res_df$baseMean, na.rm = TRUE) * 0.5,
             y = -3.5,
             color = "#619CFF",
             fontface = "bold") +
    annotate(geom = "text",
             label = paste0("n = ", n_down),
             x = max(res_df$baseMean, na.rm = TRUE) * 0.5,
             y = -3.9,
             color = "#619CFF")
  
  return(list(plot = res_gg, data = res_df))
}

## Function to create density plot
create_density_plot <- function(res_df) {
  
  ## Filter to only include differential peaks (gained or lost)
  res_diff <- res_df |>
    filter(isDE %in% c("Increased", "Decreased"))
  
  ggplot(res_diff, aes(y = log2FoldChange, fill = isDE)) +
    geom_density(
      color = NA,              # No outline
      alpha = 0.25
    ) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey40") +
    scale_fill_manual(values = c(
      "Increased" = "#F8766D",
      "Decreased" = "#619CFF"
    )) +
    ylim(c(-4, 4)) +
    xlim(c(0, 1.5)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line.x = element_blank()
    )
}

# Analysis ----------------------------------------------------------------

## Create output directories
dir.create("data/processed/cutntag/deseq2", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

## Get list of proteins
proteins <- c("CTCF", "H3K27ac", "RAD21", "YAP1")

## Run DESeq2 analysis for each protein
deseq_results <- list()
ma_plot_data <- list()

for (prot in proteins) {
  
  ## Run analysis
  deseq_results[[prot]] <- run_deseq_protein(peaks, prot)
  
  ## Use DESeq2's built-in plotMA function
  pdf(paste0("plots/", prot, "_DESeq2_MAplot.pdf"), width = 8, height = 6)
  plotMA(deseq_results[[prot]]$results, 
         ylim = c(-4, 4),
         main = paste("MA Plot:", prot, "(sorbitol vs control)"))
  dev.off()
  
  ## Extract data for summary statistics
  res_df <- as.data.frame(deseq_results[[prot]]$results) |>
    mutate(
      isDE = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Increased",
        log2FoldChange < -1 & padj < 0.05 ~ "Decreased",
        TRUE ~ "Not significant"
      )
    )
  
  ma_plot_data[[prot]] <- res_df
  
  ## ---- Custom MA plot outputs -------------------------------------------
  # Create MA plot (returns ggplot + data.frame)
  custom_ma <- create_ma_plot(deseq_results[[prot]]$results, protein_name = prot)
  
  # Save the custom MA plot
  ggsave(filename = paste0("plots/", prot, "_custom_MA.pdf"),
         plot = custom_ma$plot, width = 7, height = 6, useDingbats = FALSE)
  
  # Save the data used for the plot
  write.csv(custom_ma$data,
            file = paste0("data/processed/cutntag/deseq2/", prot, "_custom_MA_data.csv"),
            row.names = FALSE)
  
  # Density plot (right-side strip)
  dens_plot <- create_density_plot(custom_ma$data)
  ggsave(filename = paste0("plots/", prot, "_custom_density.pdf"),
         plot = dens_plot, width = 2.0, height = 6, useDingbats = FALSE)
  
  ## Concatenate peaks and DESeq results
  mcols(deseq_results[[prot]]$peaks) <- cbind(
    mcols(deseq_results[[prot]]$peaks),
    deseq_results[[prot]]$results
  )
  
  diff_peakCounts <- deseq_results[[prot]]$peaks |>
    keepStandardChromosomes(pruning.mode = "coarse")
  
  seqlevelsStyle(diff_peakCounts) <- "UCSC"
  
  ## Save as .rds
  saveRDS(diff_peakCounts,
          file = paste0("data/processed/cutntag/deseq2/diff_", prot, "_counts.rds"))
}

# Save all results --------------------------------------------------------

saveRDS(deseq_results, "data/processed/cutntag/deseq2/deseq2_results_all.rds")

# Summary Statistics ------------------------------------------------------

cat("\n\n")
cat(rep("=", 70), "\n", sep = "")
cat("SUMMARY OF DIFFERENTIAL PEAK ANALYSIS\n")
cat(rep("=", 70), "\n", sep = "")

summary_table <- lapply(names(ma_plot_data), function(prot) {
  ma_plot_data[[prot]] |>
    group_by(isDE) |>
    summarise(count = n(), .groups = "drop") |>
    mutate(protein = prot)
}) |>
  bind_rows() |>
  pivot_wider(names_from = isDE, values_from = count, values_fill = 0) |>
  mutate(
    total = rowSums(across(where(is.numeric))),
    pct_increased = round(Increased / total * 100, 1),
    pct_decreased = round(Decreased / total * 100, 1)
  )

print(summary_table)

## Save summary table
write.csv(summary_table,
          "data/processed/cutntag/deseq2/deseq2_summary_stats.csv",
          row.names = FALSE)

cat("\n\nAnalysis complete!\n")
cat("Individual plots and results saved for each protein\n")
cat("Summary saved to: data/processed/cutntag/deseq2/deseq2_summary_stats.csv\n")

## Print session info
sessionInfo()

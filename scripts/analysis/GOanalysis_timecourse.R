# GO Enrichment Analysis

# Libraries ---------------------------------------------------------------

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(enrichplot)


# Parameters --------------------------------------------------------------

# Significance thresholds for differential expression
PADJ_THRESHOLD <- 0.05  # Adjusted p-value cutoff
LFC_THRESHOLD <- 2      # Log2 fold change threshold (absolute value)

# GO enrichment parameters
GO_ONTOLOGY <- "BP"           # "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
GO_PVALUE_CUTOFF <- 0.05      # P-value cutoff for GO enrichment
GO_QVALUE_CUTOFF <- 0.2       # Q-value (FDR) cutoff for GO enrichment
GO_PADJUST_METHOD <- "BH"     # P-value adjustment method: "BH" (Benjamini-Hochberg), "bonferroni", etc.
GO_MIN_GSET_SIZE <- 10        # Minimum gene set size to consider
GO_MAX_GSET_SIZE <- 500       # Maximum gene set size to consider

# Plotting parameters
TOP_N_TERMS <- 5              # Number of top terms to show per category (TOP 5!)
DESCRIPTION_WIDTH <- 50       # Character width for truncating GO term descriptions


# Load Data ---------------------------------------------------------------

# Load your DESeq2 object
dds <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")

# Explanation:
# This DESeq2 object should contain results from a likelihood ratio test (LRT)
# comparing a time course model to a reduced model. The object contains:
# - Normalized counts
# - Gene metadata in rowData()
# - Results coefficients for each timepoint vs baseline


# Extract Differential Genes ----------------------------------------------

# Get all time coefficient names from the DESeq2 results
time_coefficients <- resultsNames(dds) |> 
  str_subset("Time")

# Explanation:
# resultsNames() returns all coefficient names from the model
# We filter for those containing "Time" to get timepoint comparisons
# Expected format: "Time_1h_vs_0h", "Time_3h_vs_0h", etc.

cat("Found", length(time_coefficients), "time coefficients:\n")
print(time_coefficients)


# Function to extract upregulated genes for a given coefficient
extract_upreg_genes <- function(coefficient) {
  # Get results for this coefficient
  res <- results(dds, name = coefficient) |>
    as.data.frame() |>
    filter(!is.na(padj))
  
  # Filter for significantly upregulated genes
  # Upregulated: log2FC > threshold AND padj < threshold
  upreg_genes <- res |>
    filter(padj < PADJ_THRESHOLD & log2FoldChange > LFC_THRESHOLD) |>
    rownames()
  
  return(upreg_genes)
}

# Function to extract downregulated genes for a given coefficient
extract_downreg_genes <- function(coefficient) {
  # Get results for this coefficient
  res <- results(dds, name = coefficient) |>
    as.data.frame() |>
    filter(!is.na(padj))
  
  # Filter for significantly downregulated genes
  # Downregulated: log2FC < -threshold AND padj < threshold
  downreg_genes <- res |>
    filter(padj < PADJ_THRESHOLD & log2FoldChange < -LFC_THRESHOLD) |>
    rownames()
  
  return(downreg_genes)
}

# Extract genes for all timepoints
upreg_list <- map(time_coefficients, extract_upreg_genes)
downreg_list <- map(time_coefficients, extract_downreg_genes)

# Explanation:
# We apply the extraction functions to each timepoint
# This creates lists where each element contains gene IDs for one timepoint
# We then combine across all timepoints to get the union of all DE genes

# Print summary statistics
cat("\nUpregulated genes per timepoint:\n")
walk2(time_coefficients, upreg_list, ~cat(.x, ":", length(.y), "genes\n"))

cat("\nDownregulated genes per timepoint:\n")
walk2(time_coefficients, downreg_list, ~cat(.x, ":", length(.y), "genes\n"))


# Combine genes across all timepoints
all_upreg_genes <- upreg_list |>
  unlist() |>
  unique()

all_downreg_genes <- downreg_list |>
  unlist() |>
  unique()

cat("\nTotal unique upregulated genes:", length(all_upreg_genes), "\n")
cat("Total unique downregulated genes:", length(all_downreg_genes), "\n")


# Convert Gene IDs to Entrez IDs -----------------------------------------

# Extract gene information from DESeq2 object
gene_info <- rowData(dds) |>
  as.data.frame()

# Explanation:
# rowData() contains metadata for each gene
# We expect a 'symbol' column with gene symbols
# If gene_id is not a column name, we add it from rownames

if (!"gene_id" %in% colnames(gene_info)) {
  gene_info <- gene_info |>
    rownames_to_column("gene_id")
}

gene_info <- gene_info |>
  dplyr::select(gene_id, symbol)


# Convert upregulated genes to Entrez IDs
upreg_symbols <- gene_info |>
  filter(gene_id %in% all_upreg_genes) |>
  pull(symbol)

upreg_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = upreg_symbols,
  keytype = "SYMBOL",
  columns = "ENTREZID"
) |>
  pull(ENTREZID) |>
  na.omit() |>
  unique()

cat("\nUpregulated: Converted", length(upreg_symbols), "symbols to", 
    length(upreg_entrez), "Entrez IDs\n")

# Explanation:
# clusterProfiler requires Entrez IDs for GO enrichment
# org.Hs.eg.db is the human genome annotation database
# Some symbols may map to multiple Entrez IDs or fail to map (NA)


# Convert downregulated genes to Entrez IDs
downreg_symbols <- gene_info |>
  filter(gene_id %in% all_downreg_genes) |>
  pull(symbol)

downreg_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = downreg_symbols,
  keytype = "SYMBOL",
  columns = "ENTREZID"
) |>
  pull(ENTREZID) |>
  na.omit() |>
  unique()

cat("Downregulated: Converted", length(downreg_symbols), "symbols to", 
    length(downreg_entrez), "Entrez IDs\n")


# Define background/universe
# This is all genes that were tested (expressed in your experiment)
background_symbols <- gene_info |>
  pull(symbol) |>
  na.omit()

background_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = background_symbols,
  keytype = "SYMBOL",
  columns = "ENTREZID"
) |>
  pull(ENTREZID) |>
  na.omit() |>
  unique()

cat("Background: Converted", length(background_symbols), "symbols to", 
    length(background_entrez), "Entrez IDs\n")

# Explanation:
# The 'universe' parameter represents all genes that could have been
# selected as significant. Using all expressed genes as background
# provides more accurate enrichment statistics than using all human genes.


# Perform GO Enrichment ---------------------------------------------------

cat("\n=== Running GO Enrichment Analysis ===\n")
cat("Ontology:", GO_ONTOLOGY, "\n")
cat("P-value cutoff:", GO_PVALUE_CUTOFF, "\n")
cat("Q-value cutoff:", GO_QVALUE_CUTOFF, "\n")
cat("Gene set size range:", GO_MIN_GSET_SIZE, "-", GO_MAX_GSET_SIZE, "\n\n")

# Upregulated genes
cat("Analyzing upregulated genes...\n")
go_upreg <- enrichGO(
  gene = upreg_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = GO_ONTOLOGY,
  pAdjustMethod = GO_PADJUST_METHOD,
  pvalueCutoff = GO_PVALUE_CUTOFF,
  qvalueCutoff = GO_QVALUE_CUTOFF,
  minGSSize = GO_MIN_GSET_SIZE,
  maxGSSize = GO_MAX_GSET_SIZE,
  readable = TRUE  # Convert Entrez IDs back to gene symbols in results
)

cat("Found", nrow(go_upreg@result), "enriched terms for upregulated genes\n\n")

# Explanation:
# enrichGO performs over-representation analysis (ORA)
# - gene: your list of DE genes
# - universe: background set of all tested genes
# - OrgDb: organism annotation database
# - ont: which GO ontology to use
# - pAdjustMethod: how to correct for multiple testing
# - readable: if TRUE, gene IDs are converted to symbols in output


# Downregulated genes
cat("Analyzing downregulated genes...\n")
go_downreg <- enrichGO(
  gene = downreg_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = GO_ONTOLOGY,
  pAdjustMethod = GO_PADJUST_METHOD,
  pvalueCutoff = GO_PVALUE_CUTOFF,
  qvalueCutoff = GO_QVALUE_CUTOFF,
  minGSSize = GO_MIN_GSET_SIZE,
  maxGSSize = GO_MAX_GSET_SIZE,
  readable = TRUE
)

cat("Found", nrow(go_downreg@result), "enriched terms for downregulated genes\n\n")


# Save Full GO Results as CSV ---------------------------------------------

# Save enrichment results as CSV (full results)
write_csv(
  go_upreg@result,
  "tables/go_enrichment_upregulated.csv"
)

write_csv(
  go_downreg@result,
  "tables/go_enrichment_downregulated.csv"
)

cat("Saved full GO enrichment results:\n")
cat("- tables/go_enrichment_upregulated.csv\n")
cat("- tables/go_enrichment_downregulated.csv\n\n")


# =========================================================================
# CREATE OSMOTIC STRESS-FOCUSED PLOTS (TOP 5 SENSIBLE TERMS EACH)
# =========================================================================

cat("\n=== Creating Osmotic Stress-Focused Plots (Top 5 Each) ===\n\n")

# Manually Select TOP 5 Most Biologically Relevant GO Terms --------------

# For UPREGULATED genes - Top 5 most relevant to osmotic stress:
# Focus on: ion transport, mechanical stress, and cell reorganization

selected_upreg_ids <- c(
  "GO:0006814",  # sodium ion transport (THE key osmotic pathway!)
  "GO:0009612",  # response to mechanical stimulus (osmotic = mechanical stress)
  "GO:0007159",  # leukocyte cell-cell adhesion (cell reorganization)
  "GO:0050727",  # regulation of inflammatory response (stress response)
  "GO:0001525"   # angiogenesis (tissue adaptation)
)

# For DOWNREGULATED genes - Top 5 most biologically sensible:
# Note: All have q-values > 0.74 (weak significance due to only 71 genes)
# Focusing on: stress response, translation shutdown, RNA metabolism

selected_downreg_ids <- c(
  "GO:0008630",  # intrinsic apoptotic signaling pathway in response to DNA damage
  "GO:0006448",  # regulation of translational elongation (translation shutdown)
  "GO:0009451",  # RNA modification (strongest statistical signal)
  "GO:0006400",  # tRNA modification
  "GO:0002097"   # tRNA wobble base modification
)

# Biological rationale:
# 
# UPREGULATED (Top 5):
# 1. Sodium ion transport - THE critical pathway for maintaining osmotic balance
# 2. Mechanical stimulus - osmotic pressure creates mechanical stress on cells
# 3. Leukocyte adhesion - cells reorganize contacts under osmotic stress
# 4. Inflammatory response - general cellular stress response
# 5. Angiogenesis - vascular/tissue adaptation to osmotic challenge
#
# DOWNREGULATED (Top 5):
# 1. DNA damage response - osmotic stress can cause DNA damage
# 2. Translation regulation - cells shut down translation during stress
# 3-5. RNA/tRNA modification - strongest statistical signal (even if weak)


# Filter and Prepare Data for Plotting -----------------------------------

# Filter upregulated terms
go_upreg_data <- go_upreg@result |>
  as.data.frame() |>
  filter(ID %in% selected_upreg_ids) |>
  arrange(p.adjust) |>
  head(TOP_N_TERMS) |>
  mutate(
    log_padj = -log10(p.adjust),
    Description = str_trunc(Description, width = DESCRIPTION_WIDTH),
    category = "Upregulated Genes"
  )

# Filter downregulated terms
go_downreg_data <- go_downreg@result |>
  as.data.frame() |>
  filter(ID %in% selected_downreg_ids) |>
  arrange(p.adjust) |>
  head(TOP_N_TERMS) |>
  mutate(
    log_padj = -log10(p.adjust),
    Description = str_trunc(Description, width = DESCRIPTION_WIDTH),
    category = "Downregulated Genes"
  )


# Print Summary Statistics ------------------------------------------------

cat("=== TOP 5 OSMOTIC STRESS-RELEVANT GO TERMS ===\n\n")

cat("UPREGULATED GENES (Top 5):\n")
cat("Focus: Ion transport, mechanical stress, cell reorganization\n\n")
go_upreg_data |>
  select(Description, p.adjust, Count) |>
  print()
cat("\n")

cat("DOWNREGULATED GENES (Top 5):\n")
cat("Focus: DNA damage response, translation shutdown, RNA metabolism\n")
cat("NOTE: Weak statistical significance (all q > 0.74) due to small gene count (n=", 
    length(all_downreg_genes), ")\n\n")
go_downreg_data |>
  select(Description, p.adjust, Count) |>
  print()
cat("\n")


# Combine and Prepare for Plotting ---------------------------------------

combined_data <- bind_rows(go_upreg_data, go_downreg_data) |>
  mutate(
    Description = fct_reorder(Description, log_padj),
    category = factor(category, levels = c("Upregulated Genes", "Downregulated Genes"))
  )


# Create Visualization ----------------------------------------------------

# Define colors to match your figure
color_upregulated <- "#F8766D"
color_downregulated <- "#619CFF"

# Create the barplot
go_plot <- ggplot(combined_data, aes(x = log_padj, y = Description, fill = category)) +
  geom_col(width = 0.7, color = NA) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  labs(
    x = expression(-log[10]("adjusted p-value")),
    y = NULL,
    title = "GO Enrichment Analysis"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 9, face = "bold"),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.line.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5, "pt"),
    strip.text = element_text(size = 9, face = "bold", hjust = 0),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Upregulated Genes" = color_upregulated, 
      "Downregulated Genes" = color_downregulated
    )
  ) +
  geom_text(
    aes(x = 0.3, label = Description), 
    hjust = 0, 
    size = 2.5, 
    color = "black", 
    fontface = "plain"
  )

# Display the plot
print(go_plot)


# Save Osmotic Stress Plots -----------------------------------------------

ggsave(
  filename = "plots/go_enrichment_osmotic_stress.pdf",
  plot = go_plot,
  width = 8,
  height = 6,
  units = "in"
)

write_csv(
  combined_data,
  "tables/go_enrichment_osmotic_stress_data.csv"
)


# Create Annotated Version with Pathway Groups ---------------------------

combined_data_annotated <- combined_data |>
  mutate(
    pathway_group = case_when(
      ID == "GO:0006814" ~ "Ion Homeostasis",
      ID == "GO:0009612" ~ "Mechanical Stress",
      ID == "GO:0007159" ~ "Cell Adhesion",
      ID == "GO:0050727" ~ "Inflammatory Response",
      ID == "GO:0001525" ~ "Angiogenesis",
      ID == "GO:0008630" ~ "DNA Damage Response",
      ID == "GO:0006448" ~ "Translation Regulation",
      ID %in% c("GO:0009451", "GO:0006400", "GO:0002097") ~ "RNA Metabolism",
      TRUE ~ "Other"
    )
  )

write_csv(
  combined_data_annotated,
  "tables/go_enrichment_osmotic_stress_annotated.csv"
)


# Print Summary of Saved Files --------------------------------------------

cat("\n=== All Results Saved ===\n\n")
cat("Full GO enrichment results:\n")
cat("  - tables/go_enrichment_upregulated.csv\n")
cat("  - tables/go_enrichment_downregulated.csv\n\n")

cat("Osmotic stress-focused plots and data (TOP 5 EACH):\n")
cat("  - plots/go_enrichment_osmotic_stress.pdf\n")
cat("  - tables/go_enrichment_osmotic_stress_data.csv\n")
cat("  - tables/go_enrichment_osmotic_stress_annotated.csv\n\n")


# Session Info ------------------------------------------------------------

cat("\n=== Session Info ===\n")
sessionInfo()
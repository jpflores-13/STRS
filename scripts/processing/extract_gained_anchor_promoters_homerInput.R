# extract_gained_anchor_promoters.R
# Identify genes whose promoters overlap with gained loop anchors

library(GenomicRanges)
library(plyranges)
library(dplyr)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

cat("========================================\n")
cat("Creating HOMER Input Files\n")
cat("========================================\n\n")

# Load differential loops -------------------------------------------------
cat("Loading differential loops...\n")
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

# Extract GAINED loops ----------------------------------------------------
gained_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                        diff_loopCounts$log2FoldChange > 0)]

cat("Number of gained loops:", length(gained_loops), "\n")

# Get all loop anchors from gained loops ----------------------------------
gained_anchors <- c(anchors(gained_loops, type = "first"),
                    anchors(gained_loops, type = "second")) |>
  unique()

# CRITICAL: Ensure chromosome names have 'chr' prefix
seqlevelsStyle(gained_anchors) <- "UCSC"

cat("Number of unique gained anchors:", length(gained_anchors), "\n")

# Define promoters by genes -----------------------------------------------
cat("Defining gene promoters...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)

# Create promoters from genes
promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 500)

# CRITICAL: Ensure chromosome names match
seqlevelsStyle(promoters_gr) <- "UCSC"

# Gene IDs are in the names
promoters_gr$gene_id <- names(promoters_gr)

# Keep only standard chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
promoters_gr <- promoters_gr[seqnames(promoters_gr) %in% standard_chrs]
gained_anchors <- gained_anchors[seqnames(gained_anchors) %in% standard_chrs]

cat("Promoters on standard chromosomes:", length(promoters_gr), "\n")
cat("Gained anchors on standard chromosomes:", length(gained_anchors), "\n")

# Find promoters overlapping gained anchors -------------------------------
cat("Finding promoters at gained anchors...\n")
overlaps <- findOverlaps(promoters_gr, gained_anchors)
gained_promoter_indices <- unique(queryHits(overlaps))

gained_promoters <- promoters_gr[gained_promoter_indices]
all_promoters <- promoters_gr

cat("Gained anchor promoters found:", length(gained_promoters), "\n")

# Export for HOMER --------------------------------------------------------
cat("\nPreparing BED files...\n")
dir.create("data/processed/cutntag/homer_input", 
           showWarnings = FALSE, 
           recursive = TRUE)

# Convert to data frames with UNIQUE IDs
gained_df <- as.data.frame(gained_promoters) |>
  dplyr::select(seqnames, start, end, gene_id, strand) |>
  mutate(
    start = start - 1,  # 0-based for BED
    score = 0,
    seqnames = as.character(seqnames),
    strand = as.character(strand),
    # Create unique ID based on coordinates
    name = paste0(seqnames, ":", start, "-", end)
  )

all_df <- as.data.frame(all_promoters) |>
  dplyr::select(seqnames, start, end, gene_id, strand) |>
  mutate(
    start = start - 1,
    score = 0,
    seqnames = as.character(seqnames),
    strand = as.character(strand),
    # Create unique ID based on coordinates
    name = paste0(seqnames, ":", start, "-", end)
  )

# Check for duplicates
cat("\nChecking for duplicate IDs...\n")
cat("Gained promoters - unique IDs:", length(unique(gained_df$name)), 
    "out of", nrow(gained_df), "\n")
cat("All promoters - unique IDs:", length(unique(all_df$name)), 
    "out of", nrow(all_df), "\n")

if (any(duplicated(gained_df$name))) {
  warning("Duplicates found in gained promoters! Removing...")
  gained_df <- gained_df[!duplicated(gained_df$name), ]
}

if (any(duplicated(all_df$name))) {
  warning("Duplicates found in all promoters! Removing...")
  all_df <- all_df[!duplicated(all_df$name), ]
}

# Write BED files (BED6 format: chr, start, end, name, score, strand)
write.table(gained_df[, c("seqnames", "start", "end", "name", "score", "strand")], 
            "data/processed/cutntag/homer_input/gained_anchor_promoters.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(all_df[, c("seqnames", "start", "end", "name", "score", "strand")], 
            "data/processed/cutntag/homer_input/all_promoters.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Save gene lists for reference -------------------------------------------
# Add gene symbols
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = gained_promoters$gene_id,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")

gained_promoters$symbol <- gene_symbols$SYMBOL[match(gained_promoters$gene_id, 
                                                     gene_symbols$ENTREZID)]

write.table(data.frame(gene_id = gained_promoters$gene_id,
                       symbol = gained_promoters$symbol),
            "data/processed/cutntag/homer_input/gained_anchor_genes.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# Verification ------------------------------------------------------------
cat("\n========================================\n")
cat("Verification\n")
cat("========================================\n")

# Read back and verify
gained_check <- read.table("data/processed/cutntag/homer_input/gained_anchor_promoters.bed",
                           header = FALSE, stringsAsFactors = FALSE)
all_check <- read.table("data/processed/cutntag/homer_input/all_promoters.bed",
                        header = FALSE, stringsAsFactors = FALSE)

cat("\nGained promoters BED:\n")
cat("  Rows:", nrow(gained_check), "\n")
cat("  Chromosomes start with 'chr':", grepl("^chr", gained_check[1,1]), "\n")
cat("  Mean region size:", round(mean(gained_check[,3] - gained_check[,2])), "bp\n")
cat("  First 3 rows:\n")
print(head(gained_check, 3))

cat("\nAll promoters BED:\n")
cat("  Rows:", nrow(all_check), "\n")
cat("  Chromosomes start with 'chr':", grepl("^chr", all_check[1,1]), "\n")
cat("  Mean region size:", round(mean(all_check[,3] - all_check[,2])), "bp\n")

# Summary -----------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("✓ Total promoters:", nrow(all_check), "\n")
cat("✓ Gained anchor promoters:", nrow(gained_check), "\n")
cat("✓ Percentage:", round(100 * nrow(gained_check) / nrow(all_check), 2), "%\n")
cat("\n")
cat("✓ Output files written to data/processed/cutntag/homer_input/\n")
cat("  - gained_anchor_promoters.bed\n")
cat("  - all_promoters.bed\n")
cat("  - gained_anchor_genes.txt\n")
cat("\n")
cat("✓ Chromosome format: UCSC style (chr1, chr2, etc.)\n")
cat("✓ All IDs are unique (no duplicates)\n")
cat("\n")
cat("Ready to run HOMER! Submit your SLURM job now.\n")
cat("========================================\n")
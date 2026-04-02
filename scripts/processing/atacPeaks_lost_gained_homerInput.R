# Identify ATAC peaks that overlap with gained and lost loop anchors

library(GenomicRanges)
library(plyranges)
library(dplyr)
library(InteractionSet)

# Load differential loops -------------------------------------------------
diff_loopCounts <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

# Extract GAINED loops ----------------------------------------------------
gained_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                        diff_loopCounts$log2FoldChange > 0)]

# Extract LOST loops ------------------------------------------------------
lost_loops <- diff_loopCounts[which(diff_loopCounts$padj < 0.1 & 
                                      diff_loopCounts$log2FoldChange < 0)]


# Get all loop anchors from gained loops ----------------------------------
gained_anchors <- c(anchors(gained_loops, type = "first"),
                    anchors(gained_loops, type = "second")) |>
  unique()

# CRITICAL: Ensure chromosome names have 'chr' prefix
seqlevelsStyle(gained_anchors) <- "UCSC"

# Get all loop anchors from lost loops ------------------------------------
lost_anchors <- c(anchors(lost_loops, type = "first"),
                  anchors(lost_loops, type = "second")) |>
  unique()

# CRITICAL: Ensure chromosome names have 'chr' prefix
seqlevelsStyle(lost_anchors) <- "UCSC"

# Load ATAC peaks ---------------------------------------------------------
atac_peaks <- plyranges::read_narrowpeaks(
  "data/processed/atac/output/peaks/YAPP_HEK_cont_0h_peaks.narrowPeak"
)

# CRITICAL: Ensure chromosome names match
seqlevelsStyle(atac_peaks) <- "UCSC"

# Keep only standard chromosomes
standard_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
atac_peaks <- atac_peaks[seqnames(atac_peaks) %in% standard_chrs]
gained_anchors <- gained_anchors[seqnames(gained_anchors) %in% standard_chrs]
lost_anchors <- lost_anchors[seqnames(lost_anchors) %in% standard_chrs]

# Find ATAC peaks overlapping gained anchors ------------------------------
gained_overlaps <- findOverlaps(atac_peaks, gained_anchors)
gained_atac_indices <- unique(queryHits(gained_overlaps))

gained_atac <- atac_peaks[gained_atac_indices]


# Find ATAC peaks overlapping lost anchors --------------------------------
lost_overlaps <- findOverlaps(atac_peaks, lost_anchors)
lost_atac_indices <- unique(queryHits(lost_overlaps))

lost_atac <- atac_peaks[lost_atac_indices]

# Export for HOMER --------------------------------------------------------
dir.create("data/processed/cutntag/homer_input", 
           showWarnings = FALSE, 
           recursive = TRUE)

# Convert to data frames for BED format
# BED format: chr, start (0-based), end, name, score, strand
gained_df <- as.data.frame(gained_atac) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)  # Convert to 0-based for BED

lost_df <- as.data.frame(lost_atac) |>
  dplyr::select(seqnames, start, end, name, score, strand) |>
  mutate(start = start - 1)

# Check for duplicates
cat("\nChecking for duplicate IDs...\n")
cat("Gained anchor ATAC peaks - unique IDs:", length(unique(gained_df$name)), 
    "out of", nrow(gained_df), "\n")
cat("Lost anchor ATAC peaks - unique IDs:", length(unique(lost_df$name)), 
    "out of", nrow(lost_df), "\n")

if (any(duplicated(gained_df$name))) {
  warning("Duplicates found in gained anchor ATAC peaks! Removing...")
  gained_df <- gained_df[!duplicated(gained_df$name), ]
}

if (any(duplicated(lost_df$name))) {
  warning("Duplicates found in lost anchor ATAC peaks! Removing...")
  lost_df <- lost_df[!duplicated(lost_df$name), ]
}

# Write BED files
write.table(gained_df, 
            "data/processed/atac/homer_input/gained_anchor_atac.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(lost_df, 
            "data/processed/atac/homer_input/lost_anchor_atac.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Verification ------------------------------------------------------------
cat("\n========================================\n")
cat("Verification\n")
cat("========================================\n")

# Read back and verify
gained_check <- read.table("data/processed/cutntag/homer_input/gained_anchor_atac.bed",
                           header = FALSE, stringsAsFactors = FALSE)
lost_check <- read.table("data/processed/cutntag/homer_input/lost_anchor_atac.bed",
                         header = FALSE, stringsAsFactors = FALSE)

cat("\nGained anchor ATAC peaks BED:\n")
cat("  Rows:", nrow(gained_check), "\n")
cat("  Chromosomes start with 'chr':", grepl("^chr", gained_check[1,1]), "\n")
cat("  Mean region size:", round(mean(gained_check[,3] - gained_check[,2])), "bp\n")
cat("  First 3 rows:\n")
print(head(gained_check, 3))

cat("\nLost anchor ATAC peaks BED:\n")
cat("  Rows:", nrow(lost_check), "\n")
cat("  Chromosomes start with 'chr':", grepl("^chr", lost_check[1,1]), "\n")
cat("  Mean region size:", round(mean(lost_check[,3] - lost_check[,2])), "bp\n")
cat("  First 3 rows:\n")
print(head(lost_check, 3))

# Summary -----------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("✓ Total ATAC peaks:", length(atac_peaks), "\n")
cat("✓ ATAC peaks at gained anchors:", nrow(gained_check), "\n")
cat("✓ ATAC peaks at lost anchors:", nrow(lost_check), "\n")
cat("✓ Gained percentage:", round(100 * nrow(gained_check) / length(atac_peaks), 2), "%\n")
cat("✓ Lost percentage:", round(100 * nrow(lost_check) / length(atac_peaks), 2), "%\n")
cat("\n")
cat("✓ Output files written to data/processed/cutntag/homer_input/\n")
cat("  - gained_anchor_atac.bed\n")
cat("  - lost_anchor_atac.bed\n")
cat("\n")
cat("✓ Chromosome format: UCSC style (chr1, chr2, etc.)\n")
cat("✓ All IDs are unique (no duplicates)\n")
cat("\n")
cat("Ready to run HOMER! Submit your SLURM job now.\n")
cat("========================================\n")
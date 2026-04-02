## =============================================================================
## EISA PART 1: BAM PROCESSING AND READ COUNTING
## =============================================================================
##
## Purpose: Count reads in exonic and intronic regions from BAM files
##
## Input:  BAM files from HISAT2 alignment
## Output: Exonic and intronic count matrices
##
## Runtime: ~30-40 minutes (the slow part!)
##
## =============================================================================

library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)

message("Starting EISA Part 1 (BAM processing): ", Sys.time())

## =============================================================================
## CONFIGURATION
## =============================================================================

## Base directories
BASE_DIR <- "/work/users/j/p/jpflores/projects/STRS"
BAM_DIR <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/align/stranded")
EISA_DIR <- file.path(BASE_DIR, "data/processed/rna/timecourse/output/EISA")

## Create output directory
dir.create(EISA_DIR, recursive = TRUE, showWarnings = FALSE)

## Annotation files
TXDB_FILE <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb"

## =============================================================================
## PART 1: DEFINE REGIONS
## =============================================================================

message("\n=== STEP 1: Loading annotations and defining exonic/intronic regions ===")

## Load TxDb
if (!file.exists(TXDB_FILE)) {
  stop("TxDb file not found: ", TXDB_FILE)
}

message("Loading pre-built TxDb from: ", TXDB_FILE)
txdb <- loadDb(TXDB_FILE)

## Get genes
genes_gr <- genes(txdb)
genes_total <- length(genes_gr)

## Keep only standard chromosomes (chr1-22, X, Y)
## Removes: contigs, patches, chrM
genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
genes_gr <- genes_gr[seqnames(genes_gr) != "chrM"]
seqlevelsStyle(genes_gr) <- "UCSC"

message(sprintf("Total genes in annotation: %d", genes_total))
message(sprintf("Genes on standard chromosomes (chr1-22, X, Y): %d (removed %d)", 
                length(genes_gr), genes_total - length(genes_gr)))

## Find non-overlapping genes
message("Finding non-overlapping genes...")
overlaps <- findOverlaps(genes_gr, genes_gr, ignore.strand = FALSE)
overlap_counts <- table(queryHits(overlaps))
nonoverlap_indices <- as.integer(names(overlap_counts[overlap_counts == 1]))
nonoverlap_genes <- genes_gr[nonoverlap_indices]

message(sprintf("Non-overlapping genes: %d (%.1f%%)", 
                length(nonoverlap_genes),
                length(nonoverlap_genes) / length(genes_gr) * 100))

## Get exons by gene
message("Extracting exons...")
exons_by_gene <- exonsBy(txdb, by = "gene")
exons_by_gene <- exons_by_gene[names(nonoverlap_genes)]

## Create SAF helper function
create_saf <- function(gr_list, extend_exons = 10) {
  if (extend_exons > 0) {
    gr_list <- endoapply(gr_list, function(x) {
      resize(x, width = width(x) + 2 * extend_exons, fix = "center")
    })
  }
  
  saf <- data.frame(
    GeneID = rep(names(gr_list), sapply(gr_list, length)),
    Chr = as.character(seqnames(unlist(gr_list))),
    Start = start(unlist(gr_list)),
    End = end(unlist(gr_list)),
    Strand = as.character(strand(unlist(gr_list)))
  )
  return(saf)
}

## Create exon SAF (extend by 10bp to avoid junction reads)
message("Creating exon annotation file...")
exon_saf <- create_saf(exons_by_gene, extend_exons = 10)
message(sprintf("Exonic regions defined: %d features", nrow(exon_saf)))

## Create intron SAF
message("Creating intron annotation file...")
create_intron_saf <- function(genes, exons_by_gene) {
  intron_list <- list()
  
  for (gene_id in names(exons_by_gene)) {
    gene_body <- genes[gene_id]
    gene_exons <- exons_by_gene[[gene_id]]
    
    ## Extend exons by 10bp
    extended_exons <- resize(gene_exons, 
                             width = width(gene_exons) + 20, 
                             fix = "center")
    
    ## Get introns
    introns <- setdiff(gene_body, extended_exons)
    
    if (length(introns) > 0) {
      intron_list[[gene_id]] <- introns
    }
  }
  
  return(create_saf(GRangesList(intron_list), extend_exons = 0))
}

intron_saf <- create_intron_saf(nonoverlap_genes, exons_by_gene)
message(sprintf("Intronic regions defined: %d features", nrow(intron_saf)))

## Save SAF files
write.table(exon_saf, 
            file.path(EISA_DIR, "exon_regions.saf"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(intron_saf, 
            file.path(EISA_DIR, "intron_regions.saf"),
            quote = FALSE, sep = "\t", row.names = FALSE)

message("✓ Region definitions saved")

## =============================================================================
## PART 2: COUNT READS
## =============================================================================

message("\n=== STEP 2: Counting reads in exonic and intronic regions ===")

## Check for existing counts
exon_counts_file <- file.path(EISA_DIR, "exon_counts.rds")
intron_counts_file <- file.path(EISA_DIR, "intron_counts.rds")

if (file.exists(exon_counts_file) && file.exists(intron_counts_file)) {
  message("✓ Count files already exist! Skipping featureCounts.")
  message("  Exon counts: ", exon_counts_file)
  message("  Intron counts: ", intron_counts_file)
  message("\nTo recount, delete these files and rerun this script.")
  quit(save = "no", status = 0)
}

## Get BAM files
bam_files <- list.files(BAM_DIR, 
                        pattern = "\\.bam$", 
                        full.names = TRUE,
                        recursive = FALSE)

if (length(bam_files) == 0) {
  stop("No BAM files found in: ", BAM_DIR)
}

message(sprintf("Found %d BAM files", length(bam_files)))

## Check if BAM files are indexed
message("Checking BAM file indices...")
bai_files <- paste0(bam_files, ".bai")
missing_indices <- !file.exists(bai_files)

if (any(missing_indices)) {
  message(sprintf("Indexing %d BAM files...", sum(missing_indices)))
  for (bam in bam_files[missing_indices]) {
    message("  Indexing: ", basename(bam))
    system(paste("samtools index", bam))
  }
}

## featureCounts parameters
IS_PAIRED_END <- TRUE
STRAND_SPECIFIC <- 2  # 2 = reverse stranded (typical for Illumina TruSeq)
N_THREADS <- 16

message("\nfeatureCounts settings:")
message("  Paired-end: ", IS_PAIRED_END)
message("  Strand-specific: ", STRAND_SPECIFIC, " (reverse)")
message("  Threads: ", N_THREADS)

## Count exonic reads
message("\n--- Counting EXONIC reads ---")
message("This will take ~15-20 minutes...")

exon_fc <- featureCounts(
  files = bam_files,
  annot.ext = exon_saf,
  useMetaFeatures = TRUE,
  allowMultiOverlap = TRUE,
  isPairedEnd = IS_PAIRED_END,
  strandSpecific = STRAND_SPECIFIC,
  countReadPairs = TRUE,
  requireBothEndsMapped = TRUE,
  nthreads = N_THREADS,
  verbose = TRUE
)

## Count intronic reads
message("\n--- Counting INTRONIC reads ---")
message("This will take ~15-20 minutes...")

intron_fc <- featureCounts(
  files = bam_files,
  annot.ext = intron_saf,
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  isPairedEnd = IS_PAIRED_END,
  strandSpecific = STRAND_SPECIFIC,
  countReadPairs = TRUE,
  requireBothEndsMapped = TRUE,
  nthreads = N_THREADS,
  verbose = TRUE
)

## Extract count matrices
exon_counts <- exon_fc$counts
intron_counts <- intron_fc$counts

## Fix column names (remove path, keep sample ID)
colnames(exon_counts) <- basename(colnames(exon_counts)) |> str_remove("\\.bam$")
colnames(intron_counts) <- basename(colnames(intron_counts)) |> str_remove("\\.bam$")

## Save counts
message("\nSaving count matrices...")
saveRDS(exon_counts, exon_counts_file)
saveRDS(intron_counts, intron_counts_file)

## Save detailed statistics
write.csv(exon_fc$stat, 
          file.path(EISA_DIR, "exon_counts_stats.csv"),
          row.names = FALSE)
write.csv(intron_fc$stat,
          file.path(EISA_DIR, "intron_counts_stats.csv"),
          row.names = FALSE)

message("✓ Count matrices saved")

## =============================================================================
## QUALITY CHECKS
## =============================================================================

message("\n=== Quality Check: Intronic Coverage ===")

intronic_pct <- colSums(intron_counts) / 
  (colSums(intron_counts) + colSums(exon_counts)) * 100

message(sprintf("Mean intronic percentage: %.1f%%", mean(intronic_pct)))
message(sprintf("Range: %.1f%% - %.1f%%", min(intronic_pct), max(intronic_pct)))
message("\nExpected ranges:")
message("  PolyA RNA: 5-15%")
message("  Total RNA: 15-30%")

if (mean(intronic_pct) < 3) {
  warning("Very low intronic coverage (<3%). EISA results may be unreliable.")
} else if (mean(intronic_pct) < 5) {
  warning("Low intronic coverage (3-5%). Consider using total RNA for better EISA signal.")
} else {
  message("✓ Intronic coverage is sufficient for EISA")
}

## =============================================================================
## COMPLETION
## =============================================================================

message("\n")
message("================================================================================")
message("                    PART 1 COMPLETE: BAM PROCESSING")
message("================================================================================")
message("\n")
message("Output files:")
message(sprintf("  1. Exon counts: %s", exon_counts_file))
message(sprintf("  2. Intron counts: %s", intron_counts_file))
message(sprintf("  3. Exon SAF: %s", file.path(EISA_DIR, "exon_regions.saf")))
message(sprintf("  4. Intron SAF: %s", file.path(EISA_DIR, "intron_regions.saf")))
message("\n")
message("Summary:")
message(sprintf("  - %d samples processed", ncol(exon_counts)))
message(sprintf("  - %d genes with exonic counts", nrow(exon_counts)))
message(sprintf("  - %d genes with intronic counts", nrow(intron_counts)))
message(sprintf("  - Mean intronic percentage: %.1f%%", mean(intronic_pct)))
message("\n")
message("Next step: Run Part 2 (DESeq2 analysis)")
message("  Rscript scripts/EISA/EISA_part2_DESeq2.R")
message("\n")
message("Completed: ", Sys.time())
message("================================================================================")
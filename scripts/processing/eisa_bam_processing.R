# ##############################################################################
# filename:    eisa_bam_processing.R
# author:      JP Flores
# project:     STRS
# date:        2026-04-02
# description: EISA Part 1 — count reads in exonic and intronic regions from
#              HISAT2-aligned BAM files using featureCounts with 16 threads;
#              runtime ~30-40 minutes
# ##############################################################################

# Libraries ----
library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)

# Parameters ----
base_dir         <- "/work/users/j/p/jpflores/projects/STRS"
bam_dir          <- file.path(base_dir, "data/processed/rna/timecourse/output/align/stranded")
eisa_dir         <- file.path(base_dir, "data/processed/rna/timecourse/output/EISA")
txdb_file        <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb"
is_paired_end    <- TRUE
strand_specific  <- 2  # 2 = reverse stranded (typical for Illumina TruSeq)
n_threads        <- 16

exon_counts_file   <- file.path(eisa_dir, "exon_counts.rds")
intron_counts_file <- file.path(eisa_dir, "intron_counts.rds")

# Data import ----
message("Starting EISA Part 1 (BAM processing): ", Sys.time())
dir.create(eisa_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== STEP 1: Loading annotations and defining exonic/intronic regions ===")

if (!file.exists(txdb_file)) {
  stop("TxDb file not found: ", txdb_file)
}

message("Loading pre-built TxDb from: ", txdb_file)
txdb <- loadDb(txdb_file)

## Get genes on standard chromosomes
genes_gr    <- genes(txdb)
genes_total <- length(genes_gr)

genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
genes_gr <- genes_gr[seqnames(genes_gr) != "chrM"]
seqlevelsStyle(genes_gr) <- "UCSC"

message(sprintf("Total genes in annotation: %d", genes_total))
message(sprintf("Genes on standard chromosomes (chr1-22, X, Y): %d (removed %d)",
                length(genes_gr), genes_total - length(genes_gr)))

## Find non-overlapping genes
message("Finding non-overlapping genes...")
overlaps       <- findOverlaps(genes_gr, genes_gr, ignore.strand = FALSE)
overlap_counts <- table(queryHits(overlaps))
nonoverlap_indices <- as.integer(names(overlap_counts[overlap_counts == 1]))
nonoverlap_genes   <- genes_gr[nonoverlap_indices]

message(sprintf("Non-overlapping genes: %d (%.1f%%)",
                length(nonoverlap_genes),
                length(nonoverlap_genes) / length(genes_gr) * 100))

## Get exons by gene
message("Extracting exons...")
exons_by_gene <- exonsBy(txdb, by = "gene")
exons_by_gene <- exons_by_gene[names(nonoverlap_genes)]

# Analysis ----

## Helper: build SAF data frame from a GRangesList
create_saf <- function(gr_list, extend_exons = 10) {
  if (extend_exons > 0) {
    gr_list <- endoapply(gr_list, function(x) {
      resize(x, width = width(x) + 2 * extend_exons, fix = "center")
    })
  }

  saf <- data.frame(
    GeneID = rep(names(gr_list), sapply(gr_list, length)),
    Chr    = as.character(seqnames(unlist(gr_list))),
    Start  = start(unlist(gr_list)),
    End    = end(unlist(gr_list)),
    Strand = as.character(strand(unlist(gr_list)))
  )
  return(saf)
}

## Helper: build intron SAF (gene body minus extended exons)
create_intron_saf <- function(genes, exons_by_gene) {
  intron_list <- list()

  for (gene_id in names(exons_by_gene)) {
    gene_body    <- genes[gene_id]
    gene_exons   <- exons_by_gene[[gene_id]]
    extended_exons <- resize(gene_exons,
                             width = width(gene_exons) + 20,
                             fix = "center")
    introns <- setdiff(gene_body, extended_exons)

    if (length(introns) > 0) {
      intron_list[[gene_id]] <- introns
    }
  }

  return(create_saf(GRangesList(intron_list), extend_exons = 0))
}

## Create exon SAF (extend by 10 bp to avoid junction reads)
message("Creating exon annotation file...")
exon_saf <- create_saf(exons_by_gene, extend_exons = 10)
message(sprintf("Exonic regions defined: %d features", nrow(exon_saf)))

## Create intron SAF
message("Creating intron annotation file...")
intron_saf <- create_intron_saf(nonoverlap_genes, exons_by_gene)
message(sprintf("Intronic regions defined: %d features", nrow(intron_saf)))

## Save SAF files
write.table(exon_saf,   file.path(eisa_dir, "exon_regions.saf"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(intron_saf, file.path(eisa_dir, "intron_regions.saf"),
            quote = FALSE, sep = "\t", row.names = FALSE)

message("Region definitions saved")

message("\n=== STEP 2: Counting reads in exonic and intronic regions ===")

## Skip if counts already exist
if (file.exists(exon_counts_file) && file.exists(intron_counts_file)) {
  message("Count files already exist! Skipping featureCounts.")
  message("  Exon counts: ", exon_counts_file)
  message("  Intron counts: ", intron_counts_file)
  message("\nTo recount, delete these files and rerun this script.")
  quit(save = "no", status = 0)
}

## Get BAM files
bam_files <- list.files(bam_dir,
                        pattern = "\\.bam$",
                        full.names = TRUE,
                        recursive = FALSE)

if (length(bam_files) == 0) {
  stop("No BAM files found in: ", bam_dir)
}

message(sprintf("Found %d BAM files", length(bam_files)))

## Index any unindexed BAM files
bai_files      <- paste0(bam_files, ".bai")
missing_indices <- !file.exists(bai_files)

if (any(missing_indices)) {
  message(sprintf("Indexing %d BAM files...", sum(missing_indices)))
  for (bam in bam_files[missing_indices]) {
    message("  Indexing: ", basename(bam))
    system(paste("samtools index", bam))
  }
}

message("\nfeatureCounts settings:")
message("  Paired-end: ", is_paired_end)
message("  Strand-specific: ", strand_specific, " (reverse)")
message("  Threads: ", n_threads)

## Count exonic reads
message("\n--- Counting EXONIC reads ---")
exon_fc <- featureCounts(
  files                 = bam_files,
  annot.ext             = exon_saf,
  useMetaFeatures       = TRUE,
  allowMultiOverlap     = TRUE,
  isPairedEnd           = is_paired_end,
  strandSpecific        = strand_specific,
  countReadPairs        = TRUE,
  requireBothEndsMapped = TRUE,
  nthreads              = n_threads,
  verbose               = TRUE
)

## Count intronic reads
message("\n--- Counting INTRONIC reads ---")
intron_fc <- featureCounts(
  files                 = bam_files,
  annot.ext             = intron_saf,
  useMetaFeatures       = TRUE,
  allowMultiOverlap     = FALSE,
  isPairedEnd           = is_paired_end,
  strandSpecific        = strand_specific,
  countReadPairs        = TRUE,
  requireBothEndsMapped = TRUE,
  nthreads              = n_threads,
  verbose               = TRUE
)

## Extract count matrices
exon_counts   <- exon_fc$counts
intron_counts <- intron_fc$counts

## Fix column names
colnames(exon_counts)   <- basename(colnames(exon_counts))   |> str_remove("\\.bam$")
colnames(intron_counts) <- basename(colnames(intron_counts)) |> str_remove("\\.bam$")

# Save outputs ----
message("\nSaving count matrices...")
saveRDS(exon_counts,   exon_counts_file)
saveRDS(intron_counts, intron_counts_file)

write.csv(exon_fc$stat,
          file.path(eisa_dir, "exon_counts_stats.csv"),
          row.names = FALSE)
write.csv(intron_fc$stat,
          file.path(eisa_dir, "intron_counts_stats.csv"),
          row.names = FALSE)

## Quality check: intronic coverage
message("\n=== Quality Check: Intronic Coverage ===")
intronic_pct <- colSums(intron_counts) /
  (colSums(intron_counts) + colSums(exon_counts)) * 100

message(sprintf("Mean intronic percentage: %.1f%%", mean(intronic_pct)))
message(sprintf("Range: %.1f%% - %.1f%%", min(intronic_pct), max(intronic_pct)))

if (mean(intronic_pct) < 3) {
  warning("Very low intronic coverage (<3%). EISA results may be unreliable.")
} else if (mean(intronic_pct) < 5) {
  warning("Low intronic coverage (3-5%). Consider using total RNA for better EISA signal.")
} else {
  message("Intronic coverage is sufficient for EISA")
}

message("\nOutput files:")
message(sprintf("  1. Exon counts:   %s", exon_counts_file))
message(sprintf("  2. Intron counts: %s", intron_counts_file))
message(sprintf("  3. Exon SAF:      %s", file.path(eisa_dir, "exon_regions.saf")))
message(sprintf("  4. Intron SAF:    %s", file.path(eisa_dir, "intron_regions.saf")))
message("\nCompleted: ", Sys.time())

sessionInfo()

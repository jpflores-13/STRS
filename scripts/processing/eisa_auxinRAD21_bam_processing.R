# ##############################################################################
# filename:    eisa_auxinRAD21_bam_processing.R
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: EISA Part 1 for AuxinRAD21 data — counts reads in exonic and
#              intronic regions from per-replicate HISAT2-aligned BAM files
#              (HCT116 mAID2-RAD21: untreated rep1/rep2, sorbitol 6h rep1/rep2,
#              sorbitol+auxin 6h rep1/rep2) using featureCounts with 16 threads.
#              Uses per-replicate BAMs in align/ (n=6) so that eisaR/edgeR can
#              estimate dispersion properly with 2 replicates per condition.
#              Outputs go to EISA/per_replicate/ to avoid overwriting the
#              existing merged-BAM EISA results.
#              NOTE: designed to be submitted as a SLURM job via
#              eisa_auxinRAD21_bam_processing.sh — do not run interactively
#              as featureCounts is memory and thread intensive.
# input:       align/ per-replicate BAM files from AuxinRAD21 bagPipes run
# output:      data/processed/rna/AuxinRAD21/output/EISA/per_replicate/exon_counts.rds
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/intron_counts.rds
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/exon_regions.saf
#              data/processed/rna/AuxinRAD21/output/EISA/per_replicate/intron_regions.saf
# ##############################################################################


# Libraries -------------------------------------------------------------------

library(GenomicFeatures)
library(GenomicRanges)
library(Rsubread)
library(stringr)
library(GenomeInfoDb)


# Parameters ------------------------------------------------------------------

project_dir <- "/work/users/j/p/jpflores/projects/STRS"

## Per-replicate BAM directory — 6 BAMs (2 reps x 3 conditions)
bam_dir <- file.path(project_dir,
                     "data/processed/rna/AuxinRAD21/output/align")

## Separate output directory from merged-BAM EISA results
eisa_dir <- file.path(project_dir,
                      "data/processed/rna/AuxinRAD21/output/EISA/per_replicate")

## GENCODE v49 TxDb
txdb_file <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.TxDb"

## featureCounts settings
is_paired_end   <- TRUE
strand_specific <- 2      # reverse stranded (Illumina TruSeq)
n_threads       <- 16

## Output paths
exon_counts_file   <- file.path(eisa_dir, "exon_counts.rds")
intron_counts_file <- file.path(eisa_dir, "intron_counts.rds")


# Setup -----------------------------------------------------------------------

message("Starting AuxinRAD21 EISA BAM processing (per-replicate): ", Sys.time())
dir.create(eisa_dir, recursive = TRUE, showWarnings = FALSE)


# Annotation ------------------------------------------------------------------

message("Loading TxDb: ", txdb_file)
if (!file.exists(txdb_file)) stop("TxDb not found: ", txdb_file)
txdb <- loadDb(txdb_file)

genes_gr    <- genes(txdb)
genes_total <- length(genes_gr)
genes_gr    <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
genes_gr    <- genes_gr[seqnames(genes_gr) != "chrM"]
seqlevelsStyle(genes_gr) <- "UCSC"

message(sprintf("Total genes: %d | Standard chromosomes: %d (removed %d)",
                genes_total, length(genes_gr),
                genes_total - length(genes_gr)))

## Non-overlapping genes only
overlaps         <- findOverlaps(genes_gr, genes_gr, ignore.strand = FALSE)
overlap_counts   <- table(queryHits(overlaps))
nonoverlap_idx   <- as.integer(names(overlap_counts[overlap_counts == 1L]))
nonoverlap_genes <- genes_gr[nonoverlap_idx]

message(sprintf("Non-overlapping genes: %d (%.1f%%)",
                length(nonoverlap_genes),
                length(nonoverlap_genes) / length(genes_gr) * 100))

exons_by_gene <- exonsBy(txdb, by = "gene")
exons_by_gene <- exons_by_gene[names(nonoverlap_genes)]


# SAF helpers -----------------------------------------------------------------

create_saf <- function(gr_list, extend_exons = 10L) {
  if (extend_exons > 0L) {
    gr_list <- endoapply(gr_list, function(x) {
      resize(x, width = width(x) + 2L * extend_exons, fix = "center")
    })
  }
  data.frame(
    GeneID = rep(names(gr_list), sapply(gr_list, length)),
    Chr    = as.character(seqnames(unlist(gr_list))),
    Start  = start(unlist(gr_list)),
    End    = end(unlist(gr_list)),
    Strand = as.character(strand(unlist(gr_list)))
  )
}

create_intron_saf <- function(genes, exons_by_gene) {
  intron_list <- list()
  for (gene_id in names(exons_by_gene)) {
    gene_body      <- genes[gene_id]
    gene_exons     <- exons_by_gene[[gene_id]]
    extended_exons <- resize(gene_exons, width = width(gene_exons) + 20L,
                             fix = "center")
    introns <- setdiff(gene_body, extended_exons)
    if (length(introns) > 0L) intron_list[[gene_id]] <- introns
  }
  create_saf(GRangesList(intron_list), extend_exons = 0L)
}

message("Building exon SAF...")
exon_saf <- create_saf(exons_by_gene, extend_exons = 10L)
message(sprintf("  Exonic features: %d", nrow(exon_saf)))

message("Building intron SAF...")
intron_saf <- create_intron_saf(nonoverlap_genes, exons_by_gene)
message(sprintf("  Intronic features: %d", nrow(intron_saf)))

write.table(exon_saf,   file.path(eisa_dir, "exon_regions.saf"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(intron_saf, file.path(eisa_dir, "intron_regions.saf"),
            quote = FALSE, sep = "\t", row.names = FALSE)


# BAM files -------------------------------------------------------------------

bam_files <- list.files(bam_dir,
                        pattern = "_sorted\\.bam$",
                        full.names = TRUE)
bam_files <- bam_files[!grepl("\\.bai$", bam_files)]

if (length(bam_files) == 0L) stop("No BAM files found in: ", bam_dir)
message(sprintf("Found %d BAM files:", length(bam_files)))
message(paste0("  ", basename(bam_files), collapse = "\n"))

## Index any un-indexed BAMs
bai_files       <- paste0(bam_files, ".bai")
missing_indices <- !file.exists(bai_files)
if (any(missing_indices)) {
  message("Indexing ", sum(missing_indices), " BAM files...")
  for (bam in bam_files[missing_indices]) {
    message("  samtools index ", basename(bam))
    system(paste("samtools index", bam))
  }
}


# featureCounts ---------------------------------------------------------------

## Wrap in if/else rather than quit() so the script never kills an interactive
## session — safe to run both interactively and as a SLURM batch job
if (file.exists(exon_counts_file) && file.exists(intron_counts_file)) {
  
  message("Count files already exist in per_replicate/ — skipping featureCounts.")
  message("  Exon:   ", exon_counts_file)
  message("  Intron: ", intron_counts_file)
  message("To recount, delete those files and rerun.")
  
} else {
  
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
  
  ## Clean column names: strip path, .bam, _sorted, and experiment prefix
  ## e.g. YAPP_HCT116_mAID2-RAD21_sorbitol_6h_1_1_sorted -> sorbitol_6h_1_1
  clean_colnames <- function(mat) {
    cn <- basename(colnames(mat)) |>
      str_remove("\\.bam$") |>
      str_remove("_sorted$") |>
      str_remove("^YAPP_HCT116_mAID2-RAD21_")
    colnames(mat) <- cn
    mat
  }
  
  exon_counts   <- clean_colnames(exon_fc$counts)
  intron_counts <- clean_colnames(intron_fc$counts)
  
  message("Final column names: ", paste(colnames(exon_counts), collapse = ", "))
  
  ## Save
  saveRDS(exon_counts,   exon_counts_file)
  saveRDS(intron_counts, intron_counts_file)
  
  write.csv(exon_fc$stat,
            file.path(eisa_dir, "exon_counts_stats.csv"),   row.names = FALSE)
  write.csv(intron_fc$stat,
            file.path(eisa_dir, "intron_counts_stats.csv"), row.names = FALSE)
  
  ## Intronic coverage QC
  intronic_pct <- colSums(intron_counts) /
    (colSums(intron_counts) + colSums(exon_counts)) * 100
  
  message("\nIntronic coverage QC:")
  message(sprintf("  Mean: %.1f%%  Range: %.1f%% - %.1f%%",
                  mean(intronic_pct), min(intronic_pct), max(intronic_pct)))
  
  if (mean(intronic_pct) < 3) {
    warning("Very low intronic coverage (<3%). EISA results may be unreliable.")
  } else if (mean(intronic_pct) < 5) {
    warning("Low intronic coverage (3-5%). Consider total RNA for better EISA signal.")
  } else {
    message("  Intronic coverage sufficient for EISA.")
  }
  
  message("\nOutputs:")
  message("  ", exon_counts_file)
  message("  ", intron_counts_file)
  message("Completed: ", Sys.time())
  
}


# Session info ----------------------------------------------------------------

sessionInfo()
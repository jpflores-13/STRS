## Validate gene compatibility with plotgardener
## Ensures genes work with TxDb.Hsapiens.UCSC.hg38.knownGene

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)

# Function to validate and clean genes for plotgardener compatibility --------

validate_genes_for_plotgardener <- function(gene_list, verbose = TRUE) {
  
  if (verbose) cat("=== PLOTGARDENER GENE COMPATIBILITY CHECK ===\n")
  
  # Load the TxDb that plotgardener uses
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  if (verbose) cat("Input genes:", length(gene_list), "\n")
  
  # Method 1: Check directly against TxDb gene symbols
  txdb_genes <- genes(txdb)
  
  # Get gene symbols from TxDb using the same org.Hs.eg.db mapping
  txdb_entrez <- mcols(txdb_genes)$gene_id
  txdb_symbols <- mapIds(org.Hs.eg.db, 
                         keys = txdb_entrez,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  # Clean up - remove NAs and create lookup table
  valid_txdb_symbols <- txdb_symbols[!is.na(txdb_symbols)]
  
  if (verbose) {
    cat("TxDb contains", length(valid_txdb_symbols), "gene symbols\n")
  }
  
  # Check which of our genes are in TxDb
  genes_in_txdb <- gene_list[gene_list %in% valid_txdb_symbols]
  genes_not_in_txdb <- gene_list[!gene_list %in% valid_txdb_symbols]
  
  if (verbose) {
    cat("✓ Compatible genes:", length(genes_in_txdb), "\n")
    cat("✗ Incompatible genes:", length(genes_not_in_txdb), "\n")
  }
  
  # Method 2: Try to rescue incompatible genes with alternative symbols
  rescued_genes <- c()
  
  if (length(genes_not_in_txdb) > 0 && verbose) {
    cat("\nAttempting to rescue incompatible genes...\n")
  }
  
  for (gene in genes_not_in_txdb) {
    # Try to get Entrez ID first
    entrez_id <- mapIds(org.Hs.eg.db, 
                        keys = gene,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
    
    if (!is.na(entrez_id)) {
      # Check if this Entrez ID is in TxDb
      if (entrez_id %in% txdb_entrez) {
        # Get the "official" symbol used in TxDb
        official_symbol <- mapIds(org.Hs.eg.db,
                                  keys = entrez_id,
                                  column = "SYMBOL", 
                                  keytype = "ENTREZID",
                                  multiVals = "first")
        
        if (!is.na(official_symbol) && official_symbol %in% valid_txdb_symbols) {
          rescued_genes <- c(rescued_genes, official_symbol)
          if (verbose) {
            cat("  Rescued:", gene, "→", official_symbol, "\n")
          }
        }
      }
    }
  }
  
  # Final compatible gene list
  final_genes <- unique(c(genes_in_txdb, rescued_genes))
  
  if (verbose) {
    cat("\n=== FINAL RESULTS ===\n")
    cat("Original genes:", length(gene_list), "\n")
    cat("Compatible genes:", length(genes_in_txdb), "\n") 
    cat("Rescued genes:", length(rescued_genes), "\n")
    cat("Final compatible genes:", length(final_genes), "\n")
    cat("Success rate:", round(length(final_genes)/length(gene_list)*100, 1), "%\n")
    
    if (length(genes_not_in_txdb) - length(rescued_genes) > 0) {
      still_incompatible <- genes_not_in_txdb[!genes_not_in_txdb %in% rescued_genes]
      cat("\nStill incompatible genes (", length(still_incompatible), "):\n")
      cat(paste(head(still_incompatible, 10), collapse = ", "))
      if (length(still_incompatible) > 10) cat(", ...")
      cat("\n")
    }
  }
  
  return(list(
    compatible_genes = final_genes,
    original_count = length(gene_list),
    compatible_count = length(final_genes),
    incompatible_genes = genes_not_in_txdb[!genes_not_in_txdb %in% rescued_genes],
    rescued_genes = rescued_genes
  ))
}

# Function to test genes with actual plotgardener calls -------------------

test_plotgardener_compatibility <- function(gene_list, n_test = 5) {
  
  cat("=== TESTING ACTUAL PLOTGARDENER COMPATIBILITY ===\n")
  
  # Test a subset of genes
  test_genes <- head(gene_list, n_test)
  
  successful_genes <- c()
  failed_genes <- c()
  
  for (gene in test_genes) {
    tryCatch({
      # Try to create pgParams - this will fail if gene not found
      p <- pgParams(assembly = "hg38",
                    resolution = 10e3,
                    gene = gene,
                    geneBuffer = 50e3)  # Small buffer for testing
      
      # If we get here, the gene was found
      successful_genes <- c(successful_genes, gene)
      cat("✓", gene, "- OK\n")
      
    }, error = function(e) {
      failed_genes <- c(failed_genes, gene)
      cat("✗", gene, "- FAILED:", e$message, "\n")
    })
  }
  
  cat("\nTest results:\n")
  cat("Successful:", length(successful_genes), "/", length(test_genes), "\n")
  
  return(list(
    successful = successful_genes,
    failed = failed_genes
  ))
}

# Function to get comprehensive gene information for debugging -------------

get_gene_info <- function(gene_symbols) {
  
  cat("=== COMPREHENSIVE GENE INFORMATION ===\n")
  
  gene_info <- map_dfr(gene_symbols, function(gene) {
    
    # Get Entrez ID
    entrez_id <- mapIds(org.Hs.eg.db,
                        keys = gene,
                        column = "ENTREZID", 
                        keytype = "SYMBOL",
                        multiVals = "first")
    
    # Get alternative symbols
    aliases <- mapIds(org.Hs.eg.db,
                      keys = gene,
                      column = "ALIAS",
                      keytype = "SYMBOL", 
                      multiVals = "list")[[1]]
    
    # Check if in TxDb
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    txdb_genes <- genes(txdb)
    in_txdb <- !is.na(entrez_id) && entrez_id %in% mcols(txdb_genes)$gene_id
    
    # Get chromosome info if available
    chr_info <- NA
    if (!is.na(entrez_id) && in_txdb) {
      gene_ranges <- genes(txdb, filter = list(gene_id = entrez_id))
      if (length(gene_ranges) > 0) {
        chr_info <- as.character(seqnames(gene_ranges))[1]
      }
    }
    
    tibble(
      symbol = gene,
      entrez_id = ifelse(is.na(entrez_id), "Not found", entrez_id),
      in_txdb = in_txdb,
      chromosome = ifelse(is.na(chr_info), "Unknown", chr_info),
      alias_count = ifelse(is.null(aliases), 0, length(aliases)),
      first_alias = ifelse(is.null(aliases) || length(aliases) == 0, 
                           "None", aliases[1])
    )
  })
  
  return(gene_info)
}

# Alternative gene selection approach using TxDb directly -----------------

get_genes_from_txdb_directly <- function(go_terms) {
  
  cat("=== ALTERNATIVE: GET GENES DIRECTLY FROM TxDb ===\n")
  
  # Load TxDb
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Get all genes from TxDb with their Entrez IDs
  all_txdb_genes <- genes(txdb)
  all_entrez_ids <- mcols(all_txdb_genes)$gene_id
  
  # For each GO term, get genes and filter to those in TxDb
  txdb_compatible_genes <- c()
  
  for (go_id in go_terms$go_id) {
    # Get genes for this GO term
    go_gene_data <- select(org.Hs.eg.db,
                           keys = go_id,
                           columns = c("ENTREZID", "SYMBOL"),
                           keytype = "GO")
    
    if (nrow(go_gene_data) > 0) {
      # Filter to genes present in TxDb
      txdb_genes_this_go <- go_gene_data |> 
        filter(ENTREZID %in% all_entrez_ids) |> 
        filter(!is.na(SYMBOL)) |> 
        pull(SYMBOL) |> 
        unique()
      
      txdb_compatible_genes <- c(txdb_compatible_genes, txdb_genes_this_go)
    }
  }
  
  final_genes <- unique(txdb_compatible_genes)
  
  cat("Genes directly from TxDb:", length(final_genes), "\n")
  
  return(final_genes)
}

# Usage example for your workflow -----------------------------------------

validate_and_clean_gene_list <- function(go_stress_genes) {
  
  cat("=== CLEANING GENE LIST FOR PLOTGARDENER ===\n\n")
  
  # Step 1: Basic validation
  validation_results <- validate_genes_for_plotgardener(go_stress_genes$SYMBOL)
  
  # Step 2: Test a few genes with actual plotgardener
  if (length(validation_results$compatible_genes) > 0) {
    cat("\n")
    test_results <- test_plotgardener_compatibility(validation_results$compatible_genes, n_test = 3)
  }
  
  # Step 3: Return cleaned gene list
  clean_genes <- validation_results$compatible_genes
  
  cat("\n=== RECOMMENDATION ===\n")
  cat("Use", length(clean_genes), "validated genes for plotgardener\n")
  
  # Save diagnostic information
  write_csv(tibble(gene_symbol = clean_genes), 
            "data/processed/plotgardener_compatible_genes.csv")
  
  if (length(validation_results$incompatible_genes) > 0) {
    write_csv(tibble(gene_symbol = validation_results$incompatible_genes),
              "data/processed/incompatible_genes.csv")
  }
  
  return(clean_genes)
}
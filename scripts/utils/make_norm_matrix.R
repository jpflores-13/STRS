# make_norm_matrix.R ----------------------------------------------------- ## Build a z-score normalized count matrix from a DESeq2 object

# Libraries ----
library(DESeq2)
library(matrixStats)

#' Build a z-score normalized count matrix for clustering
#'
#' VST-normalizes a DESeq2 object, subsets to significant genes defined by the
#' caller-supplied \code{sigLRT_genes} vector, and z-score normalizes each row.
#'
#' @param dds DESeqDataSet. A fitted DESeq2 object
#' @return matrix. Z-score normalized count matrix (genes x samples), subsetted
#'   to rows named in \code{sigLRT_genes} (must exist in the calling environment)
make_norm_matrix <- function(dds) {
  # normalize dds
  dds_norm <- vst(dds)
  # pull out only the genes in sigLRT_genes
  count_matrix <- dds_norm[rownames(dds_norm) %in% sigLRT_genes,]
  # make the data into a matrix that we can easily read
  count_matrix <- assay(count_matrix)
  # z-score normalize rows
  count_matrix_norm <- (count_matrix - rowMeans(count_matrix)) / rowSds(count_matrix + .5)
  return(count_matrix_norm)
}

# Footer ----
message("make_norm_matrix.R loaded. Available functions:")
message("  make_norm_matrix(dds)")
message("  Note: requires 'sigLRT_genes' character vector in the calling environment")
message("Example: mat <- make_norm_matrix(dds)")

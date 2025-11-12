# Make normalized count matrix for clustering and plotting
make_norm_matrix <- function(dds){
  #normalize dds
  dds_norm <- vst(dds)
  #pull out only the genes in sigLRT_genes
  count_matrix <- dds_norm[rownames(dds_norm) %in% sigLRT_genes,]
  #make the data into a matrix that we can easily read
  count_matrix <- assay(count_matrix)
  #normalize to 0 and declare to global environment
  count_matrix_norm <- (count_matrix-rowMeans(count_matrix))/rowSds(count_matrix+.5)
}
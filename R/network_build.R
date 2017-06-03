# builds a network based on correlation (default=spearman)
# return network and laplacian of the network

network_build <- function(data_train, p, corr) {

  cor.method = "spearman"
  network_corr <- cor(data_train,method = cor.method)


  num_sample = nrow(data_train) - 2 # degrees of freedom

  # T-test w/ H0: correlation == 0; H1: correlation != 0
  network_corr_t <-
    network_corr / (sqrt(abs(1 - (network_corr ^ 2)) / (num_sample))) # t-stat for T-test
  network_corr_t_p <-
    2 * (1 - pt(abs(network_corr_t), num_sample)) # p-value for T-test

  # FDR adjustment for multiple hypothesis testing
  network_p <-
    matrix(
      p.adjust(network_corr_t_p, 'fdr'), # adjusted p-values (as a vector)
      nrow = nrow(network_corr_t_p),
      ncol = ncol(network_corr_t_p)
    ) # p-values as a matrix

  # p-value filter
  network <- network_corr
  if(p>0){
    network[network_p<min(p,max(network_p))] = 0
  }

  # correlation value filter
  if(corr=='plus'){
    network[network < 0] <- 0
  }

  # graph laplacian
  require('igraph', quietly = TRUE)
  g <- graph.adjacency(network,mode = "upper",weighted = TRUE,diag = FALSE)
  b <- as.matrix(laplacian_matrix(g))

  return(list(network = network, laplacian = b))
}


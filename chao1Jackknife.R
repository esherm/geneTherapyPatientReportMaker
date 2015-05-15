getPopEstimates <- function(counts) {
  if(all(dim(counts)>2) & is.matrix(counts)) {
    # do jackknife correction if more than 1 samples/replicates are available
    pseudo <- ncol(counts)*
      estimateR(rowSums(counts))-(ncol(counts)-1)*
      sapply(1:ncol(counts), function(y) estimateR(rowSums(counts[,-y])))
    round(rowMeans(pseudo)["S.chao1"])
  } else {
    pseudo <- estimateR(rowSums(counts))
    round(pseudo["S.chao1"])
  }
}
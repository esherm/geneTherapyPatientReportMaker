getAbundanceSums <- function(splitData, abundCutoff){
  do.call(rbind, lapply(seq(length(splitData)), function(x){
    sites <- splitData[[x]]
    sites$maskedRefGeneName <- ifelse(sites$estAbundProp > abundCutoff,
                                      sites$nearest_refSeq_gene,
                                      "LowAbund")
    
    res <- aggregate(estAbundProp~maskedRefGeneName, mcols(sites), sum)
    res$Timepoint <- strsplit(names(splitData)[x], ":")[[1]][1]
    res$CellType <- strsplit(names(splitData)[x], ":")[[1]][2]
    
    res
  }))
}
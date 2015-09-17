getAbundanceThreshold <- function(sites, numGenes){
  orderedAbundances <- sites[order(-sites$estAbundProp)]
  #'unique' functionally removes anything that's duplicated, thus preserving order
  abundCutoff <- orderedAbundances$nearest_refSeq_gene==tail(head(unique(orderedAbundances$nearest_refSeq_gene), numGenes),1)
  abundCutoff <- orderedAbundances[which(abundCutoff)[1]]$estAbundProp  
  abundCutoff
}

filterLowAbund <- function(sites, abundCutoff){
  sites$maskedRefGeneName <- ifelse(sites$estAbundProp > abundCutoff,
                                    sites$nearest_refSeq_gene,
                                    "LowAbund")
  sites
}

getAbundanceSums <- function(sites, splitBy){
  splitBy <- mcols(sites)[,splitBy]
  splitSites <- split(sites, apply(as.data.frame(splitBy), 1, paste, collapse=""))
  
  do.call(rbind, lapply(splitSites, function(sites){
    res <- aggregate(estAbundProp~maskedRefGeneName, mcols(sites), sum)
    res$Timepoint <- sites[1]$Timepoint
    res$CellType <- sites[1]$CellType
    
    res
  }))
}

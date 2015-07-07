getEstimatedAbundance <- function(sites, use.sonicLength=FALSE){
  #inputs at this point have been standardized, so it's ok to use posid as PK
  if(use.sonicLength){
    estAbund.uniqueFragLen <- function(location, fragLen, replicate=NULL){
      if(is.null(replicate)){replicate <- 1}  #Need for downstream workflow
      dfr <- data.frame(location = location, fragLen = fragLen, 
                      replicate = replicate)
      dfr <- distinct(dfr)
      site.list <- split(dfr, dfr$location)
      theta <- sapply(site.list, function(x){nrow(x)})
      list(theta=theta)
    }
    estAbund <- estAbund.uniqueFragLen
  }
  
  sites$posid = paste0(seqnames(sites), strand(sites), start(flank(sites, -1, start=T)))
  dfr <- data.frame("ID"=sites$posid,
                    "fragLength"=width(sites),
                    "replicate"=sites$replicate)
  
  #sonic abundance will crash if given a vector of replicates with one value
  if(length(unique(dfr$replicate)) == 1){
    estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength)
  }else{
    estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength, dfr$replicate)
  }
  
  dereplicatedSites = granges(dereplicateSites(sites)) #instead of mcols()=NULL

  #regenerate posid which is guaranteed to match since the input data was
  #already run through the dereplicator
  dereplicatedSites$posid = paste0(seqnames(dereplicatedSites),
                                   strand(dereplicatedSites),
                                   start(flank(dereplicatedSites, -1, start=T)))
  
  dereplicatedSites$estAbund <- round(estimatedAbundances$theta) #estAbund preserves order
  dereplicatedSites$estAbundProp <- dereplicatedSites$estAbund/sum(dereplicatedSites$estAbund)
  dereplicatedSites$estAbundRank <- rank(-1*dereplicatedSites$estAbundProp, ties.method="max")
  
  dereplicatedSites
}

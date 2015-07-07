getEstimatedAbundance <- function(sites){
  #inputs at this point have been standardized, so it's ok to use posid as PK
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


getEstimatedAbundance.chris <- function(sites){
  #standardized input is converted to df and pcr dereplicated
  dfr <- as.data.frame(sites)
  dfr <- distinct(dfr)
  
  #split based on unique start position
  sites.list <- split(dfr, dfr$posid)
  
  #collect relevant info for abundance and reconstruction
  unique.sites.df <- do.call(rbind, lapply(sites.list, function(x){
    intSite <- ifelse(x$strand[1] == "+", x$start[1], x$end[1])
    breakpoint <- ifelse(x$strand[1] == "+", max(x$end), min(x$start))
    df <- data.frame(seqnames = x$seqnames[1], 
                     start = ifelse(x$strand[1] == "+", intSite, breakpoint), 
                     end = ifelse(x$strand[1] == "+", breakpoint, intSite), 
                     strand = x$strand[1], 
                     GTSP = x$GTSP[1], 
                     posid = x$posid[1], 
                     estAbund = nrow(x))
  }))
  row.names(unique.sites.df) <- NULL
  
  ranges <- IRanges(start = unique.sites.df$start, end = unique.sites.df$end)
  
  unique.sites <- GRanges(
    seqnames = unique.sites.df$seqnames,
    ranges = ranges,
    strand = unique.sites.df$strand,
    seqinfo = seqinfo(sites),
    GTSP = unique.sites.df$GTSP,
    posid = unique.sites.df$posid,
    estAbund = unique.sites.df$estAbund)
  
  #add proportion and rank to dereplicated unique grange
  unique.sites$estAbundProp <- unique.sites$estAbund/sum(unique.sites$estAbund)
  unique.sites$estAbundRank <- rank(-1*unique.sites$estAbundProp, ties.method="max")
  
  unique.sites
}
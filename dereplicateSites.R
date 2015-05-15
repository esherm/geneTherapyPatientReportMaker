dereplicateSites <- function(uniqueSites){
  if(length(uniqueSites)>0){
    sites.reduced <- reduce(flank(uniqueSites, -5, both=TRUE), with.revmap=T)
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    allSites <- uniqueSites[unlist(sites.reduced$revmap)]
    allSites <- split(allSites, Rle(values=c(1:length(sites.reduced)), lengths=sites.reduced$counts))
    allSites <- unlist(reduce(allSites, min.gapwidth=5))
    mcols(allSites) <- mcols(sites.reduced)
    allSites
  }else{
    uniqueSites
  }
}
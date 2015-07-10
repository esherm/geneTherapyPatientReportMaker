dereplicateSites <- function(sites){
  #Reduce sites which have the same starts, but loose range info
  #(no need to add a gapwidth as sites are standardized)
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced <- unlist(reduce(sites.reduced, with.revmap=TRUE))
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
  #Order original sites by revmap  
  dereplicatedSites <- sites[unlist(sites.reduced$revmap)]
  
  #Skip this step and provide similar output if length(sites) = 0
  if(length(sites) > 0){
    dereplicatedSites <- split(dereplicatedSites, Rle(values = seq(length(sites.reduced)), lengths = sites.reduced$counts))
  }  
  
  #Dereplicate reads with same standardized starts and provide the longeset width
  dereplicatedSites <- unlist(reduce(dereplicatedSites))
  mcols(dereplicatedSites) <- mcols(sites.reduced)
  
  dereplicatedSites
}
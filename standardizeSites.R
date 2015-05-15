standardizeSites <- function(unstandardizedSites){
  
  dereplicated <- dereplicateSites(unstandardizedSites)
  dereplicated$dereplicatedSiteID <- seq(length(dereplicated))
  
  #expand, keeping the newly standardized starts
  standardized <- unname(dereplicated[rep(dereplicated$dereplicatedSiteID, dereplicated$counts)])
  
  #order the original object to match
  unstandardizedSites <- unstandardizedSites[unlist(dereplicated$revmap)]
  
  #graft over the widths and metadata
  trueBreakpoints <- start(flank(unstandardizedSites, -1, start=F))
  standardizedStarts <- start(flank(standardized, -1, start=T))
  standardized <- GRanges(seqnames=seqnames(standardized),
                         ranges=IRanges(start=pmin(standardizedStarts, trueBreakpoints),
                                        end=pmax(standardizedStarts, trueBreakpoints)),
                         strand=strand(standardized))
  mcols(standardized) <- mcols(unstandardizedSites)
  
  standardized
}
dereplicateSites <- function(uniqueReads){
    #do the dereplication, but loose the coordinates
    sites.reduced <- reduce(flank(uniqueReads, -5, both=TRUE), with.revmap=T)
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    #order the unique sites as described by revmap
    dereplicatedSites <- uniqueReads[unlist(sites.reduced$revmap)]
    
    #split the unique sites as described by revmap (sites.reduced$counts came from revmap above)
    dereplicatedSites <- split(dereplicatedSites, Rle(values=seq(length(sites.reduced)), lengths=sites.reduced$counts))
        
    #do the standardization - this will pick a single starting position and
    #choose the longest fragment as ending position
    dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth=5))
    mcols(dereplicatedSites) <- mcols(sites.reduced)
    
    dereplicatedSites
}


standardizeSites <- function(unstandardizedSites){
  if(length(unstandardizedSites) > 0){
    #Known conflict with package:dplyr::count(), need to unload package if present
    isthere <- which("dplyr" == loadedNamespaces())
    if(length(isthere) > 0){detach("package:dplyr", unload = TRUE)}
    
    #Get called start values for clustering  
    unstandardizedSites$Position <- ifelse(strand(unstandardizedSites) == "+", start(unstandardizedSites), end(unstandardizedSites))
    unstandardizedSites$Break <- ifelse(strand(unstandardizedSites) == "+", end(unstandardizedSites), start(unstandardizedSites))
    unstandardizedSites$Score <- 95
    unstandardizedSites$qEnd <- width(unstandardizedSites)
    
    #Positions clustered by 5L window and best position is chosen for cluster
    standardized <- clusterSites(
      psl.rd = unstandardizedSites,
      weight = rep(1, length(unstandardizedSites)) 
    )
    
    start(standardized) <- ifelse(strand(standardized) == "+", 
                                  standardized$clusteredPosition, standardized$Break)
    end(standardized) <- ifelse(strand(standardized) == "-", 
                                standardized$clusteredPosition, standardized$Break)
    
    standardized$Position <- NULL
    standardized$Break <- NULL
    standardized$score <- NULL
    standardized$qEnd <- NULL
    standardized$clusteredPosition <- NULL
    standardized$clonecount <- NULL
    standardized$clusterTopHit <- NULL
    
    if(length(isthere) > 0){suppressMessages(library("dplyr"))}
    
    return(sort(standardized))
  }else{
    unstandardizedSites
  }
}
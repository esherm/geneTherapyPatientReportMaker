# access patient metadata for GTSP

#' @param GTSP vector of GTSP 
#' @return df with columns: Trial, GTSP, Patient, Timepoint, CellType, FragMethod, VCN.
#' @note all GTSP in th vector should be unique
get_metadata_for_GTSP <- function(GTSP, db_group) {
    stopifnot(length(GTSP) == length(unique(GTSP)))
    GTSP = paste0("^", GTSP, "$", collapse="|")
    
    GTSPDBconn <- dbConnect(MySQL(), group=db_group)
    
    query = paste0("SELECT Trial, SpecimenAccNum, Patient, Timepoint, CellType,
                           SamplePrepMethod, VCN
                   FROM gtsp
                   WHERE SpecimenAccNum
                   REGEXP ", dbQuoteString(GTSPDBconn, GTSP), ";")
    
    sets <- dbGetQuery(GTSPDBconn, query)
    dbDisconnect(GTSPDBconn)
    
    names(sets) <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", "FragMethod", "VCN")
    sets    
}

# access patient metadata for GTSP

#' @param GTSP vector of GTSP 
#' @return df with columns: Trial, GTSP, Patient, Timepoint, CellType, FragMethod, VCN.
#' @note all GTSP in th vector should be unique
get_metadata_for_GTSP <- function(GTSP) {
    #TODO: should repeated GTSP be allowed?
    stopifnot(length(GTSP) == length(unique(GTSP)))
    # TODO: get only what is needed from the DB
    GTSPDBconn <- dbConnect(MySQL(), group="specimen_management")
    sets <- dbGetQuery(GTSPDBconn, "SELECT Trial, SpecimenAccNum, Patient, Timepoint, CellType, 
                                          SamplePrepMethod, VCN
                                   FROM gtsp;")
    names(sets) <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", "FragMethod", "VCN")
    dbDisconnect(GTSPDBconn)
    sets <- sets[ sets$GTSP %in% GTSP, ]
}

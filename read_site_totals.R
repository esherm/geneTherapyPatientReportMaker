source("intSiteRetriever/intSiteRetriever.R")

#' for a given df with 2 cols: sampleName and GTSP
#' gets counts from intSite DB and collapse replicates.
#' @param df with cols: sampleName, GTSP
get_read_site_totals <- function(sampleName_GTSP) {
    reads <- get_count_per_GTSP(sampleName_GTSP, getReadCounts, "TotalReads")
    sites <- get_count_per_GTSP(sampleName_GTSP, getUniqueSiteCounts, "UniqueSites")
    merge(reads, sites)
}


#' call database_function to get counts then collapse all samples with the same GTSP
#' numbers.
#' @return df with 2 cols: GTSP, column_name
get_count_per_GTSP <- function(sampleName_GTSP, database_function, column_name) {
    sample_count <- database_function(sampleName_GTSP$sampleName)
    names(sample_count) <- c("sampleName", column_name)
    sample_count_GTSP <- merge(sample_count, sampleName_GTSP)
    aggregate(sample_count_GTSP[column_name], sample_count_GTSP['GTSP'], FUN=sum)
}

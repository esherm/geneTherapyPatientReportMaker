#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
source("intSiteRetriever/intSiteRetriever.R")
source("CancerGeneList/onco_genes.R")
source("utilities.R")

#INPUTS: csv file/table GTSP to sampleName
sampleName_GTSP = read.csv("sampleName_GTSP.csv")

dbConn <- .connectToDB(NULL)
stopifnot(all(setNameExists(sampleName_GTSP$sampleName)))

read_counts <- getReadCounts(sampleName_GTSP$sampleName)
sites_counts <- getUniqueSiteCounts(sampleName_GTSP$sampleName)

#grab all possible GTSP samples, and pick the ones we care about for this report
GTSPDBconn = dbConnect(MySQL(), group="specimen_management")
sets = dbGetQuery(GTSPDBconn, "SELECT SpecimenAccNum, Trial, CellType, Patient,
                                      Timepoint, SamplePrepMethod, VCN
                               FROM gtsp;")
dbDisconnect(GTSPDBconn)

sets = merge(sets, sampleName_GTSP, by.x="SpecimenAccNum" , by.y="GTSP")

stopifnot(length(unique(sets$Patient)) == 1) # reporst are for a single patient

sets <- merge(sets, read_counts)
sets <- merge(sets, sites_counts)

# 
refGenomes = getRefGenome(sampleName_GTSP$sampleName)
sets = merge(sets, refGenomes)
sets$timepointDay <- mdy_to_day(sets$Timepoint)

sites = getUniquePCRbreaks(sets$sampleName)

###GET MULTIHITS EVENTUALLY###




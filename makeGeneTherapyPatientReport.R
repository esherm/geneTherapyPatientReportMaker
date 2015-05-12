#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
source("oldReports/helper_functions.R")
source("../intSiteRetriever/intSiteRetriever.R")
source("../genomicHeatmapMaker/CancerGeneList/onco_genes.R")

#INPUTS: csv file/table GTSP to sampleName
sampleName_GTSP = read.csv("sampleName_GTSP.csv")
stopifnot(all(setNameExists(sampleName_GTSP$sampleName)))

#grab all possible GTSP samples, and pick the ones we care about for this report
GTSPDBconn = dbConnect(MySQL(), group="specimen_management")
sets = dbGetQuery(GTSPDBconn, "SELECT SpecimenAccNum, Trial, CellType, Patient,
                                      Timepoint, SamplePrepMethod, VCN
                               FROM gtsp;")
dbDisconnect(GTSPDBconn)

sets = merge(sets, sampleName_GTSP, by.x="SpecimenAccNum" , by.y="GTSP")
refGenomes = getRefGenome(sampleName_GTSP$sampleName)
sets = merge(sets, refGenomes)
sets$timepointDay <- timepointDays(sets$timepoint)

sites = getUniquePCRbreaks(sets$sampleName)

###GET MULTIHITS EVENTUALLY###




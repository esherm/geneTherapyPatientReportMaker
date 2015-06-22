library("methods", quietly=TRUE)
library("RMySQL", quietly = TRUE) #also loads DBI

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237")

sql <- "select trial, patient, CellType, Timepoint, SpecimenAccNum from specimen_management.gtsp"
##message(sql)
trial_pat_gtsp <- dbGetQuery(dbConn,sql)
names(trial_pat_gtsp) <- tolower(names(trial_pat_gtsp))

workDir <- commandArgs(trailingOnly=TRUE)[1]
if( is.na(workDir) ) workDir <- "."
workDir <- normalizePath(workDir, mustWork=TRUE)
gtspDir <- grep("^GTSP", list.dirs(path=workDir, full.names=FALSE), ignore.case=TRUE, value=TRUE)

gtsp <- sub("-\\d+$", "", gtspDir)

df <- data.frame(GTSP=gtsp, sampleName=gtspDir)

merged.tab <- merge(df, trial_pat_gtsp,
                    by.x="GTSP", by.y="specimenaccnum")
write.table(merged.tab, "", sep = "\t", row.names=FALSE, quote=FALSE)

q()


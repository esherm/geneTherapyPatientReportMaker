#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("markdown")
library("knitr")

source("intSiteRetriever/intSiteRetriever.R")
source("CancerGeneList/onco_genes.R")
source("utilities.R")
source("specimen_management.R")
source("read_site_totals.R")

#INPUTS: csv file/table GTSP to sampleName
sampleName_GTSP = read.csv("sampleName_GTSP.csv")
stopifnot(all(setNameExists(sampleName_GTSP$sampleName)))

read_sites_sample_GTSP <- get_read_site_totals(sampleName_GTSP)

sets <- get_metadata_for_GTSP(unique(sampleName_GTSP$GTSP))
# reports are for a single patient
stopifnot(length(unique(sets$Patient)) == 1)
# all GTSP in the database
stopifnot(nrow(sets) == length(unique(sampleName_GTSP$GTSP)))
#end INPUTS 

sets <- merge(sets, read_sites_sample_GTSP)

refGenomes = getRefGenome(sampleName_GTSP$sampleName)
# at present the whole report is done for one genome
stopifnot(length(unique(refGenomes$refGenome))==1)

sets$timepointDay <- mdy_to_day(sets$Timepoint)

sites = getUniquePCRbreaks(sets$sampleName)
###GET MULTIHITS EVENTUALLY###

# MAGIC HERE

### PRINT REPORT

#set variables for markdown report
patient = sets$Patient[1]
freeze = refGenomes[1, 2]
timepoint = sort(unique(sets$timepointDay))

cols <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN")
summaryTable <- arrange(sets,timepointDay,CellType)                
summaryTable <- summaryTable[,cols]

#end setting variables for markdown report
  
filename = "report.md"
outFilename <- gsub("\\.md",".html",filename)
options(knitr.table.format = 'html')
knit("GTSPreport.Rmd", output=filename)
markdownToHTML(filename, outFilename, extensions=c('tables'),
    options=c(markdownHTMLOptions(defaults=T),"toc"),
    stylesheet="GTSPreport.css") 




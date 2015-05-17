#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("markdown")
library("knitr")

source("intSiteRetriever/intSiteRetriever.R")
source("CancerGeneList/onco_genes.R")
source("utilities.R")
source("specimen_management.R")
source("estimatedAbundance.R")
source("chao1Jackknife.R")
source("dereplicateSites.R")
source("standardizeSites.R")
source("read_site_totals.R")

#INPUTS: csv file/table GTSP to sampleName
sampleName_GTSP <- read.csv("sampleName_GTSP.csv")
GTSPs <- unique(sampleName_GTSP$GTSP)

stopifnot(all(setNameExists(sampleName_GTSP$sampleName)))

read_sites_sample_GTSP <- get_read_site_totals(sampleName_GTSP)

sets <- get_metadata_for_GTSP(unique(sampleName_GTSP$GTSP))
# reports are for a single patient
stopifnot(length(unique(sets$Patient)) == 1)
# all GTSP in the database
stopifnot(nrow(sets) == length(unique(sampleName_GTSP$GTSP)))
#end INPUTS 

sets <- merge(sets, read_sites_sample_GTSP)

refGenomes <- getRefGenome(sampleName_GTSP$sampleName)
# at present the whole report is done for one genome
stopifnot(length(unique(refGenomes$refGenome))==1)

sets$timepointDay <- mdy_to_day(sets$Timepoint)

sites <- merge(getUniquePCRbreaks(sets$sampleName), sampleName_GTSP)

#we really don't care about seqinfo - we just want a GRange object for easy manipulation
uniqueSites.gr <- GRanges(seqnames=Rle(sites$chr),
                          ranges=IRanges(start=pmin(sites$integration, sites$breakpoint),
                                         end=pmax(sites$integration, sites$breakpoint)),
                          strand=Rle(sites$strand))
mcols(uniqueSites.gr) <- sites[,c("sampleName", "GTSP")]

#standardize sites across GTSP
standardizedReplicatedSites <- lapply(split(uniqueSites.gr, uniqueSites.gr$GTSP), function(sites){
  res <- standardizeSites(sites)
  res$replicate <- as.integer(as.factor(res$sampleName))
  res$sampleName <- NULL
  res$posid <- paste0(seqnames(res), strand(res), start(flank(res, -1, start=T)))
  res
})

#this is slow (~1.5min/sample), but would be easy to parallelize - just be
#careful about memory consumption!  sonic abundance could get 20GB+ per thread
standardizedDereplicatedSites <- lapply(standardizedReplicatedSites, function(sites){
  res <- getEstimatedAbundance(sites)
  res$GTSP <- sites[1]$GTSP
  res$posid <- paste0(seqnames(res), strand(res), start(flank(res, -1, start=T)))
  res
})

populationInfo <- lapply(GTSPs, function(GTSP){
  #can iterate through standardizedReplicatedSites and standardizedDereplicatedSites using GTSP#
  replicatedSites <- standardizedReplicatedSites[[GTSP]]
  dereplicatedSites <- standardizedDereplicatedSites[[GTSP]]

  #in preperation for cast below, this needs to be a base R data.frame
  merged <- as.data.frame(merge(mcols(replicatedSites[,c("posid", "replicate")]),
                               mcols(dereplicatedSites[,c("posid", "estAbund")])))

  replicatesByPosid <- acast(merged, formula=posid~replicate, fun.aggregate=sum, value.var="estAbund", fill=0)

  data.frame("GTSP"=GTSP,
             "S.chao1"=getPopEstimates(replicatesByPosid),
             #"Gini"=round(gini(standardizedDereplicatedSites[[GTSP]]$estAbundProp), 2),
             "Gini"=gini(standardizedDereplicatedSites[[GTSP]]$estAbundProp),
             "Replicates"=max(replicatedSites$replicate))
})

populationInfo <- do.call(rbind, populationInfo)
rownames(populationInfo) <- NULL

#standardizedSites <- unname(unlist(GRangesList(standardizedSites)))

###GET MULTIHITS EVENTUALLY###

# MAGIC HERE

### PRINT REPORT

#set variables for markdown report

patient <- sets$Patient[1]
freeze <- refGenomes[1, "refGenome"]
timepoint <- sort(unique(sets$timepointDay))

cols <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN")
summaryTable <- arrange(sets,timepointDay,CellType)
summaryTable <- summaryTable[,cols]

cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
          "Replicates", "FragMethod", "VCN", "S.chao1", "Gini")
popSummaryTable <- merge(sets,  populationInfo)
popSummaryTable <- arrange(popSummaryTable,timepointDay,CellType)
popSummaryTable <- popSummaryTable[,cols]

#end setting variables for markdown report
  
filename <- "report.md"
outFilename <- gsub("\\.md",".html",filename)
options(knitr.table.format='html')
knit("GTSPreport.Rmd", output=filename)
markdownToHTML(filename, outFilename, extensions=c('tables'),
    options=c(markdownHTMLOptions(defaults=T),"toc"),
    stylesheet="GTSPreport.css") 


#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("markdown")
library("knitr")
library("hiAnnotator")
library("ggplot2")

source("intSiteRetriever/intSiteRetriever.R")
source("CancerGeneList/onco_genes.R")
source("utilities.R")
source("specimen_management.R")
source("estimatedAbundance.R")
source("dereplicateSites.R")
source("standardizeSites.R")
source("read_site_totals.R")
source("populationInfo.R")
source("annotatedAbundanceSums.R")

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

#==========GET AND PERFORM BASIC DEREPLICATION/SONICABUND ON SITES=============
sites <- merge(getUniquePCRbreaks(sampleName_GTSP$sampleName), sampleName_GTSP)

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

standardizedReplicatedSites <- unname(unlist(GRangesList(standardizedReplicatedSites)))
mcols(standardizedReplicatedSites) <- merge(as.data.frame(mcols(standardizedReplicatedSites)),
                                           sets[,c("GTSP", "Timepoint",
                                                   "CellType", "timepointDay")])

standardizedDereplicatedSites <- unname(unlist(GRangesList(standardizedDereplicatedSites)))
mcols(standardizedDereplicatedSites) <- merge(as.data.frame(mcols(standardizedDereplicatedSites)),
                                           sets[,c("GTSP", "Timepoint",
                                                   "CellType", "timepointDay")])

#============CALCULATE POPULATION SIZE/DIVERSITY INFORMATION=================
populationInfo <- getPopulationInfo(standardizedReplicatedSites,
                                    standardizedDereplicatedSites,
                                    "GTSP")
populationInfo$Replicates <- sapply(split(standardizedReplicatedSites$replicate,
                                          standardizedReplicatedSites$GTSP),
                                    max)

#========CALCULATE POPULATION SIZE/DIVERSITY INFORMATION BY TIMEPOINT==========
timepointPopulationInfo <- getPopulationInfo(standardizedReplicatedSites,
                                             standardizedDereplicatedSites,
                                             "Timepoint")

timepointPopulationInfo$UniqueSites <- sapply(split(standardizedDereplicatedSites,
                                                    standardizedDereplicatedSites$Timepoint),
                                              length)


#==================ANNOTATE AND AGGREGATE ABUNDANCES=====================
refSeq_genes <- makeGRanges(
  getUCSCtable("refGene", "RefSeq Genes", freeze=reference_genome),
  freeze=unique(refGenomes$refGenome) #there will only be one as enforced above
)

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   feature.colnam="name2")

#used for calling sites as 'lowAbund'
abundCutoff.barplots <- 0.03

dereplicatedTimepointCellType <- split(standardizedDereplicatedSites, paste0(standardizedDereplicatedSites$Timepoint,
                                                                             ":",
                                                                             standardizedDereplicatedSites$CellType))

barplotAbunds <- getAbundanceSums(dereplicatedTimepointCellType, abundCutoff.barplots)


#tweak abundCutoff
orderedAbundances <- standardizedDereplicatedSites[order(-standardizedDereplicatedSites$estAbundProp)]
#'unique' functionally removes anything that's duplicated, thus preserving order
abundCutoff.detailed <- orderedAbundances$nearest_refSeq_gene==unique(orderedAbundances$nearest_refSeq_gene)[50]
abundCutoff.detailed <- orderedAbundances[which(abundCutoff.detailed)[1]]$estAbundProp

detailedAbunds <- getAbundanceSums(dereplicatedTimepointCellType, abundCutoff.detailed)




###GET MULTIHITS EVENTUALLY###


#set variables for markdown report

patient <- sets$Patient[1]
freeze <- refGenomes[1, "refGenome"]
timepoint <- sort(unique(sets$timepointDay))

#used for calling sites as 'lowAbund'
abundCutoff <- 0.03

cols <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN")
summaryTable <- arrange(sets,timepointDay,CellType)
summaryTable <- summaryTable[,cols]

cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
          "Replicates", "FragMethod", "VCN", "S.chao1", "Gini", "Shannon")
popSummaryTable <- merge(sets,  populationInfo)
popSummaryTable <- arrange(popSummaryTable,timepointDay,CellType)
popSummaryTable <- popSummaryTable[,cols]

timepointPopulationInfo <- melt(timepointPopulationInfo, "group")


#end setting variables for markdown report
  
filename <- "report.md"
outFilename <- gsub("\\.md",".html",filename)
options(knitr.table.format='html')
#knit("GTSPreport.Rmd", output=filename)
theme_set(theme_bw()) #for ggplot2
knit("shortGTSPReport.Rmd", output=filename)
markdownToHTML(filename, outFilename, extensions=c('tables'),
    options=c(markdownHTMLOptions(defaults=T),"toc"),
    stylesheet="GTSPreport.css") 


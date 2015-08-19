options(stringsAsFactors = FALSE)

library(argparse)

parser <- ArgumentParser(description="Gene Therapy Patient Report for Single Patient")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-s", action='store_true', help="abundance by sonicLength package (Berry, C. 2012)")
parser$add_argument("--ref_genome", default="hg18", help="reference genome used for all samples")
parser$add_argument("--sites_group", default="intsites_miseq", help="group to use for integration sites db from ~/.my.cnf")
parser$add_argument("--gtsp_group", default="specimen_management", help="group to use for specimen management GTSP db from ~/.my.cnf")
arguments <- parser$parse_args()

args <- commandArgs(trailingOnly=TRUE)

# defaults:
use.sonicLength <-  ! arguments$s
db_group_sites <- arguments$sites_group
db_group_gtsp <- arguments$gtsp_group
ref_genome <- arguments$ref_genome
#### INPUTS: csv file/table GTSP to sampleName ####
csvfile <- arguments$sample_gtsp

if( !file.exists(csvfile) ) stop(csvfile, "not found")

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="geneTherapyPatientReportMaker$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
stopifnot(file.exists(file.path(codeDir, "GTSPreport.css")))
stopifnot(file.exists(file.path(codeDir, "GTSPreport.Rmd")))

#### load up require packages + objects #### 
library("RMySQL") #also loads DBI
library("plyr")
library("dplyr")
library("stringr")
library("markdown")
library("knitr")
library("PubMedWordcloud")
library("hiAnnotator")
library("ggplot2")
library("reldist")
library("sonicLength")
library("reshape2")
library("scales")
library("intSiteRetriever")
library("BiocParallel")
library(devtools)

source(file.path(codeDir, "utilities.R"))
source(file.path(codeDir, "specimen_management.R"))
source(file.path(codeDir, "estimatedAbundance.R"))
source(file.path(codeDir, "read_site_totals.R"))
source(file.path(codeDir, "populationInfo.R"))
source(file.path(codeDir, "abundanceFilteringUtils.R"))

url <- "https://raw.githubusercontent.com/BushmanLab/intSiteCaller/master/"
source_url(paste0(url, "hiReadsProcessor.R"))
source_url(paste0(url, "standardization_based_on_clustering.R"))

#### load datasets and process them before knit #### 
message("\nReading csv from ", csvfile)
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

dbConn <- dbConnect(MySQL(), group=db_group_sites)
info <- dbGetInfo(dbConn)
dbConn <- src_sql("mysql", dbConn, info = info)

sampleName_GTSP$refGenome <- rep(ref_genome, nrow(sampleName_GTSP))

stopifnot(all(setNameExists(sampleName_GTSP, dbConn)))

read_sites_sample_GTSP <- get_read_site_totals(sampleName_GTSP, dbConn)

sets <- get_metadata_for_GTSP(unique(sampleName_GTSP$GTSP), db_group_gtsp)
# reports are for a single patient
stopifnot(length(unique(sets$Patient)) == 1)
patient <- sets$Patient[1]
# and for a single trial
stopifnot(length(unique(sets$Trial)) == 1)
trial <- sets$Trial[1]

RDataFile <- paste(trial, patient, format(Sys.Date(), format="%Y%m%d"), "RData", sep=".")

# all GTSP in the database
stopifnot(nrow(sets) == length(unique(sampleName_GTSP$GTSP)))

##end INPUTS
sets[sets$Timepoint=="NULL", "Timepoint"] = "d0"

sets <- merge(sets, read_sites_sample_GTSP)
sets$Timepoint <- sortFactorTimepoints(sets$Timepoint)

# at present the whole report is done for one genome
stopifnot(length(unique(sampleName_GTSP$refGenome))==1)
freeze <- sampleName_GTSP[1, "refGenome"]

##==========GET AND PERFORM BASIC DEREPLICATION/SONICABUND ON SITES=============
message("Fetching unique sites and estimating abundance")
dbConn <- dbConnect(MySQL(), group=db_group_sites)
info <- dbGetInfo(dbConn)
dbConn <- src_sql("mysql", dbConn, info = info)
sites <- merge(getUniquePCRbreaks(sampleName_GTSP, dbConn), sampleName_GTSP)
names(sites)[names(sites)=="position"] <- "integration"

#we really don't care about seqinfo - we just want a GRange object for easy manipulation
uniqueSites.gr <- GRanges(seqnames=Rle(sites$chr),
                          ranges=IRanges(start=pmin(sites$integration, sites$breakpoint),
                                         end=pmax(sites$integration, sites$breakpoint)),
                          strand=Rle(sites$strand))
mcols(uniqueSites.gr) <- sites[,c("sampleName", "GTSP")]

#standardize sites across all GTSPs
isthere <- which("dplyr" == loadedNamespaces()) # temp work around  of
#Known conflict with package:dplyr::count(), need to unload package if present
# unloading and reloading the package
if(length(isthere) > 0){detach("package:dplyr", unload = TRUE)}
standardizedReplicatedSites <- standardizeSites(uniqueSites.gr)
if(length(isthere) > 0){suppressMessages(library("dplyr"))}

standardizedReplicatedSites$posid <- paste0(seqnames(standardizedReplicatedSites),
                                            strand(standardizedReplicatedSites),
                                            start(flank(standardizedReplicatedSites, -1, start=T)))
standardizedReplicatedSites <- split(standardizedReplicatedSites,
                                     standardizedReplicatedSites$GTSP)
standardizedReplicatedSites <- lapply(standardizedReplicatedSites, function(x){
    x$replicate <- as.integer(as.factor(x$sampleName))
    x$sampleName <- NULL
    x
})

#this is slow (~1.5min/sample), but would be easy to parallelize - just be
#careful about memory consumption!  sonic abundance could get 20GB+ per thread
standardizedReplicatedSites <- standardizedReplicatedSites[sapply(standardizedReplicatedSites, length)>0]

standardizedDereplicatedSites <- lapply(standardizedReplicatedSites, function(sites){
    res <- getEstimatedAbundance(sites, use.sonicLength = use.sonicLength)
    res$GTSP <- sites[1]$GTSP
    res$posid <- paste0(seqnames(res), strand(res), start(flank(res, -1, start=T)))
    res
})

standardizedReplicatedSites <- prepSiteList(standardizedReplicatedSites)
standardizedDereplicatedSites <- prepSiteList(standardizedDereplicatedSites)
standardizedDereplicatedSites <- flank(standardizedDereplicatedSites, -1, start=TRUE)

unique_sites_per_GTSP <- sapply(split(standardizedDereplicatedSites,
                                      standardizedDereplicatedSites$GTSP),
                                function(x){length(unique(x$posid))})
unique_sites_per_GTSP <- data.frame("GTSP" = names(unique_sites_per_GTSP),
                                    "UniqueSites" = unique_sites_per_GTSP)
sets <- merge(sets, unique_sites_per_GTSP, by = "GTSP")

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
                                              function(x){length(unique(x$posid))})


#=======================ANNOTATE DEREPLICATED SITES==========================
#standard refSeq genes
message("Annotating unique hit sites")
refSeq_genes <- makeGRanges(
  getUCSCtable("refGene", "RefSeq Genes", freeze=freeze),
  freeze=freeze
)
save.image(RDataFile)

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   feature.colnam="name2")

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="nearest_refSeq_gene",
                                                   side="5p",
                                                   feature.colnam="name2")

standardizedDereplicatedSites <- getSitesInFeature(standardizedDereplicatedSites,
                                                   refSeq_genes,
                                                   colnam="inGene",
                                                   feature.colnam="name2")

#oncogenes
oncogenes <- scan(file= file.path(codeDir, "allonco_no_pipes.csv"), what='charactor')
oncogenes <- oncogenes[!grepl("geneName", oncogenes, ignore.case=TRUE)]

refSeq_oncogene <- refSeq_genes[toupper(refSeq_genes$name2) %in% toupper(oncogenes)]

standardizedDereplicatedSites <- getNearestFeature(standardizedDereplicatedSites,
                                                   refSeq_oncogene,
                                                   colnam="NrstOnco",
                                                   side="5p",
                                                   feature.colnam="name2")


wantedgenes <- scan(file=file.path(codeDir, "genes_adverse_event.csv"), what='charactor')
wantedgenes <- wantedgenes[!grepl("geneName", wantedgenes, ignore.case=TRUE)]

## * in transcription units
## ~ within 50kb of a onco gene 
## ! nearest is a known bad gene 
standardizedDereplicatedSites$geneMark <- ""

## ~ nearest is a known bad gene 
isNearWanted <- standardizedDereplicatedSites$nearest_refSeq_gene %in% wantedgenes 
isInWanted <- sapply( standardizedDereplicatedSites$inGene,
                     function(txt) any(unlist(strsplit(txt, ',')) %in% wantedgenes) )
standardizedDereplicatedSites$geneMark <- ifelse(
    isNearWanted | isInWanted,
    paste0(standardizedDereplicatedSites$geneMark, "!"),
    standardizedDereplicatedSites$geneMark )

## ~ within 50kb of a onco gene 
isNearOnco <- (standardizedDereplicatedSites$nearest_refSeq_gene %in% oncogenes &
               abs(standardizedDereplicatedSites$nearest_refSeq_geneDist) < 50000 )
isInOnco <- sapply( standardizedDereplicatedSites$inGene,
                   function(txt) any(unlist(strsplit(txt, ',')) %in% oncogenes) )
standardizedDereplicatedSites$geneMark <- ifelse(
    isNearOnco | isInOnco,
    paste0(standardizedDereplicatedSites$geneMark, "~"),
    standardizedDereplicatedSites$geneMark )

## * in transcription units
standardizedDereplicatedSites$geneMark <- ifelse(
    toupper(standardizedDereplicatedSites$inGene)!="FALSE",
    paste0(standardizedDereplicatedSites$geneMark, "*"),
    standardizedDereplicatedSites$geneMark)

## attach gene marks
standardizedDereplicatedSites$nearest_refSeq_gene <- paste0(
    standardizedDereplicatedSites$nearest_refSeq_gene,
    standardizedDereplicatedSites$geneMark)

#===================GENERATE EXPANDED CLONE DATAFRAMES======================
#barplots
abundCutoff.barplots <- getAbundanceThreshold(standardizedDereplicatedSites, 10)

barplotAbunds <- getAbundanceSums(filterLowAbund(standardizedDereplicatedSites,
                                                 abundCutoff.barplots),
                                  c("CellType", "Timepoint"))

barplotAbunds <- order_barplot(barplotAbunds)
CellType_order <- unique(barplotAbunds$CellType)
barplotAbunds$CellType <- factor(barplotAbunds$CellType, 
                                 levels=CellType_order)

#detailed abundance plot
abundCutoff.detailed <- getAbundanceThreshold(standardizedDereplicatedSites, 50)

detailedAbunds <- getAbundanceSums(filterLowAbund(standardizedDereplicatedSites,
                                                 abundCutoff.detailed),
                                  c("CellType", "Timepoint"))

categorySums <- sapply(split(detailedAbunds$estAbundProp,
                             detailedAbunds$maskedRefGeneName),sum)

detailedAbunds$maskedRefGeneName <- factor(detailedAbunds$maskedRefGeneName,
                                           levels=names(sort(categorySums)))
detailedAbunds$CellType <- factor(detailedAbunds$CellType,
                                  levels=CellType_order)

#================Longitudinal Behaviour===============================
longitudinal <- as.data.frame(standardizedDereplicatedSites)[,c("Timepoint",
                                                                "CellType",
                                                                "estAbundProp",
                                                                "posid")]
longitudinal$CellType <- factor(longitudinal$CellType,
                                levels=CellType_order)

has_longitudinal_data <- length(unique(longitudinal$Timepoint)) > 1


#==================DETAILED REPORTS FOR BAD ACTORS=====================
badActors <- c("LMO2", "IKZF1", "CCND2", "HMGA2", "MECOM")

badActorData <- sapply(badActors, function(badActor){
  hasBadActor <- grepl(badActor, standardizedDereplicatedSites$X5pNrstOnco)
  badActorWithin100K <- abs(standardizedDereplicatedSites$X5pNrstOncoDist) <= 100000
  standardizedDereplicatedSites[hasBadActor & badActorWithin100K]
})

badActorData <- lapply(badActorData, function(x){
  x$CellType <- factor(x$CellType, levels=CellType_order)
  x
})

#==================SET VARIABLES FOR MARKDOWN REPORT=====================
timepoint <- levels(sets$Timepoint)

cols <- c("Trial", "GTSP", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN")
summaryTable <- arrange(sets,Timepoint,CellType)
summaryTable <- summaryTable[,cols]

##cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
##          "Replicates", "FragMethod", "VCN", "S.chao1", "Gini", "Shannon")
cols <- c("Patient", "Timepoint", "CellType", "UniqueSites",
          "Replicates", "FragMethod", "VCN", "Gini", "Shannon")
popSummaryTable <- merge(sets,  populationInfo, by.x="GTSP", by.y="group")
popSummaryTable <- arrange(popSummaryTable,Timepoint,CellType)
##popSummaryTable <- popSummaryTable[,cols]

cols <- c("Trial", "GTSP", "Replicates", "Patient", "Timepoint", "CellType", 
          "TotalReads", "UniqueSites", "FragMethod", "VCN", "Gini", "Shannon")
summaryTable <- popSummaryTable[,cols]

summaryTable$VCN <- ifelse(summaryTable$VCN == 0, NA, summaryTable$VCN)
    
timepointPopulationInfo <- melt(timepointPopulationInfo, "group")

#==================Get abundance for multihit events=====================
message("Fetching multihit sites and estimating abundance")
dbConn <- dbConnect(MySQL(), group=db_group_sites)
info <- dbGetInfo(dbConn)
dbConn <- src_sql("mysql", dbConn, info = info)
sites.multi <- merge( suppressWarnings(getMultihitLengths(sampleName_GTSP, dbConn)), sampleName_GTSP)
if( nrow(sites.multi) > 0 ) {
    sites.multi <- sites.multi %>%
    group_by(multihitID) %>%
    mutate(replicate=as.integer(as.factor(sampleName)))
    
    dfr <- data.frame(ID=sites.multi$multihitID,
                      fragLength=sites.multi$length,
                      replicate=sites.multi$replicate)
    
    if(use.sonicLength){
      estAbund.uniqueFragLen <- function(location, fragLen, replicate=NULL){
        if(is.null(replicate)){replicate <- 1}  #Need for downstream workflow
        dfr <- data.frame(location = location, fragLen = fragLen, 
                          replicate = replicate)
        dfr_dist <- distinct(dfr)
        site_list <- split(dfr_dist, dfr_dist$location)
        theta <- sapply(site_list, function(x){nrow(x)})
        theta <- theta[unique(dfr$location)]
        list(theta=theta)
      }
      estAbund <- estAbund.uniqueFragLen
    }
    
    if(length(unique(dfr$replicate))==1){
        estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength)
    }else{
        estimatedAbundances <- estAbund(dfr$ID, dfr$fragLength, dfr$replicate)
    }
    
    sites.multi <- subset(sites.multi, !duplicated(multihitID))
    sites.multi$estAbund <- round(estimatedAbundances$theta[as.character(sites.multi$multihitID)])
    
    sites.multi <- merge(sites.multi, sets, by="GTSP")
    sites.multi <- sites.multi %>% group_by(Patient, Timepoint, CellType) %>%
    mutate(Rank=rank(-estAbund, ties.method="max"))

}

save.image(RDataFile)
##end setting variables for markdown report

fig.path <- paste(unique(trial), unique(patient),
                  format(Sys.Date(), format="%Y%m%d"), "Figures",
                  sep=".")

#### begin generating markdown ####
unlink(fig.path, force=TRUE, recursive=TRUE)
mdfile <- paste(unique(trial), unique(patient),
                format(Sys.Date(), format="%Y%m%d"), "md",
                sep=".")

htmlfile <- gsub("\\.md$",".html",mdfile)
options(knitr.table.format='html')
theme_set(theme_bw()) #for ggplot2
knit(file.path(codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(codeDir, "GTSPreport.css") )


message("\nReport ", htmlfile, " is generated from ", csvfile)
save.image(RDataFile)


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
## intSiteRetriever package was installed from github as follows
## git clone https://github.com/BushmanLab/intSiteRetriever.git
## cd intSiteRetriever
## git checkout remove_multihitClusterID
## R
## devtools::document()
## devtools::install()

source(file.path(codeDir, "utilities.R"))
source(file.path(codeDir, "specimen_management.R"))
source(file.path(codeDir, "estimatedAbundance.R"))
source(file.path(codeDir, "dereplicateSites.R"))
source(file.path(codeDir, "standardizeSites.R"))
source(file.path(codeDir, "read_site_totals.R"))
source(file.path(codeDir, "populationInfo.R"))
source(file.path(codeDir, "abundanceFilteringUtils.R"))

load("debug.RData")

knit(file.path(codeDir, "GTSPreport.Rmd"), output=mdfile)
markdownToHTML(mdfile, htmlfile, extensions=c('tables'),
               options=c(markdownHTMLOptions(defaults=T), "toc"),
               stylesheet=file.path(codeDir, "GTSPreport.css") )
cmd <- paste("Rscript ~/geneTherapyPatientReportMaker/printReportToPdf.R", htmlfile)
system(cmd)
pdffile <- sub("html$", "pdf", htmlfile)

system(paste("scp", htmlfile, "microb215:Sites/share/GTSPReports/bushmanlab/"))
system(paste("scp", pdffile, "microb215:Sites/share/GTSPReports/bushmanlab/"))
              

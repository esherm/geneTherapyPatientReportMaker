# geneTherapyPatientReportMaker
For a specific patient, report integration near oncogenes and potentially expanded clones across cell types and multiple time points
	
# Input
A csv file such as `sampleName_GTSP.csv` to describe the replicates and the samples.
```
head sampleName_GTSP.csv

sampleName,GTSP
GTSP0308-1,GTSP0308
GTSP0308-2,GTSP0308
GTSP0308-3,GTSP0308
GTSP0308-4,GTSP0308
GTSP0309-1,GTSP0309
GTSP0309-2,GTSP0309
GTSP0309-3,GTSP0309
GTSP0309-4,GTSP0309
```

* only `sampleName`, `GTSP` are necessary,
* there could be more columns, but the extra columns are ignored,
* the GTSPxxxx names must correspond to the same patient,
* the GTSPxxxx names must be in the `specimen_management` database,
* all sites should be computed based on one reference genome.
  
# Output
An embeded html file named `$trial.$patient.$today.html`

# Code example
```
Rscript path/to/makeGeneTherapyPatientReport.R  #assuming sampleName_GTSP.csv in current folder
Rscript path/to/makeGeneTherapyPatientReport.R sampleName_GTSP.csv
Rscript path/to/makeGeneTherapyPatientReport.R path/to/test.csv
```
	
# Note
Do NOT run multiple instances with in the same folder

# Database connfig file location

config file should be in home directory and called .my.cnf,
e.g. ~/.my.cnf

The .my.cnf format is as follows:

```
[intSitesDEV-dev]
user=YYYYYYY
password=XXXXXX
host=microb98.med.upenn.edu
port=3309
database=intsitesdev
```

# Dependencies

intSiteRetriever : https://github.com/esherm/intSiteRetriever 
(at present get the project in current folder:
```
git clone https://github.com/esherm/intSiteRetriever
```)
hiAnnotator
reldist

Cancer Gene list:

```
git clone https://github.com/anatolydryga/CancerGeneList.git
```

# Testing

Run in the R console:

```bash
library(testthat)
test_dir(".")
```

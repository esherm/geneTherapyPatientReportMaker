#### geneTherapyPatientReportMaker
For a specific patient, report integration near oncogenes and potentially expanded clones across cell types and multiple time points
	
#### Input
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

#or 
Rscript path/to/check_patient_GTSP.R pFR03

sampleName,GTSP,patient
GTSP0308-1,GTSP0308,pFR03
GTSP0308-2,GTSP0308,pFR03
GTSP0308-3,GTSP0308,pFR03
```

* only `sampleName`, `GTSP` columns are necessary and the rest are ignored,
* the `GTSPxxxx` names must correspond to the same patient,
* the `GTSPxxxx` names must be in the `specimen_management` database,
* all sites should be computed based on one reference genome.
  
#### Output
`$trial.$patient.$today.html`

#### Code example
- 1. `check_patient_GTSP.R`: get available datasets, generate a csv file for a patient 
```
 Rscript path/to/check_patient_GTSP.R                    #get all processed samples
 Rscript path/to/check_patient_GTSP.R pFR03              #get data sets for patient pFR03 and output to csv format
 Rscript path/to/check_patient_GTSP.R pFR03 > pFR03.csv  #get data sets for patient pFR03 and output to tmp.csv
```

- 2. `makeGeneTherapyPatientReport.R`: generate report for a patient from the csv file 
```
Rscript makeGeneTherapyPatientReport.R                     #read in sampleName_GTSP.csv by default
Rscript path/to/makeGeneTherapyPatientReport.R pFR03.csv   #generated above
```

#### Note
Do NOT run multiple instances within the same folder

#### Database connfig file location

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

#### Dependencies

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

#### Testing

Run in the R console:

```bash
library(testthat)
test_dir(".")
```

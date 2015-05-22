# geneTherapyPatientReportMaker
For a specific patient, report integration near oncogenes and potentially expanded clones across cell types and multiple time points

# Input

Report is done for individual patient for multiple timepoints and cell types, each corresponding to unique
GTSP. For each sample multiple preps can be used. 

Sample name corresponds to name used in integration site calling DB.

Format is (tab-separated):
```{bash}
GTSP    sampleName
```

# Checked Assumptions
* The patient report is done for a single patient.
* All sites should be computed based on one reference genome.
* All sample names and GTSPs should exist in intSite and specimen managemnet databases.


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

# DO_pipe
Cluster-based QTL analysis pipeline built around r/qtl2

Repository for keeping track of the DO_pipe code. 

The pipeline has "5" arms:

This pipeline involves "local" processing in RStudio and "remote" processing on a cluster. We use the FARM research cluster, which is an Ubuntu OS running a SLURM scheduler. 

## Pre-requisites:
- RStudio and intermediate R skills
- Adequate local storage (> 10 Gb)
- Access to FARM cluster and UNIX command line basics
- Cyberduck or knowledge of scp

## Getting Started:

Refer to this README alongside the documentation in the DO_PIPE R markdown. Subheaders of the README refer to specific code chunks in the markdown - in RStudio you can navigate between chunks using the drop down menu in the bottom left of the markdown file. Sometimes running the chunk is all that is required, but sometimes the user is prompted for input in the RStudio console. After running each chunk, the outputs should be inspected for completeness and correctness. 

Begin by creating a project folder and placing the DO_PIPE.rmd file inside. Follow the directions below. 

## Project Arms: 

### 0. Project Setup

* {r Set up Project Environment}
** Run chunk 
** In console: Designate study prefix, e.g. DOWL, DO2, DOmoms
** Inspect new project folder architecture
** Add input files and format master_phenotypes.csv according to the instructions in the markdown

### 1. Phenotype Diagnostic (local)

* {r Pheno Process}
** Run chunk 
** In console: Change subject IDs in master_pheno is desired
** Inspect pheno.csv and covar.csv in ./inputs-processed 

* {r Pheno Diag}
** Run chunk
** Inspect diagnostics plots in ./diagnostic/phenotypes
** Adjust pdf dimensions, margins, and font sizes as needed and re-run chunk

### 2. Generate Genotype Files (local)

* {r Process Final Report}
** Run chunk
** In console: Choose to append study prefix 
** Confirm (but do not open) processed_finrep.txt in ./inputs-processed

* {r Generate Genotype Files}
** Run chunk 
** Inspect genotype files in ./inputs-processed/genotypes

* {r Create Control File}
** Run chunk
** Inspect control.json file in ./inputs-processed

* {r Generate Cross2}
** Run Chunk
** Some individuals may be excluded for missing data (this is normal) 
** Inspect the cross-summary.txt in ./diagnostic

* {r Generate Covariate Object}
** Run Chunk
** In console: Select covariates to for genome scan
** In console: Enter 'covar' to inspect design matrix

### 3. Genotype Diagnostic (cluster)
### 4. Permutation Thresholds (cluster)
### 5. Genome Scan and Pathway Enrichment Analysis


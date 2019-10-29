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

*{r Setup Project Environment}*
- Run chunk 
- In console: Designate study prefix, e.g. DOWL, DO2, DOmoms
- Inspect new project folder architecture
- Add input files and format master_phenotypes.csv according to the instructions in the markdown

### 1. Phenotype Diagnostic (local)

*{r Pheno Process}*
- Run chunk 
- In console: Change subject IDs in master_pheno is desired
- Inspect pheno.csv and covar.csv in *./inputs-processed* 

*{r Pheno Diag}*
- Run chunk
- Inspect diagnostics plots in *./diagnostic/phenotypes*
- Adjust pdf dimensions, margins, and font sizes as needed and re-run chunk

### 2. Generate Genotype Files (local)

*{r Process Final Report}*
- Run chunk
- In console: Choose to append study prefix 
- Confirm (but do not open) processed_finrep.txt in *./inputs-processed*

*{r Generate Genotype Files}*
- Run chunk 
- Inspect genotype files in *./inputs-processed/genotypes*

*{r Create Control File}*
- Run chunk
- Inspect control.json file in *./inputs-processed*

*{r Generate Cross2}*
- Run Chunk
- Some individuals may be excluded for missing data (this is normal) 
- Inspect the cross-summary.txt in *./diagnostic*

*{r Generate Covariate Object}*
- Run Chunk
- In console: Select covariates to for genome scan
- In console: Enter 'covar' to inspect design matrix

### 3. Genotype Diagnostic (cluster)

This next section takes place entirely on the FARM and involves both **R** and **UNIX**. 

*Setup Project Environment and Generate GPR*

1. Create a project directory
2. From the custom folder, upload genScanInputs.sh into the project directory
3. Make the script executable and run it, filling in your study prefix and email where noted

> chmod 755 ./genScanInputs\
> ./genScanInputs -p PREFIX -e EMAIL\
> cd scripts\
> sbatch runBatch.sh

4. FARM will email you when job starts and finishes (approx. 20-60 minutes)
*While you wait*

We will be using interactive R script to run certain aspects of the analysis, this is how to open an R session on FARM:

>module load R\
>srun --mem=60000 --time=10:00:00 --pty R

To run the diagnostic, your personal R library will need the following packages installed:

- install.packages(“fst”)
- install.packages(“qtl2”)
- install.packages(“dplyr”)
- install.packages(“pryr”)
- install.packages(“tibble”) 
- install.packages(“ggpubr”)
- install.packages(“ggplot2”)
- install.packages(“ggrepel”)
- install.packages(“broman”)
- install.packages(“qtlcharts”)
- install.packages(“gridExtra”)
- install.packages(“data.table”)

To quit interactive R session, do not save the session:

> q()


5. Confirm creation of gprRaw.rds in *./outputs*

*Run Diagnostic Pipeline*

1. In the project directory ./inputs folder, add the following files from your local directory:
- prefix.rds in local *./inputs*
- allelecodes.csv in local *./inputs/MUGA*
- processed_finrep.txt in local *./inputs-processed*
- covar.rds in local *./scan-inputs*



### 4. Permutation Thresholds (cluster)
### 5. Genome Scan and Pathway Enrichment Analysis


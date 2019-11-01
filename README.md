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

### Custom Scripts
DO_pipe will call a number of custom scripts available at: https://github.com/exsquire/DO_pipe

Use the green "Clone or download" button on the righthand side to clone or download the contents of the git repository.
Place the unzipped contents in the project directory. Before you run code, your project folder should look like this:

- /custom/
- DO_PIPE.rmd
- README.md

## Project Arms: 

### 0. Project Setup

*{r libraries and packages}
- Run chunk
- If you do not have one or more packages installed, install them now. 

*{r Code Source}
- Run chunk
- If you receive an error, check that the "custom" folder is in the first level of your project folder

*{r Setup Project Environment}*
- Run chunk 
- In console: Designate study prefix, e.g. DOWL, DO2, DOmoms
- Inspect new project folder architecture
- Add input files and format master_phenotypes.csv according to the instructions in the markdown

MUGA folders can be found on the qtl2 webstite under 'Founder Genotypes and SNP Maps': https://kbroman.org/qtl2/pages/prep_do_data.html  
Mouse_genes.sqlite database file can be found here:   
https://doi.org/10.6084/m9.figshare.5280238.v6  

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
- Ensure "sex_codes" argument matches sex designations in covariate file
- e.g. sex_codes=c("f" = "Female"), if females are designated with "f"
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

This next section takes place entirely on the FARM and involves both **R** and **UNIX**. All references to the "project directory" in this section refer to the one on the cluster. 

#### *Setup Project Environment and Generate GPR*

1. Create a project directory
2. From the *./custom*, upload genScanInputs.sh into the project directory
3. From *./scan-inputs*, upload cross2 object into project directory
4. Make the script executable and run it, filling in your study prefix and email where noted - your cluster folder structure will be made for you.

> chmod 755 ./genScanInputs.sh\
> ./genScanInputs.sh -p PREFIX -e EMAIL\
> cd scripts\
> sbatch runBatch.sh

5. FARM will email you when job starts and finishes (approx. 20-60 minutes)
#### *While you wait*

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


6. Confirm creation of gprRaw.rds in *./outputs*

#### *Run Diagnostic Pipeline*

1. Within *./inputs*, upload the following files from your local directory:
- prefix.rds in local *./inputs*
- allelecodes.csv in local *./inputs/MUGA*
- processed_finrep.txt in local *./inputs-processed*
- covar.rds in local *./scan-inputs*

2. Within the project directory, make *./cleanotype* subdirectory
3. Within *./cleanotype*, make *./output* and *./scripts* subdirectories (do not misspell)
4. Upload the following files from local *./custom* to *./scripts
- sexDiag.R
- cleanCross.R
- computeScanInputs.R
- cleanotype.R
5. From *./cleanotype/scripts*, enter the following commands for an interactive R session:

> module load R\
> srun --mem=60000 --time=10:00:00 --partition=high --pty R

6. From Interactive R

> source("cleanotype.R")

Allow program to run to completion.

7. Inspect the contents of *./cleanotype/output*
- Check Perc_Miss_Plot for missing genotyping data per sample
- Check dupeDiag plot for possible duplication errors
- Check sexDiag for sex misclassification errors or XO females
- See: https://kbroman.org/qtl2/assets/vignettes/do_diagnostics.html for information on the other plots
- *./cleanotype/output/objects* holds the intermediate data files for plot generation 

8. Inspect the contents of *./outputs* in project folder, make sure you have the following files:
- apr_Clean.rds
- cross_Clean.rds
- kLOCOC_Clean.rds

9. Quit Interactive R - do not save

> q()

### 4. Permutation Thresholds (cluster)

1. Within your project folder, create a ./permutations subdirectory and upload the following custom files to it:
- theDirector.R
- DO_PIPE_scryptic.sh
- permUtil.R

theDirector pulls apr, cross, and kLOCO from *../outputs* and covar and prefix from *../inputs/* then generates controller.rds file and batchArg.txt, allowing DO_PIPE_scryptic to request adequate cluster resources for splitting the permutations among an array jobs. 

DO_PIPE_scryptic sets up the folder architecture, R code, and sbatch script for the permutation array run.

Array jobs create many output files, permUtil sticks them together.

2. From within *./permutations*, run the following on the command line to open interactive R and call theDirector

> module load R\
> srun --mem=60000 --time=10:00:00 --partition=high --pty R\

3. From Interactive R

> source("theDirector")

Note the estimated time to completion. 

> q()

Don't save. 

4. From command line, run the following to begin permutation analysis

> chmod 755 DO_PIPE_scryptic.sh\
> ./DO_PIPE_scryptic.sh\
> cd scripts\
> sbatch run.sh #For now, personal email functionality is limited.\
> squeue -u USERNAME #to check your runs

After the first round of permutations (~10-30 minutes) download one of the files from *./permutations/outputs* and inspect in RStudio.
Allow runs to complete.

5. When job is complete, run the following commands from *./permutations/scripts*

> module load R\
> srun --mem=60000 --time=05:00:00 --partition=high --pty R

6. From Interactive R

> source("permUtil.R")

permUtil will place the following files in your project *./outputs* folder:
- **fullPerm.rds**: The full permutation matrix of your run. Can be used to calculate individual significant thresholds at desired quantiles
- **permThresh.rds**: Infividual permutations thresholds at 0.63, 0.9, and 0.95 quantiles. 

If permUtil alerts you to failed runs, re-run *run.sh* from *./permutations/scripts*, submitting the failed arrayIDs
ex.) 
> sbatch --array=2,10,1337,1338 ./run.sh

See: https://slurm.schedmd.com/job_array.html for more details on array jobs.

7. Download the following files from the project *./outputs* directory to your local *./scan-inputs* directory
- apr_Clean.rds
- cross_Clean.rds
- kLOCO_Clean.rds
- fullPerm.rds
- permThresh.rds

### 5. Genome Scan and Pathway Enrichment Analysis

*{r Run Scan}*
- Run chunk

*{r Find QTL}* 
- Modify probs vector if desired. Default is c(0.63, 0.9, and 0.95) 
- Run Chunk 
- Inspect peak files in *./outputs/data* at different levels of significance

*{r Visualize all peaks}*
- Select path to desired peak file 
- We do not plot peaks on the X chromosome by default
- Plots can be inspected in *./outputs/plots/lod*
- To run other peaks, move plots to separate folder to prevent overwritting, change *peaks* and re-run chunk
- peakSummary.png shows global architecture of QTL

*{r Find Candidate Genes}*
For finding candidate genes within the QTL, we use the full set of Mouse gene annotations made available at:
https://kbroman.org/qtl2/assets/vignettes/user_guide.html#snp_association
Direct Link: 
https://doi.org/10.6084/m9.figshare.5280238.v6

We query the database, looking for all non-predicted genes within the confidence interval of the peaks.
- Run Chunk
- Candidate genes stored in *./outputs/data/peak-genes*

*{r Enrichr Query}*

A structural variant of any one of the genes within the QTL could be responsible for the variance in phenotypic expression. It could also be multiple genes working together. 

The final step in our QTL pipeline is to check if any of the genes within our QTL are involved in (aka "enriched for") any validated gene sets. This can link genes to genomes, biological pathways, metabolites, drugs, diseases, microbes, and chemical compunds. 

We access the Ma'ayan Laboratory's tool, Enrichr, using methods from the enrichR package (https://cran.r-project.org/web/packages/enrichR/index.html), which queries multiple gene databases looking for gene set enrichment. 

- Defaults to searching all Enrichr libraries
- Can modify dbs to search specific libraries or add specific libraries to the dbNonGrata vector to exclude libraries


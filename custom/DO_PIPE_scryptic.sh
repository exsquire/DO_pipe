#!/bin/bash

#echo Creating Run Directory Architecture
#Create folders if needed
if [ ! -d outputs ];
then
    mkdir outputs
fi

if [ ! -d config ];
then 
    mkdir config
fi

if [ ! -d inputs ];
then
    mkdir inputs
fi

if [ ! -d scripts ];
then
    mkdir scripts
fi
 
#Define variables to pass to batch script
LOC=$(pwd) #Location of scryptic.sh

#A single chunk is determined by the permAlloc script
PERMS=`cat batchArgs.txt | cut -d " " -f 1`
ARRAY=`cat batchArgs.txt | cut -d " " -f 2`
TIME=`cat batchArgs.txt | cut -d " " -f 3`
CORES=`cat batchArgs.txt | cut -d " " -f 4`

if [[ -z $ARRAY || -z $TIME || -z $CORES || -z $PERMS ]]; then
  echo 'one or more permutation arguments are undefined'
  exit 1
fi


BATCHS="./run.sh"

#--------------------------------------------------
#Define Input Name Variables
APR=$(ls ../outputs | grep -i apr_Clean)
CROSS=$(ls ../outputs | grep -i cross_Clean)
KLOCO=$(ls ../outputs | grep -i kLOCO_Clean)
COVAR=$(ls ../inputs | grep -i covar)
PREFIX=$(ls ../inputs | grep -i prefix)
CONTROL=$(ls ./ | grep -i control)

if [[ -z $APR || -z $CROSS || -z $KLOCO || -z $COVAR || -z $PREFIX || -z $CONTROL ]]; then
  echo 'one or more variables are undefined'
  exit 1
fi
#--------------------------------------------------
RCODE="./Rcode.R"

cat <<END > $RCODE
#Test Scryptic
library(qtl2)

#Read in inputs
apr <- readRDS("$APR")
cross <- readRDS("$CROSS")
kLOCO <- readRDS("$KLOCO")
covar <- readRDS("$COVAR")
prefix <- readRDS("$PREFIX")
ctrl <- readRDS("$CONTROL")

#Intialize array id - array id will only function as a sequential designation
args<-as.integer(unlist(strsplit(commandArgs(TRUE)," ")))
print(args)

arrayid <-args[1]
print(arrayid)

#Cores controlled from command line
#No longer a need for control file.

start <- ctrl[arrayid, 2]
stop  <- ctrl[arrayid, 3]

#Run scan1perm 
perm <- scan1perm(apr,
                  $(echo 'cross$pheno[, start:stop, drop = FALSE]'),
                  kinship = kLOCO,
                  addcovar = covar,
                  cores = 4, 
                  n_perm = $PERMS)

out <- data.frame(perm, check.names = F)

#can also save "out" data
saveRDS(out, file = paste0(prefix,"_perm_",arrayid,".rds"))
END


#---------------------------------------------------
#Make R code executable
chmod 755 $RCODE
#---------------------------------------------------
#sbatch parameters modified for Farm environment - will not work on Ceres
cat <<END > $BATCHS
#!/bin/bash -l

#SBATCH -J permArray
#SBATCH -N 1
#SBATCH -c $CORES
#SBATCH --mem-per-cpu=2500
#SBATCH --array=1-$ARRAY
#SBATCH --partition=high
#SBATCH --time=$TIME

#Email me here when job starts, ends, or sh*ts the bed
#SBATCH --mail-user=excel.que@gmail.com
#SBATCH --mail-type=ALL

#Run script out of 'scripts' folder
#Requires a project/{input,output,scripts,config} folder structure
#SBATCH -o $LOC/config/permArray-%A_%a.out
#SBATCH -e $LOC/config/permArray-%A_%a.err

#scratch designation
export SCR_DIR=$(echo '/scratch/$USER/$SLURM_JOBID/$SLURM_ARRAY_TASK_ID')

#Use Full Paths
export OUTPUT_DIR=$LOC/outputs

# Load R
module load R

# Create scratch & copy everything over to scratch
mkdir -p $(echo '$SCR_DIR')
cd $(echo '$SCR_DIR')

#Copy over everything for permutation Run
cp -p $LOC/../outputs/$APR .
cp -p $LOC/../outputs/$CROSS .
cp -p $LOC/../outputs/$KLOCO .
cp -p $LOC/../inputs/$COVAR .
cp -p $LOC/../inputs/$PREFIX .
cp -p $LOC/inputs/$CONTROL . 
cp -p $LOC/scripts/Rcode.R .

#Confirm presence of input files in scratch
echo "before srun in dir"
pwd
echo "contents"
ls -al

# SciNet Ceres automatically cleans up scratch after jobs
# This is nice! But not ubiquitous.
# Termination Signal Trap - when a job goes over its walltime or user cancels job
termTrap()
{
        echo "Termination signal sent. Clearing scratch before exiting"
        # do whatever cleanup you want here
        rm -dr $(echo '$SCR_DIR')
        exit -1
}
# associate the function "term_handler" with the TERM signal
trap 'termTrap' TERM

#Run lightweight R instance
srun R --vanilla --args "$(echo '$SLURM_ARRAY_TASK_ID')" <  ./Rcode.R

#Confirm that output made it
echo "after srun, directory"
ls -al
echo work=$(echo '$WORK_DIR')
echo scr=$(echo '$SCR_DIR')

# Copy results over
cd $(echo '$OUTPUT_DIR')
#change to output directory (now the pwd)
cp -p $(echo '$SCR_DIR')/*_perm_* .

#Routine Scratch Cleanup
rm -dr $(echo '$SCR_DIR')

echo "End of program at $(echo '`date`')"
END

#------------------------------------------------------------
chmod 755 $BATCHS

#Cleanup
mv batchArgs.txt ./inputs
mv *controller.rds ./inputs
mv Rcode.R ./scripts
mv *scryptic.sh ./scripts
mv run.sh ./scripts
mv theDirector.R ./scripts

echo Done

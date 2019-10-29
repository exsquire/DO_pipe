#!/bin/bash

#README: genScanInputs.sh is an executable that takes in the cross2 object
#and generates the compute-heavy genome scan input: gpr.
#This code should be run in the directory with the cross2 object.
#The cross2 object should be the only file in the folder with the word
#"cross" in the filename. 

module load R

#Load and check for Cross2
CROSS=$(ls | grep -i cross)

if [[ -z $CROSS ]]; then
  echo 'Abort: Cross2 variable undefined'
  exit 1
fi
#-----------------------------------------------
echo Creating Run Directory Architecture
#Create folders if needed
if [ ! -d inputs ];
then
    mkdir inputs
fi

if [ ! -d outputs ];
then
    mkdir outputs
fi

if [ ! -d config ];
then 
    mkdir config
fi

if [ ! -d scripts ];
then
    mkdir scripts
fi 
#-------------------------------------------------
#User defined variables to pass to batch script

LOC=$(pwd)

while getopts p:e: option
do
  case "${option}"
  in
  p) PREFIX=${OPTARG};;
  e) EMAIL=${OPTARG};;
  esac
done


if [ -z $PREFIX ]; then
  PREFIX=""    #Default prefix name
fi

if [ -z $EMAIL ]; then
  EMAIL="excel.que@gmail.com"
fi

#--------------------------------------------------
RCODE="./genScanInputs.R"

cat <<END > $RCODE
#Generate gpr on cluster
library(qtl2)

#Read in cross
cross <- readRDS("$CROSS")

#gpr, apr, kLOCO
gpr <- calc_genoprob(cross,$(echo 'cross$gmap'), cores = 10)

#Save to RDS
saveRDS(gpr, file = paste0("$PREFIX","_gprRaw.rds"))

END
#---------------------------------------------------
#Take any files with .rds or .R and move them into the inputs folder
echo Moving .rds and .R files to input folders
mv *.rds inputs/
chmod 755 $RCODE
mv *.R inputs/
#---------------------------------------------------
BATCHS="./runBatch.sh"
#sbatch parameters modified for Farm environment - will not work on Ceres
cat <<END > $BATCHS
#!/bin/bash -l

#SBATCH -J genScanInputs
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --partition=high
#SBATCH --time=05:00:00

#Email me here when job starts, ends, or sh*ts the bed
#SBATCH --mail-user=$EMAIL
#SBATCH --mail-type=ALL

#Run script out of 'scripts' folder
#Requires a project/{input,output,scripts,config} folder structure
#SBATCH -o $LOC/config/genScanInputs-%A_%a.out
#SBATCH -e $LOC/config/genScanInputs-%A_%a.err

#Farm scratch designation
export SCR_DIR=$(echo '/scratch/$USER/$SLURM_JOBID')

#Use Full Paths
export WORK_DIR=$LOC/inputs
export OUTPUT_DIR=$LOC/outputs

# Load R
module load R

# Create scratch & copy everything over to scratch
mkdir -p $(echo '$SCR_DIR')
cd $(echo '$SCR_DIR')

#Copy over everything for permutation Run
cp -p $(echo '$WORK_DIR')/$CROSS .
cp -p $(echo '$WORK_DIR')/genScanInputs.R .

#Confirm presence of input files in scratch
echo "before srun in dir"
pwd
echo "contents"
ls -al


# Termination Signal Trap - when a job goes over its walltime or user cancels job
termTrap()
{
        echo "Termination signal sent. Clearing scratch before exiting"
        # do whatever cleanup you want here
        rm -rf $(echo '$SCR_DIR')/*
        exit -1
}
# associate the function "term_handler" with the TERM signal
trap 'termTrap' TERM

#Run lightweight R instance
srun R --vanilla <  ./genScanInputs.R

#Confirm that output made it
echo "after srun, directory"
ls -al
echo work=$(echo '$WORK_DIR')
echo scr=$(echo '$SCR_DIR')

# Copy results over
cd $(echo '$OUTPUT_DIR')
#change to output directory (now the pwd)
cp -p $(echo '$SCR_DIR')/*.rds .

#Routine Scratch Cleanup
rm -rf $(echo '$SCR_DIR')/*

#print out scratch contents to confirm removal
ls -al $(echo '$SCR_DIR')/*

echo "End of program at $(echo '`date`')"
END

#------------------------------------------------------------
chmod 755 $BATCHS
mv $BATCHS ./scripts
echo Done.



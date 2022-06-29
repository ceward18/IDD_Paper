#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=IDD_Known
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --mem 100G
#SBATCH --partition=cpu2019,cpu2021
#SBATCH --array=916,927,924,1700,1790

####### Set environment variables ###############
module load gcc/7.3.0
module load lib/openblas/0.3.5-gnu
module load java/9.0.4
module load lib/hdf5/1.10.0.1
export PATH=/home/caitlin.ward/R/bin:$PATH
export R_LIBS=/home/caitlin.ward/R/lib64:$R_LIBS

####### Run your script #########################
Rscript model_fitting_estimatedE.R $SLURM_ARRAY_TASK_ID

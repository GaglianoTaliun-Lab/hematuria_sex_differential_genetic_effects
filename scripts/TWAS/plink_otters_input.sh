#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --array=1-19
#SBATCH --time=00:30:00
#SBATCH --job-name=plink
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=6G

module load StdEnv/2020
module load gcc
module load python/3.9
module load tabix
module load r/4.3.1
module load plink/1.9b_6.21-x86_64 r/4

export R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/hematuria_sexspecific"

chr=${SLURM_ARRAY_TASK_ID}

# ######### Extract chromosome info from plink bfiles ############
plink \
    --bfile ${project_dir}/reference_data/g1000_eur \
    --chr ${chr} \
    --make-bed \
    --maf 0.01 \
    --out ${project_dir}/reference_data/g1000_eur_chr${chr}
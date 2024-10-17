#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=4:00:00
#SBATCH --array=13
#SBATCH --job-name=bolt_reml_M
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=8

module load StdEnv/2020
module load gcc/9.3.0 
module load bolt-lmm/2.4

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
chr=${SLURM_ARRAY_TASK_ID}
sex="male"

bolt \
    --bfile=${project_dir}/data_for_bolt/ukb_chr${chr}_${sex}_imputed_pruned_1Mb \
    --phenoFile=${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
    --phenoCol=X593 \
    --covarFile=${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
    --qCovarCol=f.34.0.0 \
    --qCovarCol=f.22009.0.{1:10} \
    --bgenMinINFO=0.4 \
    --bgenMinMAF=1e-3 \
    --reml \
    --remlNoRefine \
    --numThreads=8
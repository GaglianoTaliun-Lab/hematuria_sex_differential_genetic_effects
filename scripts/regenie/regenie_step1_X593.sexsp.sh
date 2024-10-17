#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH --job-name=regenie_step1_hem
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=6G

module load StdEnv/2020
module load gcc/9.3.0
module load regenie/3.2.1

data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

group=${SLURM_ARRAY_TASK_ID}

for sex in "female" "male"
do
  regenie \
    --step 1 \
    --bed ${data_dir}/plink_data/ukb_chrs1-23_QCd_hg38 \
    --keep ${project_dir}/data_for_regenie/ukb_WB_${sex}_ageG${group}_ids.txt \
    --phenoFile ${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
    --phenoCol X593 \
    --covarFile ${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
    --covarCol f.34.0.0 \
    --covarCol f.22009.0.{1:10} \
    --bt \
    --bsize 1000 \
    --lowmem \
    --strict \
    --lowmem-prefix tmp_rg_${sex}_G${group} \
    --out ${project_dir}/results_regenie/step1_ukb_X593_${sex}_G${group}
done

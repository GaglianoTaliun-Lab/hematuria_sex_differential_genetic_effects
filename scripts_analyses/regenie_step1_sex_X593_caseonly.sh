#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=5:00:00
#SBATCH --job-name=regenie_step1_X593
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=6G

module load StdEnv/2020
module load gcc/9.3.0
module load regenie/3.2.1

data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

regenie \
  --step 1 \
  --bed ${data_dir}/plink_data/ukb_chrs1-23_QCd_hg38 \
  --keep ${project_dir}/data_for_regenie/sexcombined_cases.txt \
  --phenoFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
  --phenoCol f.22001.0.0 \
  --covarFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
  --covarCol f.34.0.0 \
  --covarCol f.22009.0.{1:10} \
  --bt \
  --bsize 1000 \
  --lowmem \
  --strict \
  --lowmem-prefix tmp_rg \
  --out ${project_dir}/results_regenie/step1_ukb_sex_X593casesonly

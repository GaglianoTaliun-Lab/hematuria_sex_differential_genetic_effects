#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=8:00:00
#SBATCH --array=1-22
#SBATCH --job-name=regenie_step2_egfr
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50G

module load StdEnv/2020
module load gcc/9.3.0
module load regenie/3.2.1

project_ukb_wes_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

chr=${SLURM_ARRAY_TASK_ID}
sex="male"

regenie --step 2 \
    --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen \
    --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
    --keep ${project_ukb_wes_dir}/plink_data/ukb_WB_${sex}_ids.txt \
    --phenoFile ${project_dir}/data_for_regenie/egfr_withcovariates_for_REGENIE.txt \
    --covarFile ${project_dir}/data_for_regenie/egfr_withcovariates_for_REGENIE.txt \
    --qt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --phenoCol logeGFR2009 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred ${project_dir}/results_regenie/step1_ukb_logeGFR2009_${sex}_pred.list \
    --threads 12 \
    --out ${project_dir}/results_regenie/step2_logeGFR2009_HRC_${sex}_regenie_chr${chr}
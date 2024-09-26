#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=8:00:00
#SBATCH --job-name=regenie_check_COGWI
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50G

module load StdEnv/2020
module load gcc/9.3.0
module load regenie/3.2.1

data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

chr=16

regenie --step 2 \
    --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen \
    --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
    --keep ${data_dir}/plink_data/ukb_WB_all_ids.txt \
    --extract ${project_dir}/data_for_regenie/snps_COGWI_validation_chr${chr}.txt \
    --phenoFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --covarFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --af-cc \
    --firth \
    --firth-se \
    --phenoCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred ${project_dir}/results_regenie/step1_ukb_sex_X593_all_pred.list \
    --threads 12 \
    --out ${project_dir}/results_regenie/step2_sex_X593_all_HRC_regenie_chr${chr}

chr="X"

regenie --step 2 \
    --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen \
    --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/chrX/ukb22828_cX_b0_v3_s486556.sample \
    --keep ${data_dir}/plink_data/ukb_WB_all_ids.txt \
    --extract ${project_dir}/data_for_regenie/snps_COGWI_validation_chr${chr}.txt \
    --phenoFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --covarFile ${data_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --af-cc \
    --firth \
    --firth-se \
    --phenoCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred ${project_dir}/results_regenie/step1_ukb_sex_X593_all_pred.list \
    --threads 12 \
    --out ${project_dir}/results_regenie/step2_sex_X593_all_HRC_regenie_chr${chr}
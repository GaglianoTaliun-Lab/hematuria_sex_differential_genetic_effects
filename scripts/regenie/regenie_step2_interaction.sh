#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=3:0:00
#SBATCH --array=2
#SBATCH --job-name=inter_regenie
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=50G

module load StdEnv/2020
module load gcc/9.3.0
module load regenie/3.2.1

project_ukb_wes_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

chr=${SLURM_ARRAY_TASK_ID}

regenie --step 2 \
    --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen \
    --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
    --keep ${project_ukb_wes_dir}/plink_data/ukb_WB_all_ids.txt \
    --phenoFile ${project_ukb_wes_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --covarFile ${project_ukb_wes_dir}/phenotype/case.control-wb.withcovariates-4regenie.txt \
    --extract ${project_dir}/data_for_regenie/significant_SNPs_for_interaction.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --no-condtl \
    --phenoCol X593 --covarCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred ${project_ukb_wes_dir}/regenie_out/step1_ukb_X593_pred.list \
    --threads 8 \
    --interaction f.22001.0.0 \
    --out ${project_dir}/results_regenie/power_test_X593_HRC_imputed_regenie_interaction_chr${chr}

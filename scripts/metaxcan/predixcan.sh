#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=0:30:00
#SBATCH --array=1-22
#SBATCH --job-name=spredixcan_sqtl
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G

module load StdEnv/2020
module load python/3.7.9
module load scipy-stack/2021a

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
chr=${SLURM_ARRAY_TASK_ID}

sex="males"

python ${project_dir}/MetaXcan/software/SPrediXcan.py \
	--gwas_file ${project_dir}/datasets/X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz \
	--snp_column SNP \
	--effect_allele_column effect_allele \
	--non_effect_allele_column other_allele \
	--beta_column BETA \
	--se_column SE \
	--pvalue_column P \
	--model_db_path ${project_dir}/build_sex_stratified_PrediXcan_weights/gtex_v8_imputed_b37_eur_Whole_Blood_${sex}_v2.db \
	--covariance ${project_dir}/build_sex_stratified_PrediXcan_weights/covariances/gtex_v8_Whole_Blood_chr${chr}_covariances_${sex}_age_sex.txt \
	--overwrite \
	--throw \
	--keep_non_rsid \
	--model_db_snp_key rsid \
	--output_file ${project_dir}/spredixcan_output/hematuria_${sex}_chr${chr}_gtex_v8_WholeBlood.csv


sex="females"

python ${project_dir}/MetaXcan/software/SPrediXcan.py \
	--gwas_file ${project_dir}/datasets/X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz \
	--snp_column SNP \
	--effect_allele_column effect_allele \
	--non_effect_allele_column other_allele \
	--beta_column BETA \
	--se_column SE \
	--pvalue_column P \
	--model_db_path ${project_dir}/build_sex_stratified_PrediXcan_weights/gtex_v8_imputed_b37_eur_Whole_Blood_${sex}_v2.db \
	--covariance ${project_dir}/build_sex_stratified_PrediXcan_weights/covariances/gtex_v8_Whole_Blood_chr${chr}_covariances_${sex}_age_sex.txt \
	--overwrite \
	--throw \
	--keep_non_rsid \
	--model_db_snp_key rsid \
	--output_file ${project_dir}/spredixcan_output/hematuria_${sex}_chr${chr}_gtex_v8_WholeBlood.csv
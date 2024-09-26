#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=12:00:00
#SBATCH --array=2-9
#SBATCH --job-name=regenie_step1_cond
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

condition_snpid=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/data_for_regenie/conditional_snps/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $8}')
sex=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/data_for_regenie/conditional_snps/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $2}')
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/data_for_regenie/conditional_snps/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $3}')

regenie \
  --step 1 \
  --bed ${data_dir}/plink_data/ukb_chrs1-23_QCd_hg38 \
  --keep ${data_dir}/plink_data/ukb_WB_${sex}_ids.txt \
  --phenoFile ${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
  --phenoCol X593 \
  --covarFile ${project_dir}/data_for_regenie/case.control-wb.withcovariates-4regenie.txt \
  --covarCol f.34.0.0 \
  --covarCol f.22009.0.{1:10} \
  --condition-file bed,${project_dir}/data_for_regenie/conditional_snps/genotype_conditional_${sex}_chr${chr}_${condition_snpid} \
  --condition-list ${project_dir}/data_for_regenie/conditional_snps/SNP_to_condition_${condition_snpid}.txt \
  --bt \
  --bsize 1000 \
  --lowmem \
  --strict \
  --lowmem-prefix tmp_rg_${sex}_${condition_snpid} \
  --out ${project_dir}/results_regenie/step1_ukb_X593_${sex}_cond_${condition_snpid}


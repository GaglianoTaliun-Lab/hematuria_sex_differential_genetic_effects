#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=1-22
#SBATCH --job-name=bliss
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4G

module load StdEnv/2023
module load r/4.3.1

export R_LIBS=~/.local/R/$EBVERSIONR/

project_dir="/home/fridald4/links/projects/def-gsarah/fridald4/col4a2_hematuria"

for pwas_model in "deCODE" "UKBPPP_EUR"
do
  for sex in "FEMALES" "MALES"
  do
    Rscript ${project_dir}/tools/BLISS/BLISS_Association.R \
      --path_sumstats ${project_dir}/data_for_bliss/X593.${sex}.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt \
      --model ${pwas_model} \
      --chr ${SLURM_ARRAY_TASK_ID} \
      --output_dir ${project_dir}/results_bliss/ \
      --output_name PWAS_BLISS_${pwas_model}_${sex}_chr${SLURM_ARRAY_TASK_ID}.txt \
      --output_augmented TRUE \
      --output_twas_fusion FALSE \
      --clean_slate FALSE
  done
done

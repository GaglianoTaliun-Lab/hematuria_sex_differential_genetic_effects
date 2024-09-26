#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --array=15
#SBATCH --job-name=bfile_bolt
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=10G

module load plink/2.00-20231024-avx2

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"
chr=${SLURM_ARRAY_TASK_ID}

# genotyped variants:
# plink2 \
#     --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen 'ref-first' \
#     --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
#     --rm-dup 'force-first' \
#     --indep-pairwise 200 50 0.9 \
#     --chr ${chr} \
#     --mac 20 \
#     --out ${project_dir}/data_for_bolt/tmp_ukb_chr${chr}_imputed

# for sex in "female" "male"
# do
#     plink2 \
#         --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${chr}_v3.bgen 'ref-first' \
#         --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
#         --keep ${data_dir}/plink_data/ukb_WB_${sex}_ids.txt \
#         --extract ${project_dir}/data_for_bolt/tmp_ukb_chr${chr}_imputed.prune.in \
#         --make-bed \
#         --threads 5 \
#         --out ${project_dir}/data_for_bolt/ukb_chr${chr}_${sex}_imputed_pruned
# done

# get 1Mb around top SNP
# chr 2: rs538848482
# chr 13: rs7323228
# chr 19: rs56254331
for sex in "male" "female"
do
    plink2 \
        --bfile ${project_dir}/data_for_bolt/ukb_chr${chr}_${sex}_imputed_pruned \
        --keep ${data_dir}/plink_data/ukb_WB_${sex}_ids.txt \
        --from-bp 45216425 \
        --to-bp 46216425 \
        --chr 15 \
        --make-bed \
        --out ${project_dir}/data_for_bolt/ukb_chr${chr}_${sex}_imputed_pruned_1Mb
done
# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Perform eGFR GWAS sex-stratified and sex-combined with the TOPMed imputed variants.

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

for sex in male female all
do
    for chr in {1..22}
    do
        run_regenie_step2="regenie --step 2 \
        --bgen ukb${data_field}_c${chr}_b0_v1.bgen \
        --sample ukb${data_field}_c${chr}_b0_v1.sample \
        --keep ukb_WB_${sex}_ids.txt \
        --phenoFile egfr_withcovariates_for_REGENIE.txt \
        --covarFile egfr_withcovariates_for_REGENIE.txt \
        --qt \
        --bsize 200 \
        --minMAC 20 \
        --minINFO 0.4 \
        --phenoCol logeGFR2009 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
        --pred step1_ukb_logeGFR2009_${sex}_pred.list \
        --threads 8 \
        --out logeGFR2009_TOPMED_${sex}_regenie_chr${chr}"

        dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
        -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
        -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
        -iin="${data_file_dir}/egfr_withcovariates_for_REGENIE.txt" \
        -iin="${data_file_dir}/step1_ukb_logeGFR2009_${sex}_pred.list" \
        -iin="${data_file_dir}/step1_ukb_logeGFR2009_${sex}_1.loco" \
        -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
    done
done

####### testing:
chr=22
sex="female"

run_regenie_step2="regenie --step 2 \
--bgen ukb${data_field}_c${chr}_b0_v1.bgen \
--sample ukb${data_field}_c${chr}_b0_v1.sample \
--keep ukb_WB_${sex}_ids.txt \
--phenoFile egfr_withcovariates_for_REGENIE.txt \
--covarFile egfr_withcovariates_for_REGENIE.txt \
--qt \
--bsize 200 \
--minMAC 20 \
--minINFO 0.4 \
--phenoCol logeGFR2009 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
--pred step1_ukb_logeGFR2009_${sex}_pred.list \
--threads 12 \
--out logeGFR2009_TOPMED_${sex}_regenie_chr${chr}"

dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
-iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
-iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
-iin="${data_file_dir}/egfr_withcovariates_for_REGENIE.txt" \
-iin="${data_file_dir}/step1_ukb_logeGFR2009_${sex}_pred.list" \
-iin="${data_file_dir}/step1_ukb_logeGFR2009_${sex}_1.loco" \
-icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes
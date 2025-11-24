# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Ranges:
# 2:227052366-227309805
# 13:110415644-110415648
# 19:41320114-41320118 #alt: 19:41319280-41319288 (for females)
# COGWI signals:
# 16:19020988-19026281
# X:75035440-75459683

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

chr="16"
range="16:19020988-19026281"

# females (using menopause status):
for group in {1..2}
do
    run_regenie_step2="regenie --step 2 \
    --bgen ukb${data_field}_c${chr}_b0_v1.bgen \
    --sample ukb${data_field}_c${chr}_b0_v1.sample \
    --keep ukb_wB_females_menopause_G${group}_ids_updated.txt \
    --range ${range} \
    --phenoFile case.control-wb.withcovariates-4regenie.txt \
    --covarFile case.control-wb.withcovariates-4regenie.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --af-cc \
    --firth \
    --firth-se \
    --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred step1_ukb_X593_menopause_status_females_G${group}_pred.list \
    --threads 8 \
    --out step2_ukb_X593_females_menopause_status_G${group}_firth_regenie_range_${range}"

    dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
    -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
    -iin="${data_file_dir}/ukb_wB_females_menopause_G${group}_ids_updated.txt" \
    -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
    -iin="${data_file_dir}/step1_ukb_X593_menopause_status_females_G${group}_pred.list" \
    -iin="${data_file_dir}/step1_ukb_X593_menopause_status_females_G${group}_1.loco" \
    -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
    --destination="/Data/results/" --brief --yes
done

# males (using groups above and below 50 yo)
for group in {1..2}
do
    run_regenie_step2="regenie --step 2 \
    --bgen ukb${data_field}_c${chr}_b0_v1.bgen \
    --sample ukb${data_field}_c${chr}_b0_v1.sample \
    --keep ukb_WB_male_ageG${group}_ids.txt \
    --range ${range} \
    --phenoFile case.control-wb.withcovariates-4regenie.txt \
    --covarFile case.control-wb.withcovariates-4regenie.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --minINFO 0.4 \
    --af-cc \
    --firth \
    --firth-se \
    --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred step1_ukb_X593_male_ageG${group}_pred.list \
    --threads 8 \
    --out step2_ukb_X593_HRC_male_ageG${group}_firth_regenie_${range}"

    dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
    -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
    -iin="${data_file_dir}/ukb_WB_male_ageG${group}_ids.txt" \
    -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
    -iin="${data_file_dir}/step1_ukb_X593_male_ageG${group}_pred.list" \
    -iin="${data_file_dir}/step1_ukb_X593_male_ageG${group}_1.loco" \
    -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
    --destination="/Data/results/" --brief --yes
done
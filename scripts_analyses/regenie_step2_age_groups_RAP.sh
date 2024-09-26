# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Ranges:
# 2:227052366-227309805
# 13:110415644-110415648
# 19:41320114-41320118 #alt: 19:41319280-41319288

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

chr=19

for sex in "female" "male"
do
    for group in {1..2}
    do
        run_regenie_step2="regenie --step 2 \
        --bgen ukb${data_field}_c${chr}_b0_v1.bgen \
        --sample ukb${data_field}_c${chr}_b0_v1.sample \
        --keep ukb_WB_${sex}_ageG${group}_ids.txt \
        --range 19:41319280-41319288 \
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
        --pred step1_ukb_X593_${sex}_G${group}_pred.list \
        --threads 8 \
        --out step2_ukb_X593_HRC_${sex}_G${group}_firth_regenie_chr${chr}"

        dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
        -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
        -iin="${data_file_dir}/ukb_WB_${sex}_ageG${group}_ids.txt" \
        -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_G${group}_pred.list" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_G${group}_1.loco" \
        -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
    done
done
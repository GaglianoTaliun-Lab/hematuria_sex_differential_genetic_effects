# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Perform regular hematuria association sex-stratified (and sex-combined) for all genome-wide significant regions from the sex-specific HRC-imputed, but using the TOPMed imputed variants.

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

for i in {1..4}
do
    chr=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/gcta/top_signals_both_sexes.txt | awk '{print $2}')
    BP=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/gcta/top_signals_both_sexes.txt | awk '{print $3}')
    from_bp=$((BP - 500000))
    to_bp=$((BP + 500000))

    for sex in male female 
    do
        run_regenie_step2="regenie --step 2 \
        --bgen ukb${data_field}_c${chr}_b0_v1.bgen \
        --sample ukb${data_field}_c${chr}_b0_v1.sample \
        --keep ukb_WB_${sex}_ids.txt \
        --phenoFile case.control-wb.withcovariates-4regenie.txt \
        --covarFile case.control-wb.withcovariates-4regenie.txt \
        --range ${chr}:${from_bp}-${to_bp} \
        --bt \
        --bsize 200 \
        --minMAC 20 \
        --minINFO 0.4 \
        --firth \
        --firth-se \
        --af-cc \
        --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
        --pred step1_ukb_X593_${sex}_pred.list \
        --threads 8 \
        --out X593_TOPMED_${sex}_regenie_firth_chr${chr}:${from_bp}-${to_bp}"

        dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
        -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
        -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
        -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_pred.list" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_1.loco" \
        -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
    done
done

# sex-combined:
run_regenie_step2="regenie --step 2 \
--bgen ukb${data_field}_c${chr}_b0_v1.bgen \
--sample ukb${data_field}_c${chr}_b0_v1.sample \
--keep ukb_WB_all_ids.txt \
--phenoFile case.control-wb.withcovariates-4regenie.txt \
--covarFile case.control-wb.withcovariates-4regenie.txt \
--range 13:110148963-110513209 \
--bt \
--bsize 200 \
--minMAC 20 \
--minINFO 0.4 \
--spa \
--phenoCol X593 --covarCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
--pred step1_ukb_X593_pred.list \
--threads 8 \
--out X593_TOPMED_regenie_chr${chr}_COL4A2"

dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
-iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
-iin="${data_file_dir}/ukb_WB_all_ids.txt" \
-iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
-iin="${data_file_dir}/step1_ukb_X593_pred.list" \
-iin="${data_file_dir}/step1_ukb_X593_1.loco" \
-icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes
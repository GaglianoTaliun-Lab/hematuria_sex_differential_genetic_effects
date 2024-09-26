# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Perform conditional hematuria association sex-stratified for all genome-wide significant regions using the TOPMed imputed variants.

# step 1: Run PLINK in RAP (TOPMed-imputed) to extract the 1Mb region and rename variant IDs, output bgen files
# step 2: get genotype (bfile) of the conditional variants from TOPMed imputed (with plink in RAP) to add as a covariate in REGENIE step 1 (to run in Béluga)
# step 3: Run REGENIE step 1 in Béluga for each sex/covariate set using the conditional variant as covariate too (in Béluga: regenie_step1_hem_cond.sh)
# step 4: Run PLINK step 2 in RAP using the output from the previous steps

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

for i in {2..9}
do
    chr=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $3}')
    BP=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $4}')
    from_bp=$((BP - 500000))
    to_bp=$((BP + 500000))
    SNP_ID=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $8}')
    sex=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $2}')

    run_plink2="plink2 \
        --bgen ukb${data_field}_c${chr}_b0_v1.bgen 'ref-last' \
        --sample ukb${data_field}_c${chr}_b0_v1.sample \
        --from-bp ${from_bp} \
        --to-bp ${to_bp} \
        --chr ${chr} \
        --set-missing-var-ids @:# \
        --export bgen-1.2 'bits=8' \
        --out TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}"

    dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
        -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
        -icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/cond_bgen/" --brief --yes
done

for i in {2..9}
do
    chr=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $3}')
    BP=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $4}')
    from_bp=$((BP - 500000))
    to_bp=$((BP + 500000))
    SNP_ID=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $8}')
    sex=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $2}')

    run_plink2="plink2 \
        --bgen TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen 'ref-last' \
        --sample TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample \
        --keep ukb_WB_${sex}_ids.txt \
        --extract SNP_to_condition_${SNP_ID}.txt \
        --chr ${chr} \
        --make-bed \
        --out genotype_conditional_${sex}_chr${chr}_${SNP_ID}"

    dx run swiss-army-knife -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen" \
        -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample" \
        -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
        -iin="${data_file_dir}/SNP_to_condition_${SNP_ID}.txt" \
        -icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
done

for i in {2..9}
do
    chr=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $3}')
    BP=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $4}')
    from_bp=$((BP - 500000))
    to_bp=$((BP + 500000))
    SNP_ID=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $8}')
    sex=$(sed -n ${i}p /Users/frida/Documents/research-projects/col4a2_hematuria/regenie/data_for_regenie/top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt | awk '{print $2}')

    run_regenie_step2="regenie --step 2 \
        --bgen TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen \
        --sample TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample \
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
        --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
        --condition-list SNP_to_condition_${SNP_ID}.txt \
        --pred step1_ukb_X593_${sex}_cond_${SNP_ID}_pred.list \
        --threads 8 \
        --out step2_X593_TOPMED_${sex}_regenie_cond_${SNP_ID}"

    dx run swiss-army-knife -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen" \
        -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample" \
        -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
        -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_cond_${SNP_ID}_pred.list" \
        -iin="${data_file_dir}/SNP_to_condition_${SNP_ID}.txt" \
        -iin="${data_file_dir}/step1_ukb_X593_${sex}_cond_${SNP_ID}_1.loco" \
        -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
done

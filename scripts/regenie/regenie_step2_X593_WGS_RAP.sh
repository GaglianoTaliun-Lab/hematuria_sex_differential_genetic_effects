# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Sanity check for the chr X cogwi signal

wgs_bgen="/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, BGEN format [500k release]/"
data_field="24309"
data_file_dir="/Data"

run_regenie_step2="regenie --step 2 \
    --bgen ukb${data_field}_c16_b0_v1.bgen \
    --sample ukb${data_field}_c16_b0_v1.sample \
    --keep sexcombined_cases.txt \
    --range 16:19009708-19076053 \
    --phenoFile case.control-wb.withcovariates-4regenie.txt \
    --covarFile case.control-wb.withcovariates-4regenie.txt \
    --bt \
    --bsize 200 \
    --minMAC 20 \
    --firth \
    --firth-se \
    --af-cc \
    --phenoCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
    --pred step1_ukb_sex_X593casesonly_pred.list \
    --threads 16 \
    --out step2_sex_X593caseonly_WGS_regenie_chr16"

dx run swiss-army-knife -iin="${wgs_bgen}/ukb${data_field}_c16_b0_v1.bgen" \
-iin="${wgs_bgen}/ukb${data_field}_c16_b0_v1.bgen.bgi" \
-iin="${wgs_bgen}/ukb${data_field}_c16_b0_v1.sample" \
-iin="${data_file_dir}/sexcombined_cases.txt" \
-iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
-iin="${data_file_dir}/step1_ukb_sex_X593casesonly_pred.list" \
-iin="${data_file_dir}/step1_ukb_sex_X593casesonly_1.loco" \
-icmd="${run_regenie_step2}" --tag="Step2" --priority high --instance-type "mem1_hdd1_v2_x8" \
--destination="/Data/results/" --brief --yes
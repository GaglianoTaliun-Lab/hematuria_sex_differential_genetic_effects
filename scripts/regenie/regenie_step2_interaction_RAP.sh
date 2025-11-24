# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

# Ranges:
# 2:227052366-227309805
# splicing variant: rs11898094 = GRCh38:2:227060058
# secondary variant COL4A3 = GRCh38: 2:227309802
# chr2 for HRC top variant: 2:226595230-226595238
# 13:110415644-110415648
# 19:41320114-41320118

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

chr=2

run_regenie_step2="regenie --step 2 \
--bgen ukb${data_field}_c${chr}_b0_v1.bgen \
--sample ukb${data_field}_c${chr}_b0_v1.sample \
--keep ukb_WB_all_ids.txt \
--range 2:227309800-227309804 \
--phenoFile case.control-wb.withcovariates-4regenie.txt \
--covarFile case.control-wb.withcovariates-4regenie.txt \
--bt \
--bsize 200 \
--minMAC 20 \
--minINFO 0.4 \
--phenoCol X593 --covarCol f.22001.0.0 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
--interaction f.22001.0.0 \
--pred step1_ukb_X593_pred.list \
--threads 8 \
--out X593_TOPMED_regenie_chr${chr}_interaction"

dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
-iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
-iin="${data_file_dir}/ukb_WB_all_ids.txt" \
-iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
-iin="${data_file_dir}/step1_ukb_X593_pred.list" \
-iin="${data_file_dir}/step1_ukb_X593_1.loco" \
-icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes
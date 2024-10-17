# 1) plink2 to convert rom bgen to bfile
# 2) plink to compute r
# 3) plink to compute Dprime

imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"
chr=13

run_plink2="plink2 \
--bgen ukb${data_field}_c${chr}_b0_v1.bgen 'ref-first' \
--sample ukb${data_field}_c${chr}_b0_v1.sample \
--extract 'bed1' bed1_file_for_plink_ld.txt \
--keep ukb_WB_all_ids.txt \
--make-bed \
--out ukb${data_field}_chr${chr}_COL4A2"

dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
-iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
-iin="${data_file_dir}/ukb_WB_all_ids.txt" \
-iin="${data_file_dir}/bed1_file_for_plink_ld.txt" \
-icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes

run_plink="plink \
--bfile ukb${data_field}_chr${chr}_COL4A2 \
--r square \
--make-just-bim \
--out COL4A2_LD_TOPMed_imputed"

dx run swiss-army-knife -iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.bed" \
-iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.bim" \
-iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.fam" \
-icmd="${run_plink}" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes

run_plink="plink \
--bfile ukb${data_field}_chr${chr}_COL4A2 \
--r2 dprime \
--make-just-bim \
--out COL4A2_LD_Dprime_TOPMed_imputed"

dx run swiss-army-knife -iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.bed" \
-iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.bim" \
-iin="/Data/results/ukb${data_field}_chr${chr}_COL4A2.fam" \
-icmd="${run_plink}" --instance-type "mem1_ssd1_v2_x16" \
--destination="/Data/results/" --brief --yes


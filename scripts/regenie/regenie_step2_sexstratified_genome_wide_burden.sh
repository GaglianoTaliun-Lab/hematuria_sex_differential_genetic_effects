# NEED TO MODIFY THE step_1*.list output from Regenie, so that the path to the LOCO file matches in REP (i.e., delete path if LOCO file is in REP Data folder).

exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
data_field="23158"
data_file_dir="/Data"

for chr in {1..22}
do

  for sex in "male" "female"
  do

    run_regenie_step2="regenie --step 2 \
      --bed ukb${data_field}_c${chr}_b0_v1 \
      --keep ukb_WB_${sex}_ids.txt \
      --phenoFile case.control-wb.withcovariates-4regenie.txt \
      --covarFile case.control-wb.withcovariates-4regenie.txt \
      --bt \
      --exclude ukb${data_field}_500k_OQFE.90pct10dp_qc_variants.txt \
      --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
      --pred step1_ukb_X593_${sex}_pred.list \
      --anno-file regenie_lof_hc_annotation.tsv \
      --set-list list_variants_sets_genome_wide.txt \
      --mask-def mask_LoF.txt \
      --aaf-bins 0.01 \
      --write-mask-snplist \
      --check-burden-files \
      --firth \
      --threads 5 \
      --out X593_WES_regenie_genome_wide_${sex}_c${chr}"

    dx run swiss-army-knife -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.bed" \
      -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.bim" \
      -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.fam" \
      -iin="${exome_file_dir}/helper_files/ukb${data_field}_500k_OQFE.90pct10dp_qc_variants.txt" \
      -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
      -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
      -iin="${data_file_dir}/step1_ukb_X593_${sex}_pred.list" \
      -iin="${data_file_dir}/step1_ukb_X593_${sex}_1.loco" \
      -iin="${data_file_dir}/regenie_lof_hc_annotation.tsv" \
      -iin="${data_file_dir}/list_variants_sets_genome_wide.txt" \
      -iin="${data_file_dir}/mask_LoF.txt" \
      -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
      --destination="/Data/results/" --brief --yes

  done
done

for chr in "X"
do

  for sex in "male" "female"
  do

    run_regenie_step2="regenie --step 2 \
      --bed ukb${data_field}_c${chr}_b0_v1 \
      --keep ukb_WB_${sex}_ids.txt \
      --phenoFile case.control-wb.withcovariates-4regenie.txt \
      --covarFile case.control-wb.withcovariates-4regenie.txt \
      --bt \
      --exclude ukb${data_field}_500k_OQFE.90pct10dp_qc_variants.txt \
      --phenoCol X593 --covarCol f.34.0.0 --covarCol f.22009.0.{1:10} \
      --pred step1_ukb_X593_${sex}_pred.list \
      --anno-file regenie_lof_hc_annotation.tsv \
      --set-list list_variants_sets_genome_wide.txt \
      --mask-def mask_LoF.txt \
      --aaf-bins 0.01 \
      --write-mask-snplist \
      --check-burden-files \
      --firth \
      --threads 5 \
      --out X593_WES_regenie_genome_wide_${sex}_c${chr}"

    dx run swiss-army-knife -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.bed" \
      -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.bim" \
      -iin="${exome_file_dir}/ukb${data_field}_c${chr}_b0_v1.fam" \
      -iin="${exome_file_dir}/helper_files/ukb${data_field}_500k_OQFE.90pct10dp_qc_variants.txt" \
      -iin="${data_file_dir}/ukb_WB_${sex}_ids.txt" \
      -iin="${data_file_dir}/case.control-wb.withcovariates-4regenie.txt" \
      -iin="${data_file_dir}/step1_ukb_X593_${sex}_pred.list" \
      -iin="${data_file_dir}/step1_ukb_X593_${sex}_1.loco" \
      -iin="${data_file_dir}/regenie_lof_hc_annotation.tsv" \
      -iin="${data_file_dir}/list_variants_sets_genome_wide.txt" \
      -iin="${data_file_dir}/mask_LoF.txt" \
      -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16" \
      --destination="/Data/results/" --brief --yes

  done
done

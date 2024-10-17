imp_file_dir="/Bulk/Imputation/Imputation from genotype (TOPmed)/"
data_field="21007"
data_file_dir="/Data"

chr=2
from_bp=227052357
to_bp=227309812

run_plink2="plink2 \
    --bgen ukb${data_field}_c${chr}_b0_v1.bgen 'ref-first' \
    --sample ukb${data_field}_c${chr}_b0_v1.sample \
    --keep ukb_WB_all_ids.txt \
    --from-bp ${from_bp} \
    --to-bp ${to_bp} \
    --chr ${chr} \
    --set-missing-var-ids @:# \
    --export bgen-1.2 'bits=8' \
    --out TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}"

dx run swiss-army-knife -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.bgen" \
    -iin="${imp_file_dir}/ukb${data_field}_c${chr}_b0_v1.sample" \
    -iin="/Data/ukb_WB_all_ids.txt" \
    -icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
    --destination="/Data/results/" --brief --yes

for variant in "2:227277511" "2:227060058" "2:227077364" "2:227309802"
do
    run_plink2="plink2 \
        --bgen TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen 'ref-first' \
        --sample TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample \
        --snp ${variant} \
        --freq \
        --out MAF_TOPMed_imputed_wB_chr${variant}"

    dx run swiss-army-knife -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen" \
        -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample" \
        -icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
done

for variant in "2:227277511" "2:227060058" "2:227077364"
do
    run_plink2="plink2 \
        --bgen TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen 'ref-first' \
        --sample TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample \
        --ld 2:227309802 ${variant} \
        --out LD_TOPMed_imputed_wB_chr2:227309802_chr${variant}"

    dx run swiss-army-knife -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.bgen" \
        -iin="/Data/results/TOPMed_imputed_chr${chr}_${from_bp}-${to_bp}.sample" \
        -icmd="${run_plink2}" --instance-type "mem1_ssd1_v2_x16" \
        --destination="/Data/results/" --brief --yes
done
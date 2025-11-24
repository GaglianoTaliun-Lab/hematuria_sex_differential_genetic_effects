# Description: run colocalization analysis between sex-stratified NEPTUNE eQTLs and sex-stratified hematuria
# as a follow-up of TWAS signals
# Narval

# Packages -------------------------------------------------------

library(coloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/hematuria_sexspecific"
base_dir <- "/home/fridald4/projects/def-gsarah"

args <- commandArgs(TRUE)
tissue <- as.character(args[1])
sex <- as.character(args[2])

cat("Tissue = ", tissue,"\n")
cat("Sex = ", sex,"\n")

# Coloc priors
p1 = 1e-04
p2 = 1e-04
p12 = 5e-06

coloc_results_summ = list()
coloc_results_res = list()
results_names <- array()

test_loci <- read.table(here(project_dir, "data_for_coloc", stringr::str_c("locus_to_test_GRCh37_", tissue, "_", sex, ".txt")), sep = "\t", header = TRUE)

eqtls_N <- read.table(here(project_dir, "data_for_coloc", "neptune_sample_sizes.txt"), sep = "\t", header = TRUE) %>%
    filter(., tissue == tissue) %>%
    filter(., sex == sex) %>%
    dplyr::select(N)

if (sex == "FEMALE") {
    hematuria_case_ratio = 7061/213033
} else if (sex == "MALE") {
    hematuria_case_ratio = 9774/177666
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read datasets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# hematuria GWAS
cat("Reading hematuria GWAS...\n")
gwas <- fread(here(project_dir, "hematuria_sumstats", stringr::str_c("X593.", sex, "S.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"))) %>%
    mutate(MAF = case_when(
        EAF > 0.5 ~ 1-EAF,
        EAF <= 0.5 ~ EAF),
        rsid = SNP, SNP = stringr::str_c(CHR, "_", BP, "_", other_allele, "_", effect_allele), varbeta = SE^2) %>%
    filter(., MAF > 0.01 & MAF < 1) %>%
    filter(., !is.na(SNP)) %>%
    distinct(., SNP, .keep_all = TRUE) %>%
    dplyr::select(CHR, BP, SNP, BETA, varbeta, P, MAF, N)
cat("Done.\n")

# eQTLs allele frequencies (to get MAF)
cat("Reading frequency files...\n")

tissue_lower = tolower(tissue)

if (sex == "FEMALE") {
    freq <- fread(here(base_dir, "neptune-sex-specific-eQTL", stringr::str_c(tissue_lower, "_freq_bySex.tsv"))) %>%
        mutate(MAF_female = case_when(
            alt_freq_female > 0.5 ~ 1-alt_freq_female,
            alt_freq_female <= 0.5 ~ alt_freq_female
    ))
maf_freq <- freq %>% dplyr::select(SNP = chr_pos_ref_alt, MAF = MAF_female)

} else if (sex == "MALE") {
    freq <- fread(here(base_dir, "neptune-sex-specific-eQTL", stringr::str_c(tissue_lower, "_freq_bySex.tsv"))) %>%
        mutate(MAF_male = case_when(
            alt_freq_male > 0.5 ~ 1-alt_freq_male,
            alt_freq_male <= 0.5 ~ alt_freq_male
    ))
maf_freq <- freq %>% dplyr::select(SNP = chr_pos_ref_alt, MAF = MAF_male)
}

rm(freq)
cat("Done.\n")

# eQTLs sumstats
cat("Reading eQTL summary statistics...\n")
eqtls <- fread(here(base_dir, "neptune-sex-specific-eQTL", stringr::str_c("MatrixEQTL_" , tissue, "_", sex, ".csv.gz"))) %>%
    left_join(., maf_freq, by = "SNP") %>%
    filter(., MAF > 0.01 & MAF < 1) %>%
    filter(., !is.na(SNP)) %>%
    distinct(., SNP, .keep_all = TRUE) %>%
    separate(SNP, c("CHR", "BP", "A2", "A1"), sep = "_", remove = FALSE) %>%
    mutate(N = eqtls_N$N[1], Zscore = sign(beta) * abs(qnorm(`p-value`/2)), SE = beta/Zscore, varbeta = SE^2) %>%
    dplyr::select(CHR, BP, SNP, BETA = beta, varbeta, P = `p-value`, MAF, N)
cat("Done.\n")

for (i in 1:nrow(test_loci)){

    chromosome = test_loci$CHR[i]
    start_pos = test_loci$START[i]
    end_pos = test_loci$END[i]
    ensembl_id = test_loci$ENSEMBL_ID[i]

    # subset for a region only:
    gwas_region <- gwas %>%
        filter(., CHR == chromosome & BP >= start_pos & BP <= end_pos)

    # subset for a region only:
    eqtls_region <- eqtls %>%
        filter(., CHR == chromosome & BP >= start_pos & BP <= end_pos)

    write.table(gwas_region, here(project_dir, "data_for_coloc", stringr::str_c("regional_plot_data_gwas_", tissue, "_", sex, "_", ensembl_id, ".txt")), sep = "\t", row.names = F, quote = F)
    write.table(eqtls_region, here(project_dir, "data_for_coloc", stringr::str_c("regional_plot_data_eqtls_", tissue, "_", sex, "_", ensembl_id, ".txt")), sep = "\t", row.names = F, quote = F)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run coloc
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                            snp = gwas_region$SNP,
                                            beta = gwas_region$BETA,
                                            varbeta = gwas_region$varbeta,
                                            pvalues = gwas_region$P,
                                            MAF = gwas_region$MAF,
                                            N = gwas_region$N,
                                            s = hematuria_case_ratio),
                            dataset2 = list(type = "quant",
                                            snp = eqtls_region$SNP,
                                            beta = eqtls_region$BETA,
                                            varbeta = eqtls_region$varbeta,
                                            pvalues = eqtls_region$P,
                                            MAF = eqtls_region$MAF,
                                            N = eqtls_region$N),
                            p1 = p1, p2 = p2, p12 = p12
    )

        coloc_results_summ[[i]] <- coloc_results$summary
        coloc_results_res[[i]] <- coloc_results$results
        results_names[i] <- ensembl_id

}

coloc_results_summ <- setNames(coloc_results_summ, nm = results_names)
coloc_results_res <- setNames(coloc_results_res, nm = results_names)

coloc_results_summ %>% dplyr::bind_rows(., .id = "Ensembl_ID") %>%
    write.table(., here(project_dir, "results_coloc", stringr::str_c("coloc_summary_", tissue, "_", sex, ".txt")), sep = "\t", row.names = F, quote = F)

coloc_results_res %>% dplyr::bind_rows(., .id = "Ensembl_ID") %>%
    write.table(., here(project_dir, "results_coloc", stringr::str_c("coloc_all_", tissue, "_", sex, ".txt")), sep = "\t", row.names = F, quote = F)

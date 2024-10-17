# eGFR sex by genotype interaction analysis

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(R.utils)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
dbSNP144 = SNPlocs.Hsapiens.dbSNP144.GRCh38

# Main --------------------------------------

topmed_files_m <- list.files(here(project_dir, "regenie", "results"), pattern = "X593_TOPMED_male_regenie_firth_*", full.names = T) %>%
  stringr::str_subset(., "_X593.regenie")
topmed_files_f <- list.files(here(project_dir, "regenie", "results"), pattern = "X593_TOPMED_female_regenie_firth*", full.names = T) %>%
  stringr::str_subset(., "_X593.regenie")

topmed_files <- c(topmed_files_m, topmed_files_f)

topmed_sumstats <- lapply(topmed_files, function(x) read.table(x, sep = " ", header = T) %>%
                            filter(., !is.na(LOG10P)))
topmed_names <- topmed_files %>% stringr::str_remove("/Users/frida/Documents/research-projects/col4a2_hematuria/regenie/results/") %>%
  stringr::str_remove("_X593.regenie")

top_snps <- data.frame(dataset = topmed_names, CHR = seq(1, length(topmed_names)), BP = seq(1, length(topmed_names)), A0 = seq(1, length(topmed_names)), A1 = seq(1, length(topmed_names)), Freq_A1 = seq(1, length(topmed_names)), beta = seq(1, length(topmed_names)), se = seq(1, length(topmed_names)), pval = seq(1, length(topmed_names)))

for (i in 1:length(topmed_names)) {
  
  data <- arrange(topmed_sumstats[[i]], -LOG10P)
  top_snps[i,2] <- data$CHROM[1]
  top_snps[i,3] <- data$GENPOS[1]
  top_snps[i,4] <- data$ALLELE1[1]
  top_snps[i,5] <- data$ALLELE0[1]
  top_snps[i,6] <- 1-data$A1FREQ[1]
  top_snps[i,7] <- data$BETA[1] * -1
  top_snps[i,8] <- data$SE[1]
  top_snps[i,9] <- data$LOG10P[1]
  
}

top_snps <- top_snps %>%
  mutate(SNP_ID = stringr::str_c(CHR,"_",BP,"_",A1,"_",A0), SNP_ID2 = stringr::str_c(CHR,":",BP,":",A1,":",A0), P = 10^-pval) %>%
         separate(., dataset, into = c(NA, NA, "sex", NA, NA, NA), sep = "_", remove = F) %>%
  colochelpR::convert_loc_to_rs(., dbSNP144)
write.table(top_snps %>% dplyr::select(-pval), here(project_dir, "regenie", "results", "conditional_TOPMed_imputed", "top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)

top_snps <- top_snps %>%
  dplyr::select(dataset, sex, CHR, BP, A0 = A1, A1 = A0, pval, SNP_ID, SNP_ID2)
write.table(top_snps, here(project_dir, "regenie", "data_for_regenie", "top_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed.txt"), sep = "\t", row.names = F, quote = F)

for (i in 1:nrow(top_snps)) {
  write.table(top_snps$SNP_ID2[i] %>% stringr::str_extract(., "[:digit:]+:[:digit:]+"), here(project_dir, "regenie", "data_for_regenie", stringr::str_c("SNP_to_condition_", top_snps$SNP_ID[i], ".txt")), sep = "\t", row.names = F, col.names = F, quote = F)
}


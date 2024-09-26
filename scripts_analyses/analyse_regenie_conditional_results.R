# analyse regenie conditional results with topmed imputation

# Load Packages -----------------------------

library(here)
library(stringr)
library(data.table)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
dbSNP144 = SNPlocs.Hsapiens.dbSNP144.GRCh38

# Main --------------------------------------

step2_files <- list.files(here(project_dir, "regenie", "results", "conditional_TOPMed_imputed"), pattern = "*.regenie", full.names = T) %>%
  stringr::str_subset(., ".log", negate = TRUE)
step2_names <- list.files(here(project_dir, "regenie", "results", "conditional_TOPMed_imputed"), pattern = "*.regenie", full.names = F) %>%
  stringr::str_subset(., ".log", negate = TRUE) %>%
  stringr::str_extract(., "_[:lower:]+_") %>%
  stringr::str_remove_all(., "_")

step2_results <- setNames(lapply(step2_files, function(x) read.table(x, sep = " ", header = T) %>%
                          mutate(CHR = CHROM, BP = GENPOS, P = 10^-LOG10P) %>%
                          colochelpR::convert_loc_to_rs(., dbSNP144) %>%
                          filter(., P < 5e-08)
                        ), nm = step2_names) %>%
  dplyr::bind_rows(., .id = "sex") %>%
  mutate(A0 = ALLELE1, A1 = ALLELE0, Freq_A1 = 1-A1FREQ, BETA = BETA*-1,
         SNP = case_when(
           is.na(SNP) ~ stringr::str_c(CHR,":",BP),
           !is.na(SNP) ~ SNP
         )) %>%
  dplyr::select(sex, CHR, BP, SNP, ALLELE0 = A0, ALLELE1 = A1, A1FREQ = Freq_A1, INFO, N, BETA, SE, CHISQ, P)

write.table(step2_results, here(project_dir, "regenie", "results", "conditional_TOPMed_imputed", "top_conditional_SNPs_GRCh38_X593_males_and_females_TOPMed_imputed_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)



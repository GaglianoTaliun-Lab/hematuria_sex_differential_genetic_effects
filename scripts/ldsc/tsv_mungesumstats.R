# create tsv file for munge_sumstats
# it removes the N column (which causes an error in munge_sumstats.py script)
# it removes p-values==NA
# it removes variants without rsid
# keep columns: SNP, effect_allele, other_allele, BETA, SE, P 

#------------------------------------------------------- Packages ------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

#------------------------------------------------------- Arguments -----------------------------------------------------------

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

sumstats_files <- c(stringr::str_c(project_dir,"/datasets/", "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"),
                    stringr::str_c(project_dir,"/datasets/", "GCST90301699.h.DBP_females.tsv.gz"),
                    stringr::str_c(project_dir,"/datasets/", "GCST90301702.h.PP_females.tsv.gz"),
                    stringr::str_c("/home/fridald4/projects/def-gsarah", "/Aragam_2022_CARDIoGRAM_CAD_GWAS/", 
                      "GWAMA.CAD.INTERMEDIATE.SENSITIVITY.UKBB.PLUS.EXTRA.STUDIES.QCed.txt.gz"),
                    stringr::str_c(project_dir,"/datasets/", "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"),
                    stringr::str_c(project_dir, "/results_regenie/", "eGFR_male_HRC_imputed_regenie_allchrs.regenie"),
                    stringr::str_c(project_dir, "/results_regenie/", "eGFR_female_HRC_imputed_regenie_allchrs.regenie"))

#------------------------------------------------------- Main ----------------------------------------------------------------

# 1) Hematuria females:
fread(sumstats_files[1]) %>%
  filter(., !is.na(P)) %>%
      filter(., SNP != "") %>%
      mutate(P = as.numeric(P)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP,
    A1 = effect_allele,
    A2 = other_allele,
    BETA,
    SE,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_hematuria_females.tsv"), sep = "\t", row.names = F, quote = F)

# 2) DBP females:
fread(sumstats_files[2]) %>%
  filter(., !is.na(p_value)) %>%
      filter(., rsid != "") %>%
      mutate(P = as.numeric(p_value)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP = rsid,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_DBP_females.tsv"), sep = "\t", row.names = F, quote = F)

# 3) PP females:
fread(sumstats_files[3]) %>%
  filter(., !is.na(p_value)) %>%
      filter(., rsid != "") %>%
      mutate(P = as.numeric(p_value)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP = rsid,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_PP_females.tsv"), sep = "\t", row.names = F, quote = F)

# 4) CAD females:
fread(sumstats_files[4]) %>%
  filter(., !is.na(female_p_value)) %>%
      filter(., rsid_ukb != "") %>%
      mutate(P = as.numeric(female_p_value)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP = rsid_ukb,
    A1 = reference_allele,
    A2 = other_allele,
    BETA = female_beta,
    SE = female_se,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_CAD_females.tsv"), sep = "\t", row.names = F, quote = F)

# 5) Hematuria males:
fread(sumstats_files[5]) %>%
  filter(., !is.na(P)) %>%
      filter(., SNP != "") %>%
      mutate(P = as.numeric(P)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP,
    A1 = effect_allele,
    A2 = other_allele,
    BETA,
    SE,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_hematuria_males.tsv"), sep = "\t", row.names = F, quote = F)

# 6) eGFRcrea males:
fread(sumstats_files[6]) %>%
  mutate(P = 10^-(LOG10P)) %>%
  filter(., !is.na(P)) %>%
      filter(., ID != "") %>%
      mutate(P = as.numeric(P)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP = ID,
    A1 = ALLELE1,
    A2 = ALLELE0,
    BETA,
    SE,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_eGFRcrea_males.tsv"), sep = "\t", row.names = F, quote = F)

# 7) eGFRcrea females:
fread(sumstats_files[7]) %>%
  mutate(P = 10^-(LOG10P)) %>%
  filter(., !is.na(P)) %>%
      filter(., ID != "") %>%
      mutate(P = as.numeric(P)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
  dplyr::select(
    SNP = ID,
    A1 = ALLELE1,
    A2 = ALLELE0,
    BETA,
    SE,
    P
  ) %>%
  write.table(., here(project_dir,"data_for_ldsc", "tmp_eGFRcrea_females.tsv"), sep = "\t", row.names = F, quote = F)
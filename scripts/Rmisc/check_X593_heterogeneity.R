# Cancer codes in UKB -
# urological cancers in males ICD9:
#   187, 1874, 1878, 185, 1859, 186, 1869, 1860
# gynecological cancers in females ICD9:
#   179, 1799, 1838, 1839, 1840, 180, 1808, 1809, V762, 1844

# urological cancers in males ICD10: 
#   C61 = malignant neoplasm of prostate
#   C62 = malignant neoplasm of testis
#   C60 = malignant neoplasm of penis
# gynecological cancer in females ICD10: 
#   C56 + C796 = malignant neoplasm of ovary/secondary malignant neoplasm of ovary
#   C55 = malignant neoplasm of uterus
#   C52 = malignant neoplasm of vagina
#   C53 = malignant neoplasm of cervix uterine
#   C51 = malignant neoplasm of vulva

# https://biobank.ndph.ox.ac.uk/~bbdatan/CancerSummaryReport.html

# Kidney stones in UKB -
# Extracted the ICD9 (‘calculus of kidney and ureter [592] and calculus of lower urinary tract [594]’) 
# and ICD10 (‘urolithiasis’ [N20-N23]) codes). 

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ----------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

cancer_males <- c("C61", "C62", "C60", 187, 1874, 1878, 185, 1859, 186, 1869, 1860)
cancer_females <- c("C56", "C55", "C52", "C53", "C51", 179, 1799, 1838, 1839, 1840, 180, 1808, 1809, "V762", 1844)
cancer_all <- c(cancer_males, cancer_females)

# Read files ---------------------------------

ICD9 <- fread(here(project_dir, "cancer_UKB_RAP", "cancer_ICD9_f40013.txt"))
ICD10 <- fread(here(project_dir, "cancer_UKB_RAP", "cancer_ICD10_f40006.txt"), na.strings = "")
# dates <- fread(here(project_dir, "cancer_UKB_RAP", "cancer_dates_f40005.txt"))

ICD9_ks <- fread(here(project_dir, "renal_codes_UKB_RAP", "urolithiasis_ICD9_incidence.txt")) %>%
  dplyr::filter(., urolithiasis == 1)
ICD10_ks <- fread(here(project_dir, "renal_codes_UKB_RAP", "urolithiasis_ICD10_incidence.txt")) %>%
  dplyr::filter(., urolithiasis == 1)

X593 <- fread(here(project_dir, "case.control-wb.withcovariates-4regenie.txt"))

ICD10_codes <- fread(here(project_dir, "cancer_UKB_RAP", "coding_ICD10.tsv")) %>%
  dplyr::select(cancer_code_ICD10 = coding, meaning_ICD10 = meaning)
ICD9_codes <- fread(here(project_dir, "cancer_UKB_RAP", "coding_ICD9.tsv")) %>%
  dplyr::select(cancer_code_ICD9 = coding, meaning_ICD9 = meaning)

exclude_ids <- fread("/Users/frida/Documents/research-projects/misc/ukb/w66222_20260310.csv")

# Main ---------------------------------------

X593_cases <- X593 %>%
  dplyr::filter(., X593 == 1)

sex_code <- X593 %>%
  dplyr::select(eid = IID, f.22001.0.0)
# females = 0, males = 1

### cancer codes:
ICD9_593cases <- ICD9 %>%
  left_join(., sex_code, by = "eid") %>%
  dplyr::filter(eid %in% X593_cases$IID) %>%
  dplyr::filter(!if_all(c(p40013_i0:p40013_i13), is.na)) %>%
  filter(if_any(everything(), ~ .x %in% cancer_all)) %>%
  mutate(across(.cols = starts_with("p40013_"),
                .fns = ~ case_when(
                  !(. %in% cancer_all) ~ NA,
                  TRUE ~ as.character(.)
                ))) %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  mutate(cancer_code_ICD9 = coalesce(!!!select(., starts_with("p40013_")))) %>%
  dplyr::select(eid, f.22001.0.0, cancer_code_ICD9)

ICD10_593cases <- ICD10 %>%
  left_join(., sex_code, by = "eid") %>%
  dplyr::filter(eid %in% X593_cases$IID) %>%
  dplyr::filter(!if_all(c(p40006_i0:p40006_i20), is.na)) %>%
  filter(if_any(everything(), ~ .x %in% cancer_all)) %>%
  mutate(across(.cols = starts_with("p40006_"),
                .fns = ~ case_when(
                  !(. %in% cancer_all) ~ NA,
                  TRUE ~ as.character(.)
                ))) %>%
  dplyr::select(where(~ !all(is.na(.))))  %>%
  mutate(cancer_code_ICD10 = coalesce(!!!select(., starts_with("p40006_")))) %>%
  dplyr::select(eid, f.22001.0.0, cancer_code_ICD10)

X593_cases_all <- full_join(ICD9_593cases, ICD10_593cases, by = c("eid", "f.22001.0.0"))

counts_codes <- count(X593_cases_all, cancer_code_ICD9, cancer_code_ICD10) %>%
  left_join(., ICD10_codes, by = "cancer_code_ICD10") %>%
  left_join(., ICD9_codes, by = "cancer_code_ICD9") %>%
  dplyr::select("cancer_code_ICD9", "meaning_ICD9", "cancer_code_ICD10", "meaning_ICD10", "n")

fwrite(counts_codes, here(project_dir, "cancer_UKB_RAP", "counts_sex_specific_cancer_incidences_among_hematuria_cases.tsv"), sep = "\t", na = "NA")

### kidney stones codes:
ICD9_ks_593cases <- ICD9_ks %>%
  left_join(., sex_code, by = "eid") %>%
  dplyr::filter(eid %in% X593_cases$IID)

ICD10_ks_593cases <- ICD10_ks %>%
  left_join(., sex_code, by = "eid") %>%
  dplyr::filter(eid %in% X593_cases$IID)

overlap <- ICD10_593cases %>%
  full_join(., ICD9_593cases, by = c("eid", "f.22001.0.0")) %>%
  full_join(., ICD10_ks_593cases, by = c("eid", "f.22001.0.0")) %>%
  full_join(., ICD9_ks_593cases, by = c("eid", "f.22001.0.0", "urolithiasis"))

# in total, there are 1,452 hematuria cases that have a urolithiasis ICD code. Of those, 1,086 are males.

# ----------------------- create a new hematuria phenotype/covariate file for Regenie:

### cancer codes:
ICD9_nocancer_keep <- ICD9 %>%
  dplyr::rename(IID = eid) %>%
  dplyr::right_join(., X593, by = "IID") %>%
  mutate(across(.cols = starts_with("p40013_"),
                .fns = ~ case_when(
                  !(. %in% cancer_all) ~ NA,
                  TRUE ~ as.character(.)
                ))) %>%
  mutate(cancer_code_ICD9 = coalesce(!!!select(., starts_with("p40013_")))) %>%
  dplyr::filter(., is.na(cancer_code_ICD9)) 

#  dplyr::select(FID, IID, X593, f.34.0.0, f.22001.0.0, f.22009.0.1, f.22009.0.2, f.22009.0.3, f.22009.0.4, 
#                f.22009.0.5, f.22009.0.6, f.22009.0.7, f.22009.0.8, f.22009.0.9, f.22009.0.10)

ICD10_nocancer_keep <- ICD10 %>%
  dplyr::rename(IID = eid) %>%
  dplyr::right_join(., X593, by = "IID") %>%
  mutate(across(.cols = starts_with("p40006_"),
                .fns = ~ case_when(
                  !(. %in% cancer_all) ~ NA,
                  TRUE ~ as.character(.)
                ))) %>%
  mutate(cancer_code_ICD10 = coalesce(!!!select(., starts_with("p40006_")))) %>%
  dplyr::filter(., is.na(cancer_code_ICD10))

### final: keep ICD10 and ICD9 no cancer inds., remove inds. with kidney stones, remove withdrawn inds (2026):

X593_homog <- X593 %>%
  dplyr::filter(IID %in% ICD9_nocancer_keep$IID) %>%
  dplyr::filter(IID %in% ICD10_nocancer_keep$IID) %>%
  dplyr::filter(!(IID %in% ICD9_ks$eid)) %>%
  dplyr::filter(!(IID %in% ICD10_ks$eid)) %>%
  dplyr::filter(!(IID %in% exclude_ids$V1))

fwrite(X593_homog, here(project_dir, "case.control-wb.withcovariates-homogenous-4regenie.txt"), sep = "\t", na = "NA")



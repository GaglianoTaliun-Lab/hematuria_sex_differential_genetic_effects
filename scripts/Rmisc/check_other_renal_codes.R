
library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ----------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

# Read files ---------------------------------

ICD9_renal <- fread(here(project_dir, "renal_codes_UKB_RAP", "other_renal_codes_ICD9_incidence.txt")) %>%
  dplyr::rename(IID = eid)
ICD10_renal <- fread(here(project_dir, "renal_codes_UKB_RAP", "other_renal_codes_ICD10_incidence.txt")) %>%
  dplyr::rename(IID = eid)

X593 <- fread(here(project_dir, "case.control-wb.withcovariates-4regenie.txt"))

# Main ---------------------------------------

X593_cases <- X593 %>%
  dplyr::filter(., X593 == 1)

sex_code <- X593 %>%
  dplyr::select(IID, f.22001.0.0)
# females = 0, males = 1

ICD9_renal_X593cases <- ICD9_renal %>%
  dplyr::filter(., IID %in% X593_cases$IID) %>%
  left_join(., sex_code, by = "IID")

ICD10_renal_X593cases <- ICD10_renal %>%
  dplyr::filter(., IID %in% X593_cases$IID) %>%
  left_join(., sex_code, by = "IID")

ICD9_renal_all <- ICD9_renal %>%
  dplyr::filter(., IID %in% X593$IID) %>%
  left_join(., X593, by = "IID")

ICD10_renal_all <- ICD10_renal %>%
  dplyr::filter(., IID %in% X593$IID) %>%
  left_join(., X593, by = "IID")


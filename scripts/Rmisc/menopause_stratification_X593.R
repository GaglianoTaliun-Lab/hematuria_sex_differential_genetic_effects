# Get menopause status to stratify for X593 females

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ---------------------------------

project_dir = "/home/fridald4/links/projects/def-gsarah/fridald4/col4a2_hematuria"

# Read field files --------------------------------------

# hematuria status
X593 <- read.table(here(project_dir, "data_for_regenie", "case.control-wb.withcovariates-4regenie.txt"), sep = "\t", header = TRUE)

# white British female IDs
wb_females <- read.table(here(project_dir, "data_for_regenie","ukb_WB_female_ids.txt"), sep = "\t", header = F) %>%
    dplyr::rename(FID = V1, IID = V2)

# age
age <- read.table(here(project_dir, "f21022.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid, age = f.21022.0.0)

# filter hematuria status to keep only females and cases
X593_females <- inner_join(X593, wb_females, by = c("FID", "IID")) %>%
  left_join(., age, by = "FID") 
  
X593_females_cases <- X593_females %>%
  filter(., X593 == 1)

# age at menopause
age_menopause <- read.table(here(project_dir, "f3581.tsv"), sep = "\t", header = T) %>%
  dplyr::rename(FID = eid)

# menopause status
menopause_status <- read.table(here(project_dir, "f2724_all_instances.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid)

# age at diagnosis X593 (oldest time stamp - select min. age at diagnosis)
age_diagnosis <- read.table(here(project_dir, "data_for_regenie", "age_at_diagnosis_cases_X593_updated.txt"), sep = "\t", header = T) %>%
  dplyr::select(FID = IID, age_at_time_stamp = age_at_X593_diagnosis)

# Main ----------------------------------------------------

## Age at menopause - 4 instances max. registered per participant - get the most recent one
age_menopause_spread <- age_menopause %>%
  tidyr::gather(instance, age_menopause, 2:5) %>%
  filter(., age_menopause > 1) %>% # removes values of -1 and -3
  tidyr::spread(., instance, age_menopause)

age_menop_inst4 <- filter(age_menopause_spread, !is.na(X3581.3.0)) %>%
  dplyr::select(FID, X3581.latest = X3581.3.0)

age_menop_inst3 <- filter(age_menopause_spread, !is.na(X3581.2.0)) %>%
  filter(., !(FID %in% age_menop_inst4$FID)) %>%
  dplyr::select(FID, X3581.latest = X3581.2.0)

age_menop_inst2 <- filter(age_menopause_spread, !is.na(X3581.1.0)) %>%
  filter(., !(FID %in% c(age_menop_inst4$FID, age_menop_inst3$FID))) %>%
  dplyr::select(FID, X3581.latest = X3581.1.0)

age_menop_inst1 <- filter(age_menopause_spread, !is.na(X3581.0.0)) %>%
  filter(., !(FID %in% c(age_menop_inst4$FID, age_menop_inst3$FID, age_menop_inst2$FID))) %>%
  dplyr::select(FID, X3581.latest = X3581.0.0)

age_menopause_latest <- rbind(age_menop_inst1, age_menop_inst2, age_menop_inst3, age_menop_inst4) %>%
  dplyr::rename(f.3581 = X3581.latest)

# menopause status - merge all instances (values of 3 and 2 answered 'not sure', 1 answered they have had menopause and 0 they haven't had)
menopause_status <- menopause_status %>%
  mutate(f.2724.final = case_when(
    f.2724.0.0 == 3 | f.2724.1.0 == 3 | f.2724.2.0 == 3 | f.2724.3.0 == 3 ~ NA,
    f.2724.0.0 == 2 | f.2724.1.0 == 2 | f.2724.2.0 == 2 | f.2724.3.0 == 2 ~ NA,
    f.2724.0.0 == 1 | f.2724.1.0 == 1 | f.2724.2.0 == 1 | f.2724.3.0 == 1 ~ 1,
    f.2724.0.0 == 0 | f.2724.1.0 == 0 | f.2724.2.0 == 0 | f.2724.3.0 == 0 ~ 0,
    TRUE ~ NA
  )) %>%
  filter(., FID %in% wb_females$FID) %>%
  dplyr::select(FID, f.2724.final)

# Merge with hematuria status - use age at oldest hematuria time stamp and menopause status
all <- X593_females %>%
  left_join(., age_diagnosis, by = "FID") %>%
  left_join(., age_menopause_latest, by = "FID") %>%
  left_join(., menopause_status, by = "FID") %>%
  # if participant reported age at menopause (f.3581), but their menopause status (f.2724.final) = 0, change menopause status to 1.
  # if participant reported age at menopause (f.3581) values of -1 (prefer not to answer) or -3 (do not know), change menopause status to NA.
  mutate(f.2724.final = case_when(
    f.3581 > 0 & f.2724.final == 0 ~ 1,
    f.3581 < 0 ~ NA,
    TRUE ~ f.2724.final
  ))

# fix menopause status based on the age at oldest ICD9/ICD10 time stamp (age_at_X593)
# in case the diagnosis was done before menopause, but participant has since reported age at menopause
all <- all %>%
  mutate(menopause_status_X593 = case_when(
    is.na(f.2724.final) ~ NA,
    f.2724.final == 0 ~ 0,
    X593 == 0 ~ f.2724.final,
    f.2724.final == 1 & X593 == 1 & (f.3581 <= age_at_time_stamp) == TRUE ~ 1,
    f.2724.final == 1 & X593 == 1 & (f.3581 > age_at_time_stamp) == TRUE ~ 0
  ))

filter(all, menopause_status_X593 == 0) %>%
  dplyr::select(FID, IID) %>%
  write.table(here(project_dir, "data_for_regenie", "ukb_wB_females_menopause_G1_ids_updated.txt"), sep = "\t", row.names = F, quote = F)

filter(all, menopause_status_X593 == 1) %>%
  dplyr::select(FID, IID) %>%
  write.table(here(project_dir, "data_for_regenie", "ukb_wB_females_menopause_G2_ids_updated.txt"), sep = "\t", row.names = F, quote = F)

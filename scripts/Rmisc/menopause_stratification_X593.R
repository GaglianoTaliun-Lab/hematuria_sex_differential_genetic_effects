# Get menopause status to stratify for X593 females

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ---------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
data_dir="/home/fridald4/projects/def-gsarah/fridald4/ukb-wes-hematuria"

# Read field files --------------------------------------

# hematuria status
X593 <- read.table(here(project_dir, "data_for_regenie", "case.control-wb.withcovariates-4regenie.txt"), sep = "\t", header = TRUE)

# white British female IDs
wb_females <- read.table(here(data_dir,"plink_data","ukb_WB_female_ids.txt"), sep = "\t", header = F) %>%
    dplyr::rename(FID = V1, IID = V2)

# filter hematuria status to keep only females and cases
X593_females <- inner_join(X593, wb_females, by = c("FID", "IID")) %>%
  left_join(., age, by = "FID") 
  
X593_females_cases <- X593_females %>%
  filter(., X593 == 1)

# age
age <- read.table(here(project_dir, "f21022.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid, age = f.21022.0.0)

# age at menopause
age_menopause <- read.table(here(project_dir, "f3581.tsv"), sep = "\t", header = T) %>%
  dplyr::rename(FID = eid)

# time stamp for ICD10 - keep female cases only
time_stamp_ICD10 <- read.table(here(project_dir, "f41280.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid) %>%
  filter(., FID %in% X593_females_cases$FID)
time_stamp_ICD10[2:214] <- lapply(time_stamp_ICD10[2:214], as.Date) 

# time stamp for ICD9 - keep female cases only
time_stamp_ICD9 <- read.table(here(project_dir, "f41281.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid) %>%
  filter(., FID %in% X593_females_cases$FID)
time_stamp_ICD9[2:48] <- lapply(time_stamp_ICD9[2:48], as.Date)

# Main ----------------------------------------------------

## Time stamps for ICD10 and ICD9 - get most recent time stamp per participant
time_stamp_ICD10 <- time_stamp_ICD10 %>%
  rowwise() %>% mutate(max = max(c_across(f.41280.0.0:f.41280.0.212), na.rm = TRUE)) %>%
  separate(max, c("year_time_stamp_ICD10", "month", "day"), sep = "-") %>%
  dplyr::select(FID, year_time_stamp_ICD10) %>% drop_na(year_time_stamp_ICD10)

time_stamp_ICD9 <- time_stamp_ICD9 %>%
  rowwise() %>% mutate(max = max(c_across(f.41281.0.0:f.41281.0.46), na.rm = TRUE)) %>%
  separate(max, c("year_time_stamp_ICD9", "month", "day"), sep = "-") %>%
  dplyr::select(FID, year_time_stamp_ICD9) %>% drop_na(year_time_stamp_ICD9)

time_stamp_ICD <- inner_join(time_stamp_ICD10, time_stamp_ICD9, by = "FID") %>%
  mutate(year_time_stamp_ICD = case_when(
    year_time_stamp_ICD9 > year_time_stamp_ICD10 ~ year_time_stamp_ICD9,
    year_time_stamp_ICD9 <= year_time_stamp_ICD10 ~ year_time_stamp_ICD10
  )) %>%
  dplyr::select(FID, year_time_stamp_ICD)

## Age at menopause - 4 instances max. registered per participant - get the most recent one
age_menopause_spread <- age_menopause %>%
  tidyr::gather(instance, age_menopause, 2:5) %>%
  filter(., age_menopause > 1) %>% # removes values of -1 and -3
  tidyr::spread(., instance, age_menopause)

age_men_inst4 <- filter(age_menopause_spread, !is.na(X3581.3.0)) %>%
  dplyr::select(FID, age_menopause = X3581.3.0)
  
age_men_inst3 <- filter(age_menopause_spread, !is.na(X3581.2.0)) %>%
  filter(., !(FID %in% age_men_inst4$FID)) %>%
  dplyr::select(FID, age_menopause = X3581.2.0)

age_men_inst2 <- filter(age_menopause_spread, !is.na(X3581.1.0)) %>%
  filter(., !(FID %in% c(age_men_inst4$FID, age_men_inst3$FID))) %>%
  dplyr::select(FID, age_menopause = X3581.1.0)

age_men_inst1 <- filter(age_menopause_spread, !is.na(X3581.0.0)) %>%
  filter(., !(FID %in% c(age_men_inst4$FID, age_men_inst3$FID, age_men_inst2$FID))) %>%
  dplyr::select(FID, age_menopause = X3581.0.0)

age_menopause_latest <- rbind(age_men_inst1, age_men_inst2, age_men_inst3, age_men_inst4)

# Merge with hematuria status - get age at most recent time stamp and menopause status
all <- X593_females %>%
  left_join(., time_stamp_ICD, by = "FID") %>%
  mutate(age_at_time_stamp = as.numeric(year_time_stamp_ICD) - as.numeric(f.34.0.0)) %>%
  left_join(., age_menopause_latest) %>%
  mutate(menopause_status =
    case_when(
      is.na(age_menopause) ~ 0,
      !is.na(age_menopause) ~ 1,
      !is.na(age_at_time_stamp) & age_at_time_stamp < age_menopause ~ 0,
      !is.na(age_at_time_stamp) & age_at_time_stamp >= age_menopause ~ 1
    ))

filter(all, menopause_status == 0) %>%
  dplyr::select(FID, IID) %>%
  write.table(here(project_dir, "data_for_regenie", "ukb_wB_females_menopause_G1_ids.txt"), sep = "\t", row.names = F, quote = F)

filter(all, menopause_status == 1) %>%
  dplyr::select(FID, IID) %>%
  write.table(here(project_dir, "data_for_regenie", "ukb_wB_females_menopause_G2_ids.txt"), sep = "\t", row.names = F, quote = F)

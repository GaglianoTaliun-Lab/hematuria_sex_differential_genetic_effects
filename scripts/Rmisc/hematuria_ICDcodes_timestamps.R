# following extracting the ICD9 and ICD10 codes for hematuria in UKB RAP, as well as the time stamps,
# this script extracts the time stamps exclusively for the hematuria ICD codes, and creates a table
# with the oldest time stamp, representing the year of the first hematuria diagnosis for each case.

library(tidyverse)
library(here)
library(data.table)

project_dir="/Users/frida/Documents/research-projects/col4a2_hematuria"

# Get the number of instances - this will correspond to the number of columns for each instance

############ Load files:
hematuria_cases <- read.table(here(project_dir, "case.control-wb.withcovariates-4regenie.txt"), 
                              sep = "\t", header = T) %>%
  dplyr::filter(X593 == 1)

ICD10 <- read.table(here(project_dir, "hematuria_ICD10_timestamp.txt"), sep = "\t", header = T) %>%
  dplyr::filter(., eid %in% hematuria_cases$IID) %>%
  dplyr::select(-hematuria)

ICD9 <- read.table(here(project_dir, "hematuria_ICD9_timestamp.txt"), sep = "\t", header = T) %>%
  dplyr::filter(., eid %in% hematuria_cases$IID) %>%
  dplyr::select(-hematuria)

############ For ICD10:

max_instances_ICD10  <- ICD10$p41270 %>%
  replace_na("") %>% 
  str_count("\\,") %>% # Count the number of commas
  max(na.rm = TRUE) + 1 

# Perform the separation based on the comma operator
icd10_instances <- ICD10 %>%
  separate(
    p41270, 
    into = paste0("ICD10_a", 0:(max_instances_ICD10 - 1)), 
    sep = "\\,", 
    fill = "right"
  ) 

############ Extract timestamps for specific ICD codes (R31, R31.9, R31.0):

# Function to extract timestamps for specific ICD codes
extract_timestamps_for_codes <- function(icd_row, timestamp_row, target_codes) {
  # Get ICD code columns (ICD10_a0, ICD10_a1, etc.)
  icd_cols <- grep("^ICD10_a", names(icd_row), value = TRUE)
  # Get timestamp columns (p41280_a0, p41280_a1, etc.)
  timestamp_cols <- grep("^p41280_a", names(timestamp_row), value = TRUE)

  # Remove brackets and trim whitespace from ICD codes for matching
  icd_values_cleaned <- unlist(icd_row[icd_cols]) %>%
    str_remove_all("[\\[\\]]") %>%
    str_trim()

  # Initialize result list
  timestamps <- list()
  # For each target code
  for (code in target_codes) {
    # Find positions where this code appears (after removing brackets)
    positions <- which(icd_values_cleaned == code)
    # Get corresponding timestamps
    if (length(positions) > 0) {
      # positions are 1-indexed, so subtract 1 to get the correct timestamp column
      # (ICD10_a0 corresponds to p41280_a0, etc.)
      ts_values <- timestamp_row[timestamp_cols[positions]]
      timestamps[[code]] <- unlist(ts_values, use.names = FALSE)
      } else {
        timestamps[[code]] <- NA
      }
  }
  return(timestamps)
}

# Define target ICD codes
target_codes <- c("R31", "R31.9", "R31.0")

# Apply the function to each row
icd10_with_timestamps <- icd10_instances %>%
  rowwise() %>%
  mutate(
    R31_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "R31")$R31),
    R31.9_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "R31.9")$`R31.9`),
    R31.0_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "R31.0")$`R31.0`)
    ) %>%
  ungroup()

X593_ICD10_timestamps <- icd10_with_timestamps %>%
  dplyr::select(IID = eid,
                R31_timestamps,
                R31.9_timestamps,
                R31.0_timestamps)

############ For ICD9:

max_instances_ICD9  <- ICD9$p41271 %>%
  replace_na("") %>% 
  str_count("\\,") %>% # Count the number of commas
  max(na.rm = TRUE) + 1 

# Perform the separation based on the comma operator
icd9_instances <- ICD9 %>%
  separate(
    p41271, 
    into = paste0("ICD9_a", 0:(max_instances_ICD9 - 1)), 
    sep = "\\,", 
    fill = "right"
  ) 

############ Extract timestamps for specific ICD codes (5997, 59970, 59971):

# Function to extract timestamps for specific ICD codes
extract_timestamps_for_codes <- function(icd_row, timestamp_row, target_codes) {
  # Get ICD code columns (ICD9_a0, ICD9_a1, etc.)
  icd_cols <- grep("^ICD9_a", names(icd_row), value = TRUE)
  # Get timestamp columns (p41281_a0, p41281_a1, etc.)
  timestamp_cols <- grep("^p41281_a", names(timestamp_row), value = TRUE)
  
  # Remove brackets and trim whitespace from ICD codes for matching
  icd_values_cleaned <- unlist(icd_row[icd_cols]) %>%
    str_remove_all("[\\[\\]]") %>%
    str_trim()
  
  # Initialize result list
  timestamps <- list()
  # For each target code
  for (code in target_codes) {
    # Find positions where this code appears (after removing brackets)
    positions <- which(icd_values_cleaned == code)
    # Get corresponding timestamps
    if (length(positions) > 0) {
      # positions are 1-indexed, so subtract 1 to get the correct timestamp column
      # (ICD9_a0 corresponds to p41281_a0, etc.)
      ts_values <- timestamp_row[timestamp_cols[positions]]
      timestamps[[code]] <- unlist(ts_values, use.names = FALSE)
    } else {
      timestamps[[code]] <- NA
    }
  }
  return(timestamps)
}

# Define target ICD codes
target_codes <- c("5997", "59970", "59971")

# Apply the function to each row
icd9_with_timestamps <- icd9_instances %>%
  rowwise() %>%
  mutate(
    x5997_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "5997")$`5997`),
    x59970_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "59970")$`59970`),
    x59971_timestamps = list(extract_timestamps_for_codes(pick(everything()), pick(everything()), "59971")$`59971`)
  ) %>%
  ungroup()

X593_ICD9_timestamps <- icd9_with_timestamps %>%
  dplyr::select(IID = eid,
                x5997_timestamps,
                x59970_timestamps,
                x59971_timestamps)

# join ICD9 and ICD10 time stamps:
X593_timestamps <- full_join(X593_ICD10_timestamps, X593_ICD9_timestamps, by = "IID")

# convert to date format:
X593_timestamps[2:7] <- lapply(X593_timestamps[2:7], function(col) {
  lapply(col, as.Date)
})

# get the oldest and most recent:
X593_timestamps2 <- X593_timestamps %>%
  rowwise() %>% mutate(max = max(c_across(R31_timestamps:x59971_timestamps), na.rm = TRUE),
                       min = min(c_across(R31_timestamps:x59971_timestamps), na.rm = TRUE)) %>%
  separate(max, c("max_year_time_stamp", "month_max", "day_max"), sep = "-") %>%
  separate(min, c("min_year_time_stamp", "month_min", "day_mmin"), sep = "-") %>%
  drop_na(max_year_time_stamp) %>%
  drop_na(min_year_time_stamp)

X593_out <- X593_timestamps2 %>% dplyr::select(IID, min_year_time_stamp)

age_at_X593 <- left_join(hematuria_cases, X593_out, by = "IID") %>%
  mutate(age_at_X593_diagnosis = as.numeric(min_year_time_stamp) - as.numeric(f.34.0.0)) %>%
  dplyr::select(IID, age_at_X593_diagnosis)

write.table(age_at_X593, here(project_dir, "age_sex_analyses", "age_at_diagnosis_cases_X593_updated.txt"), sep = "\t", row.names = F, quote = F)








# Estradiol in the UKB, includes the following columns:
# X593 = hematuria cases (1) and controls (0)
# f.22001.0.0 = sex (all are females (0))
# f.34.0.0 = year of birth 
# f.22009.0.1 to f.22009.0.10 = PCs 1-10
# f.30800.0.0 = estradiol levels
# f.30801.year = year of estradiol measurement
# age_at_30800 = age of participant when estradiol measurement was taken
# f.30806.bin = estradiol detected (1) or not detected (0)
# f.3581 = reported age at menopause (multiple instances, so kept the most recent one)
# f.2724.final = menopause status (has had menopause = 1; no menopause = 0)
# year_time_stamp_ICD = ICD9/ICD10 time stamp for hematuria cases
# age_at_X593 = age of participant in their most recent ICD9/ICD10 time stamp
# menopause_status_estradiol = menopause status (based on estradiol date)
# menopause_status_X593 = menopause status (based on hematuria date)
# logeGFR2009 = estimated glomerular filtration rate (eGFR) 2009 equation
# uACR = albumin/creatinine ratio in urine (uACR)
# columns 26-35 (named with rsID) = minor allele counts/dosages

# Note 1: based the age at estradiol assessment only by substracting: year of assessment - year of birth
# (i.e., assuming everyone's birthday is the first day of the year).

# Note 2: NAs for menopause status are either those who preferred not to answer, responded 'not sure', or true NAs.

###################################################################################################
# Load Packages -----------------------------
###################################################################################################

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)

###################################################################################################
# Arguments ---------------------------------
###################################################################################################

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

###################################################################################################
# Read tables ---------------------------------
###################################################################################################

# X593 phenotype
hematuria_table <- read.table(here(project_dir, "regenie", "data_for_regenie", "case.control-wb.withcovariates-4regenie.txt"), sep = "\t", header = TRUE)
hematuria_cases <- filter(hematuria_table, X593 == 1)

# wB UKB females
females_wb <- read.table(here(project_dir, "regenie", "data_for_regenie", "ukb_WB_female_ids.txt"), sep = "\t", header = FALSE) %>%
  dplyr::rename(FID = V1, IID = V2)

# ukb withdraws
withdraws <- read.table(here(project_dir, "estradiol", "w66222_20241217.csv"), sep = "\t", header = FALSE)

# estradiol levels
estradiol <- read.table(here(project_dir, "estradiol", "f30800.txt"), sep = "\t", header = TRUE) %>%
  dplyr::rename(IID = f.eid) %>%
  mutate(f.30801.0.0 = as.Date(f.30801.0.0, '%Y-%m-%d'), f.30801.year = as.numeric(format(f.30801.0.0,'%Y'))) %>%
  dplyr::select(-f.30801.0.0)

# detection of estradiol
estradiol_det <- read.table(here(project_dir, "estradiol", "f30806.txt"), sep = "\t", header = TRUE) %>%
  dplyr::rename(IID = f.eid) %>%
  filter(., f.30806.0.0 %in% c(1,2,4)) %>%
  mutate(., f.30806.bin = case_when(
    f.30806.0.0 == 1 ~ "1",
    f.30806.0.0 %in% c(2,4) ~ "0"
  )) %>%
  dplyr::select(-f.30806.0.0)

# reported age at menopause
age_menopause <- read.table(here(project_dir, "f3581.tsv"), sep = "\t", header = T, fill = TRUE) %>%
  dplyr::rename(FID = eid)

#  have had menopause (yes/no)
menopause_status <- read.table(here(project_dir, "estradiol", "f2724_all_instances.txt"), sep = "\t", header = T, fill = TRUE) %>%
  dplyr::rename(FID = f.eid)

# age at hematuria diagnosis
age_diagnosis <- read.table(here(project_dir, "age_at_diagnosis_cases_X593_updated.txt"), sep = "\t", header = T) %>%
  dplyr::select(FID = IID, age_at_X593_diagnosis)

# eGFR-creatinine based (2009 equation) collected at baseline
egfr_all <- read.table(here(project_dir, "egfr", "egfr_withcovariates_for_REGENIE_updated.txt"), sep = "\t", header = TRUE) 

egfr <- egfr_all %>%
  dplyr::select(FID, logeGFR2009, wins_eGFR2009)

# calculate eGFR-cystatin C based (mg/L)
# 1) 2012 equation [https://www.kidney.org/ckd-epi-cystatin-c-equation-2012]
# 2) 100/cystatin C (as Moumita proposed)
cystatin_C <- read.table(here(project_dir, "f30720.txt"), sep = "\t", header = TRUE) %>%
  dplyr::rename(FID = f.eid) %>%
  left_join(., egfr_all, by = "FID") %>%
  mutate(eGFR_Moumita = 100/f.30720.0.0,
         eGFR_2012eq = case_when(
           f.22001.0.0 == 0 ~ (133 * pmin(f.30720.0.0/0.8, 1)^-0.499 * pmax(f.30720.0.0/0.8, 1)^-1.328 * (0.996)^f.21003.0.0 * 0.932),
           f.22001.0.0 == 1 ~ (133 * pmin(f.30720.0.0/0.8, 1)^-0.499 * pmax(f.30720.0.0/0.8, 1)^-1.328 * (0.996)^f.21003.0.0)
         ),
         wins_eGFR_2012eq = case_when(
           eGFR_2012eq <= 14.9999 ~ 15,
           eGFR_2012eq > 15.0001 ~ eGFR_2012eq
         ),
         ln_eGFR_2012eq = log(wins_eGFR_2012eq)) %>%
  dplyr::select(FID, eGFR_Moumita, wins_eGFR_2012eq, ln_eGFR_2012eq)

# 32,946 NA in Moumita's-based
# 113,694 NA in equation-based

# ---------------------------------------------------------------

p1 <- ggplot(cystatin_C, aes(x = eGFR_Moumita)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "eGFR (100/cystatin_C)", y = "Count", title = "") +
  geom_vline(xintercept = mean(cystatin_C$eGFR_Moumita, na.rm = TRUE), colour = "red") +
  annotate("text", label = str_c("N obs: ", filter(cystatin_C, !is.na(eGFR_Moumita)) %>% nrow()),
           x = 200, y = 40000, size = 4, colour = "black")

p2 <- ggplot(cystatin_C, aes(x = wins_eGFR_2012eq)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "eGFR (2012 equation)", y = "Count", title = "") +
  geom_vline(xintercept = mean(cystatin_C$wins_eGFR_2012eq, na.rm = TRUE), colour = "red") +
  annotate("text", label = str_c("N obs: ", filter(cystatin_C, !is.na(wins_eGFR_2012eq)) %>% nrow()),
           x = 130, y = 30000, size = 4, colour = "black")

# p1+p2

# ggsave(here(project_dir, "estradiol", "comparison_eGFRcys_UKB.jpg"))

# write.table(cystatin_C, here(project_dir, "estradiol", "eGFRcys_UKB.txt"), sep = "\t", row.names = F, quote = F)

# > summary(cystatin_C$eGFR_Moumita)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  13.36  101.94  112.87  113.28  124.38  338.98   32946 
# > summary(cystatin_C$eGFR_2012eq)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 5.06   77.03   88.83   87.93  100.44  164.13  113694

# ---------------------------------------------------------------

# uACR values computed as (albumin in mg/L)/(creatinine in μmol/L/1000) from urine samples collected at baseline.
# Albumin in numerator is UK Biobank field 30500 (“microalbumin in urine”), and individuals marked as having urinary 
# albumin < 6.7 mg/L (field 30505) should be assigned as having a urinary albumin value of 6.7 mg/L.
# Creatinine field for denominator: 30510. *note the units is μmol/L - so need to divide by 1000

# there are two instances for each, choose the first one, if it is NA then choose the second one, if available.

uacr <- read.table(here(project_dir, "estradiol", "ukb_phenotypes_for_uACR.txt"), sep = "\t", header = TRUE, na.strings = "NA", fill = TRUE) %>%
  mutate(f.30500 = case_when(
    !is.na(f.30500.0.0) ~ f.30500.0.0,
    is.na(f.30500.0.0) & !is.na(f.30505.0.0) ~ 6.7,
    is.na(f.30500.0.0) & is.na(f.30505.0.0) ~ NA
  ),
  f.30510 = case_when(
    !is.na(f.30510.0.0) ~ f.30510.0.0/1000,
    is.na(f.30510.0.0) ~ NA
  ),
  uACR = f.30500 / f.30510) %>%
  dplyr::select(FID = f.eid, uACR)

# genotypes for variants of interest
genotypes_files <- setNames(list.files(here(project_dir, "estradiol"), pattern = "*.raw", full.names = TRUE),
                            nm = list.files(here(project_dir, "estradiol"), pattern = "*.raw", full.names = FALSE) %>% 
                              stringr::str_extract("rs[:digit:]+"))

genotypes_list <- lapply(genotypes_files, function(x) read.table(x, sep = "\t", header = TRUE) %>%
                           dplyr::select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE))

###################################################################################################
# Main ----------------------------------------
###################################################################################################

## Age at menopause - 4 instances max. registered per participant - get the most recent one for each participant
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

## menopause status - merge all instances (values of 3 and 2 answered 'not sure', 1 answered they have had menopause and 0 they haven't had)
menopause_status <- menopause_status %>%
  mutate(f.2724.final = case_when(
    f.2724.0.0 == 3 | f.2724.1.0 == 3 | f.2724.2.0 == 3 | f.2724.3.0 == 3 ~ NA,
    f.2724.0.0 == 2 | f.2724.1.0 == 2 | f.2724.2.0 == 2 | f.2724.3.0 == 2 ~ NA,
    f.2724.0.0 == 1 | f.2724.1.0 == 1 | f.2724.2.0 == 1 | f.2724.3.0 == 1 ~ 1,
    f.2724.0.0 == 0 | f.2724.1.0 == 0 | f.2724.2.0 == 0 | f.2724.3.0 == 0 ~ 0,
    TRUE ~ NA
  )) %>%
  filter(., FID %in% females_wb$FID) %>%
  dplyr::select(FID, f.2724.final)

# Genotypes for variants of interest - join list into a data frame
genotypes_joint = genotypes_list[[1]]
for (i in 2:length(genotypes_list)){
  genotypes_joint <- inner_join(genotypes_joint, genotypes_list[[i]], by = "FID")
}
# of the 10 variants, the genotype count/dosage corresponds to the major allele, except for rs28807105 (checked manually), 
# so need to change the genotypes to get the minor allele counts/dosages
genotypes_joint[2:10] <- lapply(genotypes_joint[2:10], function(x) 2 - x) 

rsids <- colnames(genotypes_joint[2:11]) %>%
  stringr::str_extract("rs[:digit:]+")
rsids <- c("FID", rsids)
colnames(genotypes_joint) <- rsids

###################################################################################################
# Merge tables and further processing --------------------------------------------
###################################################################################################

# from the hematuria table, keep only wB females and remove participants who withdrew:
all <- filter(hematuria_table, IID %in% females_wb$IID) %>%
  filter(., !(IID %in% withdraws$V1)) %>%
  # include estradiol levels and estradiol detection (binary variable):
  left_join(., estradiol, by = "IID") %>%
  mutate(age_at_30800 = f.30801.year - f.34.0.0) %>%
  left_join(., estradiol_det, by = "IID") %>%
  # include age at menopause and menopause status:
  left_join(., age_menopause_latest, by = "FID") %>%
  left_join(., menopause_status, by = "FID") %>%
  # if participant reported age at menopause (f.3581), but their menopause status (f.2724.final) = 0, change menopause status to 1.
  # if participant reported age at menopause (f.3581) values of -1 (prefer not to answer) or -3 (do not know), change menopause status to NA.
  mutate(f.2724.final = case_when(
    f.3581 > 0 & f.2724.final == 0 ~ 1,
    f.3581 < 0 ~ NA,
    TRUE ~ f.2724.final
  )) %>%
  left_join(., age_diagnosis, by = "FID") %>%
    dplyr::rename(age_at_X593 = age_at_X593_diagnosis)

# fix menopause status based on the age at estradiol assessment (X3581.latest)
# in case the assessment was done before menopause, but participant has since reported age at menopause
all <- all %>%
  mutate(menopause_status_estradiol = case_when(
    is.na(f.2724.final) ~ NA,
    f.2724.final == 0 ~ 0,
    f.2724.final == 1 & (f.3581 <= age_at_30800) == TRUE ~ 1,
    f.2724.final == 1 & (f.3581 > age_at_30800) == TRUE ~ 0
  ),
  # fix menopause status based on the age at oldest ICD9/ICD10 time stamp (age_at_X593)
  # in case the diagnosis was done before menopause, but participant has since reported age at menopause
  menopause_status_X593 = case_when(
    is.na(f.2724.final) ~ NA,
    f.2724.final == 0 ~ 0,
    X593 == 0 ~ f.2724.final,
    f.2724.final == 1 & X593 == 1 & (f.3581 <= age_at_X593) == TRUE ~ 1,
    f.2724.final == 1 & X593 == 1 & (f.3581 > age_at_X593) == TRUE ~ 0
  )
  )

# add eGFR and uACR measurements:
all <- all %>%
  left_join(., egfr, by = "FID") %>%
  left_join(., uacr, by = "FID") %>%
  left_join(., cystatin_C, by = "FID")

# add genotypes (dosages)
all <- all %>%
  left_join(., genotypes_joint, by = "FID") %>%
  dplyr::select(FID,
                IID,
                X593,
                f.34.0.0,
                f.22009.0.1,
                f.22009.0.2,
                f.22009.0.3,
                f.22009.0.4,
                f.22009.0.5,
                f.22009.0.6,
                f.22009.0.7,
                f.22009.0.8,
                f.22009.0.9,
                f.22009.0.10,
                f.30800.0.0,
                f.30806.bin,
                menopause_status_estradiol,
                menopause_status_X593,
                logeGFR2009,
                eGFR2009 = wins_eGFR2009,
                uACR,
                eGFRcys_MB = eGFR_Moumita,
                eGFRcys_2012eq = wins_eGFR_2012eq,
                rs2145409,
                rs35997018,
                rs75770066,
                rs71181755,
                rs16991615,
                rs774021038,
                rs2454949,
                rs34929649,
                rs72823349,
                rs28807105
  )

final_table <- all %>%
  mutate(eGFR2009 = case_when(
    eGFR2009 > mean(all$eGFR2009, na.rm = T) + 5*sd(all$eGFR2009, na.rm = T) ~ NA,
    eGFR2009 <= mean(all$eGFR2009, na.rm = T) + 5*sd(all$eGFR2009, na.rm = T) ~ eGFR2009
  ),
  eGFRcys_MB = case_when(
    eGFRcys_MB > mean(all$eGFRcys_MB, na.rm = T) + 5*sd(all$eGFRcys_MB, na.rm = T) ~ NA,
    eGFRcys_MB <= mean(all$eGFRcys_MB, na.rm = T) + 5*sd(all$eGFRcys_MB, na.rm = T) ~ eGFRcys_MB
  ),
  eGFRcys_2012eq = case_when(
    eGFRcys_2012eq > mean(all$eGFRcys_2012eq, na.rm = T) + 5*sd(all$eGFRcys_2012eq, na.rm = T) ~ NA,
    eGFRcys_2012eq <= mean(all$eGFRcys_2012eq, na.rm = T) + 5*sd(all$eGFRcys_2012eq, na.rm = T) ~ eGFRcys_2012eq
  ))

###################################################################################################
# Write table --------------------------------------------
###################################################################################################

write.table(final_table, here(project_dir, "estradiol", "final_estradiol_table_clean_updated.txt"), sep = "\t", row.names = F, quote = F)

###################################################################################################
# Extra analyses
###################################################################################################

diffs <- all %>%
  mutate(diff_estradiol_X593 = case_when(
    menopause_status_estradiol == menopause_status_X593 ~ "+",
    menopause_status_estradiol != menopause_status_X593 ~ "-"
  ),
  diff_estradiol_f2724 = case_when(
    menopause_status_estradiol == f.2724.final ~ "+",
    menopause_status_estradiol != f.2724.final ~ "-"
  ),
  diff_X593_f2724 = case_when(
    menopause_status_X593 == f.2724.final ~ "+",
    menopause_status_X593 != f.2724.final ~ "-"
  )
  ) %>%
  dplyr::select(FID, X593, diff_estradiol_X593, diff_estradiol_f2724, diff_X593_f2724)

# There are 889 mismatching values of menopause status defined between estradiol and hematuria.

p1 <- ggplot(final_table, aes(x = eGFR2009)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "eGFR (creatinine)", y = "Count", title = "")

p2 <- ggplot(final_table, aes(x = eGFRcys_MB)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "eGFR (cystatin MB)", y = "Count", title = "")

p3 <- ggplot(final_table, aes(x = eGFRcys_2012eq)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "eGFR (cystatin equation)", y = "Count", title = "")

p1 + p2 + p3
ggsave(here(project_dir, "estradiol", "comparison_eGFR_all_UKB.jpg"))

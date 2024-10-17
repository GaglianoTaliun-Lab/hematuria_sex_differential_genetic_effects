# obtain sample IDs stratified by age and sex for REGENIE

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
cbbPalette <- c("#CC79A7","#56B4E9", "#E69F00", "#D55E00", "#009E73", "#0072B2", "#F0E442")

# Main --------------------------------------

# males = 1
# females = 0
hem <- read.table(here(project_dir, "case.control-wb.withcovariates-4regenie.txt"), sep = "\t", header = T)

wb_males <- read.table(here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ids.txt"), sep = "\t", header = F) %>%
  dplyr::rename(FID = V1,
         IID = V2)
wb_females <- read.table(here(project_dir, "regenie", "data_for_regenie", "ukb_WB_female_ids.txt"), sep = "\t", header = F) %>%
  dplyr::rename(FID = V1,
         IID = V2)

wb <- rbind(wb_males, wb_females)

age <- read.table(here(project_dir, "f21022.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid, age = f.21022.0.0)

# self-reported sex (males = 1; females = 0)
gender <- read.table(here(project_dir, "f31.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid, gender = f.31.0.0)

hem_all <- inner_join(hem, wb, by = c("FID", "IID")) %>%
  left_join(., age, by = "FID") %>%
  left_join(., gender, by = "FID")

# ----------- Group by ages:

# males below 55
wb_males_G1 <- inner_join(wb_males, age, by = "FID") %>%
  filter(., age < 55) %>%
  dplyr::select(FID, IID)
# males above 55
wb_males_G2 <- inner_join(wb_males, age, by = "FID") %>%
  filter(., age >= 55) %>%
  dplyr::select(FID, IID)

# females below 55
wb_females_G1 <- inner_join(wb_females, age, by = "FID") %>%
  filter(., age < 55) %>%
  dplyr::select(FID, IID)
# females above 55
wb_females_G2 <- inner_join(wb_females, age, by = "FID") %>%
  filter(., age >= 55) %>%
  dplyr::select(FID, IID)

# write files for REGENIE --keep:
write.table(wb_males_G1, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ageG1_ids.txt"), sep = "\t", row.names = F, quote = F)
write.table(wb_males_G2, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ageG2_ids.txt"), sep = "\t", row.names = F, quote = F)
write.table(wb_females_G1, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_female_ageG1_ids.txt"), sep = "\t", row.names = F, quote = F)
write.table(wb_females_G2, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_female_ageG2_ids.txt"), sep = "\t", row.names = F, quote = F)

# join with hematuria data:
hem_wb_males_G1 <- right_join(hem, wb_males_G1, by = c("FID", "IID"))
hem_wb_males_G2 <- right_join(hem, wb_males_G2, by = c("FID", "IID"))
hem_wb_females_G1 <- right_join(hem, wb_females_G1, by = c("FID", "IID"))
hem_wb_females_G2 <- right_join(hem, wb_females_G2, by = c("FID", "IID"))

# write.table(hem_wb_males_G1, here(project_dir, "regenie", "data_for_regenie", "case.control-X593-males_wb_G1.withcovariates.txt"), sep = "\t", row.names = F, quote = F)
# write.table(hem_wb_males_G2, here(project_dir, "regenie", "data_for_regenie", "case.control-X593-males_wb_G2.withcovariates.txt"), sep = "\t", row.names = F, quote = F)
# write.table(hem_wb_females_G1, here(project_dir, "regenie", "data_for_regenie", "case.control-X593-females_wb_G1.withcovariates.txt"), sep = "\t", row.names = F, quote = F)
# write.table(hem_wb_females_G2, here(project_dir, "regenie", "data_for_regenie", "case.control-X593-females_wb_G2.withcovariates.txt"), sep = "\t", row.names = F, quote = F)

# -------- Include age at menopause:

age_menopause <- read.table(here(project_dir, "f3581.tsv"), sep = "\t", header = T) %>%
  dplyr::rename(FID = eid) %>%
  inner_join(hem_all, age_menopause, by = "FID") %>%
  dplyr::select(FID, f.22001.0.0, X3581.0.0, X3581.1.0, X3581.2.0, X3581.3.0, X593, gender)

# 19 males who reported age at menopause (all had female sex self-reported)
# exclude them
males <- age_menopause %>% 
  filter(f.22001.0.0 == 1) %>% filter(!is.na(X3581.0.0) | !is.na(X3581.1.0) | !is.na(X3581.2.0) | !is.na(X3581.3.0)) %>%
  dplyr::select(FID, f.22001.0.0, X3581.0.0, X3581.1.0, X3581.2.0, X3581.3.0, X593, gender)

age_menopause_females <- age_menopause %>%
  filter(f.22001.0.0 == 0)

age_menopause_long <- age_menopause_females %>%
  tidyr::gather(instance, age_menopause, 3:6) %>%
  filter(., !is.na(age_menopause)) %>% # removes unreported instances
  filter(., age_menopause > 1) # removes values of -1 and -3

duplicate_counts <- age_menopause_long %>%
  add_count(FID) %>%
  filter(n > 1) %>%
  distinct() %>%
  group_by(FID) %>%
  mutate(diff = abs(age_menopause-lag(age_menopause, default = first(age_menopause))))

# participants with only one instance recorded (N = 118,905):
one_instance_counts <- age_menopause_long %>%
  add_count(FID) %>%
  filter(n == 1)

ggplot(duplicate_counts, aes(x = diff)) +
  geom_bar() +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5, size = 2) +
  labs( x = "Absolute difference in reported age at menopause across instances", y = "Frequency") +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank())
ggsave(here(project_dir, "age_sex_analyses", "frequency_differences_age_at_menopause_across_instances.png"), width = 15, height = 15, units = "cm")

# --------- Include time stamps for ICD codes

time_stamp_ICD10 <- read.table(here(project_dir, "f41280.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid) %>%
  filter(., FID %in% wb_females$FID)
time_stamp_ICD10[2:214] <- lapply(time_stamp_ICD10[2:214], as.Date)

time_stamp_ICD9 <- read.table(here(project_dir, "f41281.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid) %>%
  filter(., FID %in% wb_females$FID)
time_stamp_ICD9[2:48] <- lapply(time_stamp_ICD9[2:48], as.Date)

# assuming the most recent time stamp corresponds to their hematuria status, get most recent time stamp:
time_stamp_ICD10_max <- time_stamp_ICD10 %>%
  rowwise() %>% mutate(max = max(c_across(f.41280.0.0:f.41280.0.212), na.rm = TRUE)) %>%
  separate(max, c("year_time_stamp_ICD10", "month", "day"), sep = "-") %>%
  dplyr::select(FID, year_time_stamp_ICD10) %>% drop_na(year_time_stamp_ICD10)

time_stamp_ICD9_max <- time_stamp_ICD9 %>%
  rowwise() %>% mutate(max = max(c_across(f.41281.0.0:f.41281.0.46), na.rm = TRUE)) %>%
  separate(max, c("year_time_stamp_ICD9", "month", "day"), sep = "-") %>%
  dplyr::select(FID, year_time_stamp_ICD9) %>% drop_na(year_time_stamp_ICD9)

time_stamp_ICD <- inner_join(time_stamp_ICD10, time_stamp_ICD9, by = "FID") %>%
  rowwise() %>% mutate(year_time_stamp_ICD = max(c_across(year_time_stamp_ICD9, year_time_stamp_ICD10)))

# join with YOB, and get the participants' age at the most recent time stamp
yob <- hem %>% filter(., FID %in% wb_females$FID) %>%
  dplyr::select(FID, f.34.0.0)
left_join(time_stamp_ICD, yob, by = "FID") %>%
  mutate(age_at_time_stamp = year_time_stamp_ICD - f.34.0.0)

# separate either below or above 50 (mean age at menopause) - or either with age at menopause

# --------- Plots

ggplot(hem_full_wb, aes(x = age, fill = as.factor(sex))) +
  geom_histogram(binwidth = 1, colour = "black") +
  labs(x="Age (years)", y = "Count") +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust =0.5),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y  = element_blank()
  ) +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  geom_vline(aes(xintercept=55),
             color="blue", linetype="dashed", size=1) +
  guides(fill=guide_legend(title = "Group"))
ggsave(here(project_dir, "age_sex_analyses", "histogram_wB_age_bysex_vline55.png"), width=10, height=8)


hem_f_wb <- hem_full_wb %>%
  filter(., f.22001.0.0 == 0)
hem_m_wb <- hem_full_wb %>%
  filter(., f.22001.0.0 == 1)

ggplot(hem_f_wb, aes(x = age, fill = as.factor(X593))) +
  geom_bar(colour = "black", position = "fill") +
  labs(x="Age (years)", y = "Proportion") +
  scale_fill_manual(values=cbbPalette) +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust =0.5),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y  = element_blank()
  ) +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  geom_vline(aes(xintercept=55),
             color="blue", linetype="dashed", size=1) +
  guides(fill=guide_legend(title = "Hematuria (females)"))
ggsave(here(project_dir, "age_sex_analyses", "histogram_wB_females_age_case_proportion.png"), width=10, height=8)

ggplot(hem_m_wb, aes(x = age, fill = as.factor(X593))) +
  geom_bar(colour = "black", position = "fill") +
  labs(x="Age (years)", y = "Proportion") +
  scale_fill_manual(values=cbbPalette) +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust =0.5),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.text.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y  = element_blank()
  ) +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  geom_vline(aes(xintercept=55),
             color="blue", linetype="dashed", size=1) +
  guides(fill=guide_legend(title = "Hematuria (males)"))
ggsave(here(project_dir, "age_sex_analyses", "histogram_wB_males_age_case_proportion.png"), width=10, height=8)

# --------- Logistic regression

model <- glm(X593 ~ age + f.22001.0.0 + age*f.22001.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = hem_all)
summary(model)
sink(here(project_dir, "age_sex_analyses", "logistic_model_sex_interaction.txt"))
summary(model)
sink()

model <- glm(X593 ~ age + gender + age*gender + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = hem_all)
summary(model)
sink(here(project_dir, "age_sex_analyses", "logistic_model_gender_interaction.txt"))
summary(model)
sink()

model <- glm(X593 ~ age + gender + age*gender + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = not_conc)
summary(model)
sink(here(project_dir, "age_sex_analyses", "logistic_model_gender_interaction.txt"))
summary(model)
sink()
model <- glm(X593 ~ age + f.22001.0.0 + age*f.22001.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = not_conc)
summary(model)
sink(here(project_dir, "age_sex_analyses", "logistic_model_sex_interaction.txt"))
summary(model)
sink()


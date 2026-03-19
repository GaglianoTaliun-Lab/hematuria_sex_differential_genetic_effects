# obtain sample IDs stratified by age and/or menopause status and/or sex for REGENIE

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

# females
hem_f <- hem_all %>%
  dplyr::filter(., f.22001.0.0 == 0)
# males
hem_m <- hem_all %>%
  dplyr::filter(., f.22001.0.0 == 1)

# get mean and SE for age by sex:
mean(hem_f$age)
sd(hem_f$age)
mean(hem_m$age)
sd(hem_m$age)

# ----------- Group by ages (used for COGWI male stratified analysis - testing only two significant signals in chr 16 and X):

# males below 50
wb_males_G1 <- inner_join(wb_males, age, by = "FID") %>%
  filter(., age < 50) %>%
  dplyr::select(FID, IID)

# males equal or above 50
wb_males_G2 <- inner_join(wb_males, age, by = "FID") %>%
  filter(., age >= 50) %>%
  dplyr::select(FID, IID)

# sex-combined below 50 and keep cases only
wb_sexcombined_G1 <- hem_all %>%
  filter(., age < 50) %>%
  filter(., X593 == 1) %>%
  dplyr::select(FID, IID)

# sex-combined equal or above 50 and keep cases only
wb_sexcombined_G2 <- hem_all %>%
  filter(., age >= 50) %>%
  filter(., X593 == 1) %>%
  dplyr::select(FID, IID)

# write files for REGENIE --keep:
# write.table(wb_males_G1, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ageG1_ids.txt"), sep = "\t", row.names = F, quote = F)
# write.table(wb_males_G2, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ageG2_ids.txt"), sep = "\t", row.names = F, quote = F)
# write.table(wb_sexcombined_G1, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_sexcombined_cases_ageG1_ids.txt"), sep = "\t", row.names = F, quote = F)
# write.table(wb_sexcombined_G2, here(project_dir, "regenie", "data_for_regenie", "ukb_WB_sexcombined_cases_ageG2_ids.txt"), sep = "\t", row.names = F, quote = F)


# --------- Plots

hem_all_plot <- hem_all %>%
  mutate(sex = case_when(
    f.22001.0.0 == 1 ~ "male",
    f.22001.0.0 == 0 ~ "female"
  ))

ggplot(hem_all_plot, aes(x = age, fill = as.factor(sex))) +
  geom_histogram(binwidth = 1, colour = "black") +
  labs(x="Age (years)", y = "Count") +
  theme(legend.title = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust =0.5),
        axis.line.x = element_line(linewidth = 0.6),
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
  geom_vline(aes(xintercept=50),
             color="red", linetype="dashed", linewidth=1) +
  scale_fill_manual(values = c("female" = "#DDAA33", "male" = "#BB5566"),
                     name = "Sex",
                     labels = c("female" = "Females", "male" = "Males"))

ggsave(here(project_dir, "age_sex_analyses", "histogram_wB_age_bysex_vline50.png"), width=10, height=8, dpi = 300)


hem_f_wb <- hem_all %>%
  filter(., f.22001.0.0 == 0)
hem_m_wb <- hem_all %>%
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
  geom_vline(aes(xintercept=50),
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
  geom_vline(aes(xintercept=50),
             color="blue", linetype="dashed", size=1) +
  guides(fill=guide_legend(title = "Hematuria (males)"))
ggsave(here(project_dir, "age_sex_analyses", "histogram_wB_males_age_case_proportion.png"), width=10, height=8)

# --------- Logistic regression tests for age*sex interactions:

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

# model <- glm(X593 ~ age + gender + age*gender + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = not_conc)
# summary(model)
# sink(here(project_dir, "age_sex_analyses", "logistic_model_gender_interaction.txt"))
# summary(model)
# sink()
# model <- glm(X593 ~ age + f.22001.0.0 + age*f.22001.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10, family=binomial(link='logit'), data = not_conc)
# summary(model)
# sink(here(project_dir, "age_sex_analyses", "logistic_model_sex_interaction.txt"))
# summary(model)
# sink()


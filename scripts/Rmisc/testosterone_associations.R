# Description - association between testosterone levels in the UKB and hematuria.
# testosterone units are nmol/L and levels were assayed with an immunoassay 
# There are three extracted variables:
### 1) assay date: f30851
### 2) levels: f30850
### 3) detectable (binary): f30856

###################################################################################################
# Load Packages -----------------------------
###################################################################################################

library(here)
library(stringr)
library(tidyverse)
library(data.table)
library(patchwork)
library(grid)
library(gridExtra)
library(forcats)
library(vctrs)
library(AER)

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

# wB UKB males
males_wb <- read.table(here(project_dir, "regenie", "data_for_regenie", "ukb_WB_male_ids.txt"), sep = "\t", header = FALSE) %>%
  dplyr::rename(FID = V1, IID = V2)

# ukb withdraws
withdraws <- read.table(here(project_dir, "w66222_20250818.csv"), sep = "\t", header = FALSE) %>% dplyr::rename(IID = V1)

# age at hematuria diagnosis
age_diagnosis <- read.table(here(project_dir, "age_at_diagnosis_cases_X593_updated.txt"), sep = "\t", header = T) %>%
  dplyr::select(IID, age_at_X593_diagnosis)

# testosterone levels
testosterone <- read.table(here(project_dir, "testosterone", "testosterone_fields_RAP.txt"), sep = "\t", header = TRUE) %>%
  dplyr::rename(IID = eid) %>%
  mutate(f30851 = as.Date(p30851_i0, '%Y-%m-%d'), f30851.year = as.numeric(format(f30851,'%Y')),
         f30856 = case_when(
           p30856_i0 == 1 ~ 1,
           p30856_i0 == 2 ~ 0,
           p30856_i0 == 3 ~ NA,
           p30856_i0 == 4 ~ 0,
           p30856_i0 == 5 ~ NA
         ),
         f30850_tobit = case_when(
           f30856 == 1 ~ p30850_i0,
           f30856 == 0 ~ 0.34,
           is.na(f30856) ~ NA
         )) %>%
  dplyr::select(IID,
                f30850 = p30850_i0,
                f30851.year, 
                f30856,
                f30850_tobit)

# merge tables (females == 0; males == 1):
all <- left_join(hematuria_table, age_diagnosis, by = "IID") %>%
  left_join(., testosterone, by  = "IID") %>%
  filter(., !(IID %in% withdraws$IID))


###################################################################################################
# Descriptive plots ---------------------------------
###################################################################################################

testos_f <- ggplot(all %>% dplyr::filter(., f.22001.0.0 == 0) %>% dplyr::filter(., !is.na(f30850)), aes(x = f30850)) +
  geom_histogram(binwidth = 1, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "testosterone in females (nmol/L)", y = "Count (N = 175,151)", title = "")

testos_m <- ggplot(all %>% dplyr::filter(., f.22001.0.0 == 1) %>% dplyr::filter(., !is.na(f30850)), aes(x = f30850)) +
  geom_histogram(binwidth = 1, colour = "black", fill = "white") +
  theme_classic(base_size=14) +
  labs(x= "testosterone in males (nmol/L)", y = "Count (N = 177,267)", title = "")

testos_all <- ggplot(all %>% dplyr::filter(., !is.na(f30850)), aes(x = f30850, fill=as.factor(f.22001.0.0))) +
  geom_histogram(binwidth = 1, colour = "black", alpha = 0.4, position = "identity") +
  theme_classic(base_size=14) +
  scale_fill_manual(values=c("#69b3a2", "#404080"), labels = c("Females", "Males"), breaks = c(0, 1)) +
  guides(fill = guide_legend(title = "Sex")) +
  labs(x= "testosterone levels (nmol/L)", y = "Count (N = 352,418)", title = "")

ggsave(here(project_dir, "testosterone", "testosterone_distribution_f30850.jpg"),
       plot = testos_all,
       width = 6, height = 7, dpi = 300)

figure_all <- grid.arrange(testos_f, testos_m, nrow = 1, widths = c(5, 5))

ggsave(here(project_dir, "testosterone", "testosterone_distribution_f30850_by_sex.jpg"),
       plot = figure_all,
       width = 15, height = 6, dpi = 300)

###################################################################################################
# Associations with hematuria ---------------------------------
###################################################################################################

# ----------- quantitative testosterone using a tobit regression

# model the effect of quantitative testosterone on hematuria (including an interaction effect of hematuria*sex)
# Note: given that it has a bimodal distribution, it is better to separate by sex and use ln for females.
tobit_model <- tobit(formula = f30850_tobit ~ X593 + f.22001.0.0 + f.34.0.0 + (X593*f.22001.0.0) + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                     left = 0.34, right = Inf,
                     data = all)
sink(here(project_dir, "testosterone", "results_tobit_model_X593_all.txt"))
summary(tobit_model)
sink()

# stratify by sex

# males:
tobit_model_m <- tobit(formula = f30850_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                     left = 0.34, right = Inf,
                     data = all %>% dplyr::filter(., f.22001.0.0 == 1))
sink(here(project_dir, "testosterone", "results_tobit_model_X593_males.txt"))
summary(tobit_model_m)
sink()

# females:
tobit_model_f <- tobit(formula = log(f30850_tobit) ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                       left = 0.34, right = Inf,
                       data = all %>% dplyr::filter(., f.22001.0.0 == 0))
sink(here(project_dir, "testosterone", "results_log_tobit_model_X593_females.txt"))
summary(tobit_model_f)
sink()

# stratify by sex and by age groups
all <- dplyr::filter(., !is.na(f30851.year)) %>%
  mutate(age_at_testost_assay = f30851.year - f.34.0.0,
    age_groups_testost_assay = case_when(
      age_at_testost_assay < 50 ~ 1,
      age_at_testost_assay >= 50 ~ 2,
    is.na(age_at_testost_assay) ~ NA
  ))

# perform tobit test by sex and by age group:

# males
## group 1 (< 50 years):
tobit_model_m1 <- tobit(formula = f30850_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                       left = 0.34, right = Inf,
                       data = all %>% dplyr::filter(., f.22001.0.0 == 1 & age_groups_testost_assay == 1))
sink(here(project_dir, "testosterone", "results_tobit_model_X593_males_group1.txt"))
summary(tobit_model_m1)
sink()

# males
## group 2 (>= 50 years):
tobit_model_m2 <- tobit(formula = f30850_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                        left = 0.34, right = Inf,
                        data = all %>% dplyr::filter(., f.22001.0.0 == 1 & age_groups_testost_assay == 2))
sink(here(project_dir, "testosterone", "results_tobit_model_X593_males_group2.txt"))
summary(tobit_model_m2)
sink()

# females
## group 1 (< 50 years):
tobit_model_f1 <- tobit(formula = log(f30850_tobit) ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                       left = 0.34, right = Inf,
                       data = all %>% dplyr::filter(., f.22001.0.0 == 0 & age_groups_testost_assay == 1))
sink(here(project_dir, "testosterone", "results_log_tobit_model_X593_females_group1.txt"))
summary(tobit_model_f1)
sink()

# females
## group 2 (>= 50 years):
tobit_model_f2 <- tobit(formula = log(f30850_tobit) ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                        left = 0.34, right = Inf,
                        data = all %>% dplyr::filter(., f.22001.0.0 == 0 & age_groups_testost_assay == 2))
sink(here(project_dir, "testosterone", "results_log_tobit_model_X593_females_group2.txt"))
summary(tobit_model_f2)
sink()

# ----------- binary testosterone using a logistic regression

count(all, f.22001.0.0, f30856)

# Females w/ detected levels: 175,151
# Females w/out detected levels: 33,115
# Males w/ detected levels: 177,267
# Males w/out detected levels: 205

model <- glm(X593 ~ f30856 + f.22001.0.0 + f.34.0.0 + (f30856*f.22001.0.0) + f.22009.0.1 + f.22009.0.2 + 
               f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
               f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
             data = all, family=binomial(link='logit'))
sink(here(project_dir, "testosterone", "results_binary_model_X593_all.txt"))
summary(model)
sink()

model_m <- glm(X593 ~ f30856 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + 
               f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
               f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
             data = all %>% dplyr::filter(., f.22001.0.0 == 1), family=binomial(link='logit'))
sink(here(project_dir, "testosterone", "results_binary_model_X593_males.txt"))
summary(model_m)
sink()

model_f <- glm(X593 ~ f30856 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + 
                 f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
                 f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
               data = all %>% dplyr::filter(., f.22001.0.0 == 0), family=binomial(link='logit'))
sink(here(project_dir, "testosterone", "results_binary_model_X593_females.txt"))
summary(model_f)
sink()

# stratify ALL by age groups:
model_g1 <- glm(X593 ~ f30856 + f.22001.0.0 + f.34.0.0 + (f30856*f.22001.0.0) + f.22009.0.1 + f.22009.0.2 + 
               f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
               f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
             data = all %>% dplyr::filter(., age_groups_testost_assay == 1), family=binomial(link='logit'))
sink(here(project_dir, "testosterone", "results_binary_model_X593_all_group1.txt"))
summary(model_g1)
sink()

model_g2 <- glm(X593 ~ f30856 + f.22001.0.0 + f.34.0.0 + (f30856*f.22001.0.0) + f.22009.0.1 + f.22009.0.2 + 
                  f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
                  f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                data = all %>% dplyr::filter(., age_groups_testost_assay == 2), family=binomial(link='logit'))
sink(here(project_dir, "testosterone", "results_binary_model_X593_all_group2.txt"))
summary(model_g2)
sink()

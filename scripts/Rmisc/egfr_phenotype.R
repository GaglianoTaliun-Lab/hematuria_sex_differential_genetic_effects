# Description: compute eGFR phenotype for Regenie in UKB RAP

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

egfr <- read.table(here(project_dir, "egfr", "egfr-withcovariates.txt"), sep = " ", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute eGFR phenotypes

# creatinine: f.30700.0.0 (in umol/L)
# cystatin C: f.30720.0.0 (in mg/L) - need to obtain the values

# sex: f.22001.0.0, where 0 = females and 1 = males

# age at assessment day: f.21003.0.0

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. compute eGFR(creatinine) based on the 2021 equation (uses mg/dL):
# https://www.kidney.org/content/ckd-epi-creatinine-equation-2021
# to convert from umol/L to mg/dL multiply by 0.0113 (source: Nephrology Dialysis Transplantation 2004, DOI: 10.1093/ndt/gfh1030, Suppl 2)
# and compute the natural logarithm of winsorized eGFR creatinine (2009) as in Stanzick et al. 2021 (See Methods of the paper)

k_F=0.7
k_M=0.9
alpha_F_2021=-0.241
alpha_M_2021=-0.302
alpha_F_2009=-0.329
alpha_M_2009=-0.411

egfr <- egfr %>%
  mutate(FID = IID) %>%
  mutate(f.30700.0.0_mgL = f.30700.0.0*0.0113,
         egfr_cr_2021 =
           case_when(
             f.22001.0.0 == 0 ~ round(142 * pmin(f.30700.0.0_mgL/k_F, 1)^alpha_F_2021 * pmax(f.30700.0.0_mgL/k_F, 1)^-1.2 * 0.9938^f.21003.0.0 * 1.012, 2),
             f.22001.0.0 == 1 ~ round(142 * pmin(f.30700.0.0_mgL/k_M, 1)^alpha_M_2021 * pmax(f.30700.0.0_mgL/k_M, 1)^-1.2 * 0.9938^f.21003.0.0, 2)
           ),
         logeGFR2009 =
           case_when(
             ckdepi <= 14.9999 ~ log(15),
             ckdepi > 15.0001 ~ log(ckdepi)
           )
  ) 

egfr_out <- egfr %>% 
  dplyr::select(
    FID,
    IID,
    f.34.0.0,
    f.22001.0.0,
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
    f.21003.0.0,
    f.30700.0.0,
    f.30700.0.0_mgL,
    egfr_cr_2021,
    logeGFR2009
  )

# visualize phenotype distributions
ggplot(egfr, aes(x = ckdepi)) +
  geom_histogram(alpha = 7/10, color = "black", fill = "white") +
  labs( x = "eGFR creatinine", y = "Frequency") +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12))
ggsave(here(project_dir, "egfr", "histogram_eGFR_creatinine_2009.png"), width = 15, height = 15, units = "cm")

ggplot(egfr, aes(x = logeGFR2009)) +
  geom_histogram(alpha = 7/10, color = "black", fill = "white") +
  labs( x = "ln(eGFR creatinine - winsorized)", y = "Frequency") +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12))
ggsave(here(project_dir, "egfr", "histogram_ln_eGFR_creatinine_2009_winsorized_at_15.png"), width = 15, height = 15, units = "cm")

# compare between the eGFR already in file and the new one
corr <- cor.test(egfr$ckdepi, egfr$egfr_cr_2021, method = "pearson")

egfr <- egfr %>%
  mutate(sex = 
           case_when(
             f.22001.0.0 == 0 ~ "female",
             f.22001.0.0 == 1 ~ "male"
           ))

ggplot(egfr, aes(x = ckdepi, y = egfr_cr_2021)) +
  geom_point(size = 2, alpha = 5/10, aes(colour = sex)) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  labs( x = "eGFR_creatinine_2009", y = "eGFR_creatinine_2021") +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12))
ggsave(here(project_dir, "egfr", "corr_eGFR_creatinine_2009_2021.png"), width = 15, height = 15, units = "cm")

# 2. compute eGFR(cystatine C) based on the 2012 equation:
# https://www.kidney.org/content/ckd-epi-cystatin-c-equation-2012

# 3. Compute eGFR(cystatine C) based on Moumita's recommendation: 
# 100/cystatin_C (i.e. field f.30720.0.0)
  # egfr <- egfr %>%
  #  mutate(egfr_cys_simple = 100/f.30720.0.0)

write.table(egfr_out, here(project_dir, "egfr", "egfr_withcovariates_for_REGENIE.txt"), sep = "\t", row.names = F, quote = F)

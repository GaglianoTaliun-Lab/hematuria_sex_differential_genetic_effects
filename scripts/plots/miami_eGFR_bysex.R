# Create Manhattan plots of the sex-specific eGFR GWASs
# SNP, CHR, POS, PVALUE (in that order, header name does not matter)

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(R.utils)

# Arguments ---------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
source(here(project_dir, "scripts", "ggmirror_FLD.R"))

# Main --------------------------------------

# Read summary statistics for males
males <- fread(here(project_dir, "results_regenie", "eGFR_male_HRC_imputed_regenie_allchrs.regenie")) %>%
    mutate(SNP = as.character(ID), CHR = CHROM, BP = as.numeric(GENPOS), P = 10^-LOG10P) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Read summary statistics for females
females <- fread(here(project_dir, "results_regenie", "eGFR_female_HRC_imputed_regenie_allchrs.regenie")) %>%
    mutate(SNP = as.character(ID), CHR = CHROM, BP = as.numeric(GENPOS), P = 10^-LOG10P) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Plot
p <- gmirror_v2(top=females, bottom=males, tline=5e-08, bline=5e-08,
toptitle="Females", bottomtitle="Males", 
highlight_p = c(5e-08,5e-08), highlighter="red", freey = FALSE)

# Save plot
ggsave(filename = here(project_dir, "sex_specific_eGFRcrea_miami_HRCimputed.png"),
       plot = p, units="in", height=7, width=12, dpi=300)
# Create qqplots of the sex-specific eGFR GWASs

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(R.utils)
library(qqman)

# Arguments ---------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

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
jpeg(filename = here(project_dir, "results_regenie", "qqplot_eGFRcrea_males.jpg"), width = 15, height = 15, units = "cm", res = 300)
qq(males$P)
dev.off()

jpeg(filename = here(project_dir, "results_regenie", "qqplot_eGFRcrea_females.jpg"), width = 15, height = 15, units = "cm", res = 300)
qq(females$P)
dev.off()
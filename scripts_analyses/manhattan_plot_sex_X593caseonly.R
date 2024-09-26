# Create Manhattan plots of the case only sex interaction for hematuria
# SNP, CHR, POS, PVALUE (in that order, header name does not matter)

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# Arguments ---------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
source(here(project_dir, "scripts", "manhattan_function.R"))

# Main --------------------------------------

# Read summary statistics

sumstats <- fread(here(project_dir, "results_regenie", "sex_X593caseonly_HRC_regenie_allchrs.regenie")) %>%
    mutate(SNP = as.character(ID), CHR = CHROM, BP = as.numeric(GENPOS), P = 10^-LOG10P) %>%
        select(
            SNP,
            CHR,
            BP,
            P)

ylim_max = 1+(-log10(min(sumstats$P)))
gg.manhattan(df = sumstats, threshold = 5e-8, colours = c("#996699", "#990099"), ylims = c(0,ylim_max), title = "Sex effect in hematuria cases")

# Save plot
ggsave(here(project_dir, "sex_effect_X593caseonly.png"), width = 50, height = 30, units = "cm")
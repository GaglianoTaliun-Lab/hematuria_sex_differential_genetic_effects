# Create Manhattan plots of the sex-specific burden tests
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

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
source(here(project_dir, "scripts", "ggmirror_FLD.R"))

# Main --------------------------------------

# Read summary statistics for males
males <- fread(here(project_dir, "burden", "regenie_burden_males_firth.txt")) %>%
  mutate(SNP = as.character(ID), BP = as.numeric(GENPOS), P = as.numeric(PVAL)) %>%
  select(
    SNP,
    CHR = CHROM,
    BP,
    P)

# Read summary statistics for females
females <- fread(here(project_dir, "burden", "regenie_burden_females_firth.txt")) %>%
  mutate(SNP = as.character(ID), BP = as.numeric(GENPOS), P = as.numeric(PVAL)) %>%
  select(
    SNP,
    CHR = CHROM,
    BP,
    P)

# Plot
p <- gmirror_v2(top=females, bottom=males, tline=3.29e-06, bline=3.29e-06,
                toptitle="Females", bottomtitle="Males", 
                highlight_p = c(3.29e-06,3.29e-06), highlighter="forestgreen", annotate_snp = "COL4A4", freey = FALSE)

# Save plot
ggsave(filename = here(project_dir, "burden", "sex_specific_burden_test_X593.png"),
       plot = p, units="in", height=7, width=12, dpi=300)

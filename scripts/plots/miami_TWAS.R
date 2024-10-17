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
males <- fread(here(project_dir, "spredixcan_output", "hematuria_males_gtex_v8_WholeBlood_mapped.tsv")) %>%
    mutate(SNP = as.character(gene_name), CHR = chromosome, BP = as.numeric(start_position), P = pvalue) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Read summary statistics for females
females <- fread(here(project_dir, "spredixcan_output", "hematuria_females_gtex_v8_WholeBlood_mapped.tsv")) %>%
    mutate(SNP = as.character(gene_name), CHR = chromosome, BP = as.numeric(start_position), P = pvalue) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Plot
p <- gmirror_v2(top=females, bottom=males, tline=0.05/nrow(females), bline=0.05/nrow(males),
toptitle="Females", bottomtitle="Males", 
highlight_p = c(0.05/nrow(females),0.05/nrow(males)), highlighter="red", freey = FALSE)

# Save plot
ggsave(filename = here(project_dir, "spredixcan_output", "sex_specific_X593_TWAS_GTEx_WholeBlood.png"),
       plot = p, units="in", height=7, width=12, dpi=300)
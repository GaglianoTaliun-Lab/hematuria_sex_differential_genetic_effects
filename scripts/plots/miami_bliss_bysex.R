# Create Manhattan plots of the sex-specific hematuria PWASs
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

args <- commandArgs(TRUE)
pwas_model <- as.character(args[1])

project_dir = "/home/fridald4/links/projects/def-gsarah/fridald4/col4a2_hematuria"
source(here(project_dir, "scripts", "ggmirror_FLD.R"))

# Main --------------------------------------

# Read summary statistics for males
males <- fread(here(project_dir, "results_bliss", stringr::str_c("PWAS_BLISS_", pwas_model, "_MALES.txt.finished"))) %>%
    mutate(SNP = as.character(protein), CHR = chromosome, BP = as.numeric(start), P = as.numeric(p)) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Read summary statistics for females
females <- fread(here(project_dir, "results_bliss", stringr::str_c("PWAS_BLISS_", pwas_model, "_FEMALES.txt.finished"))) %>%
mutate(SNP = as.character(protein), CHR = chromosome, BP = as.numeric(start), P = as.numeric(p)) %>%
select(
    SNP,
    CHR,
    BP,
    P)

pval_males=0.05/nrow(males)
pval_females=0.05/nrow(females)

# Plot
p <- gmirror_v2(top=females, bottom=males, tline=pval_females, bline=pval_males,
toptitle="Females", bottomtitle="Males", 
highlight_p = c(pval_females, pval_males), highlighter="red", freey = FALSE)

# Save plot
ggsave(filename = here(project_dir, "results_bliss", stringr::str_c("sex_specific_X593_PWAS_BLISS_", pwas_model, ".png")),
       plot = p, units="in", height=7, width=12, dpi=300)

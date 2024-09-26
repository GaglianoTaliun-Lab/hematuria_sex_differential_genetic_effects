# Create Manhattan plots of the sex-specific hematuria GWASs
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
males <- fread(here(project_dir, "datasets", "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")) %>%
    mutate(SNP = as.character(SNP), CHR, BP = as.numeric(BP), P = as.numeric(P)) %>%
    select(
        SNP,
        CHR,
        BP,
        P)

# Read summary statistics for females
females <- fread(here(project_dir, "datasets", "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")) %>%
mutate(SNP = as.character(SNP), CHR, BP = as.numeric(BP), P = as.numeric(P)) %>%
select(
    SNP,
    CHR,
    BP,
    P)

# Plot
p <- gmirror_v2(top=females, bottom=males, tline=5e-08, bline=5e-08,
toptitle="Females", bottomtitle="Males", 
highlight_p = c(5e-08,5e-08), highlighter="purple", freey = FALSE)

# Save plot
ggsave(filename = here(project_dir, "sex_specific_X593_miami_HRCimputed.png"),
       plot = p, units="in", height=7, width=12, dpi=300)
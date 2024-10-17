# Create Manhattan plots of the SDE summary results
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

sumstats <- fread(here(project_dir, "SDE", "SDE_Z_scores_genome_wide.tsv")) %>%
mutate(SNP = as.character(SNP), CHR, BP = as.numeric(BP), P = as.numeric(z_score_pval)) %>%
select(
    SNP,
    CHR,
    BP,
    P)

head(sumstats,5)
tail(sumstats,5)

ylim_max = 1+(-log10(min(sumstats$P)))
gg.manhattan(df = sumstats, threshold = 1e-6, colours = c("azure3", "azure4"), ylims = c(0,ylim_max), title = "SDE effects")

# Save plot
ggsave(here(project_dir, "SDE", "SDE_manhattan_plot.png"), width = 50, height = 30, units = "cm")
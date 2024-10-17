# Description: estimate sexually dimorphic effects in hematuria sex-specific GWAS

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read datasets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hem_m <- fread(here(project_dir, "datasets", "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")) %>%
    dplyr::select(
        CHR,
        BP,
        SNP,
        effect_allele,
        other_allele,
        BETA_M = BETA,
        SE_M = SE,
        P_M = P
    ) %>%
    dplyr::distinct(., SNP, .keep_all = TRUE)

hem_f <- fread(here(project_dir, "datasets", "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")) %>%
    dplyr::select(
        CHR,
        BP,
        SNP,
        effect_allele,
        other_allele,
        BETA_F = BETA,
        SE_F = SE,
        P_F = P
    ) %>%
    dplyr::distinct(., SNP, .keep_all = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate Z-score statistics
hematuria_sumstats <- inner_join(hem_m, hem_f, by = c("SNP", "CHR", "BP", "effect_allele", "other_allele"))

nrow(hematuria_sumstats)
head(hematuria_sumstats)

hematuria_sumstats <- hematuria_sumstats %>%
    mutate(z_score = (BETA_F - BETA_M)/(sqrt(SE_F^2 + SE_M^2)),
           z_score_pval = 2*pnorm(q=abs(z_score), lower.tail = FALSE))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(hematuria_sumstats, here(project_dir, "SDE", "SDE_Z_scores_genome_wide.tsv"), sep = "\t", row.names = F, quote = F)

hematuria_sumstats %>%
    filter(., z_score_pval <= 1e-05) %>%
    write.table(., here(project_dir, "SDE", "SDE_Z_scores_significant_1e-05.tsv"), sep = "\t", row.names = F, quote = F)

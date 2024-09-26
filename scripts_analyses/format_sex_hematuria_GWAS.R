# process GWAS summary statistics (sex-specific hematuria) after SAIGE GWAS

# Summary statistics must have the following columns:
# CHR, BP, SNP, A1, A2, A2_AF, BETA, SE, P, N

# summary statistics are in GRCh37 already
# effect allele:
#   hematuria sex-stratified = Allele2

# Packages -------------------------------------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
out_ext=".tsv"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria male-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 9,774 cases and 177,666 controls = 187,440

files_males <- list.files(here(project_dir, "datasets", "X593-male"), pattern = "X593.ukb*", full.names = TRUE)

hematuria_m_list <- lapply(files_males, function(x) fread(x) %>%
    filter(., AC_Allele2 > 20 & AC_Allele2 < 374860) %>% # MAC > 20
    # filter(., AF_Allele2 > 0.001 & AF_Allele2 < 0.999) %>%
    filter(., imputationInfo > 0.4) %>%
    dplyr::select(
      CHR,
      BP = POS,
      SNP = rsid,
      other_allele = Allele1,
      effect_allele = Allele2, # effect allele
      EAF = AF_Allele2,
      BETA,
      SE,
      P = p.value,
      N
    )
)

hematuria_m <- rbindlist(hematuria_m_list)

fwrite(hematuria_m, file = here(project_dir, "datasets", "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"), sep = "\t")

rm(hematuria_m_list)
rm(hematuria_m)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria female-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 7,061 cases and 213,033 controls = 220,094

files_females <- list.files(here(project_dir, "datasets", "X593-female"), pattern = "X593.ukb*", full.names = TRUE)

hematuria_f_list <- lapply(files_females, function(x) fread(x) %>%
    filter(., AC_Allele2 > 20 & AC_Allele2 < 440168) %>% # MAC > 20
    # filter(., AF_Allele2 > 0.001 & AF_Allele2 < 0.999) %>%
    filter(., imputationInfo > 0.4) %>%
    dplyr::select(
      CHR,
      BP = POS,
      SNP = rsid,
      other_allele = Allele1,
      effect_allele = Allele2, # effect allele
      EAF = AF_Allele2,
      BETA,
      SE,
      P = p.value,
      N
    )
)

hematuria_f <- rbindlist(hematuria_f_list)

fwrite(hematuria_f, file = here(project_dir, "datasets", "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"), sep = "\t")
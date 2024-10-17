# process GWAS summary statistics (eGFR sex-specific from UKB) colocalization of COL4A1/2 region on chr 13

# Summary statistics (quantitative traits) must have the following columns:
# snp (chr:bp_ref_alt), beta, varbeta, pvalues, MAF, N

# Region is limited to chr 13 bp: 110,801,310 - 250kb … 111,165,556 + 100kb (GRCh37)

# Packages -------------------------------------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

home_dir="/home/fridald4/projects/def-gsarah"
project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
out_ext=".txt"

# COL4A1/2 BP region to extract (GRCh37)
flanking_kb=100

min_bp = 110801310-(flanking_kb*1000)
max_bp = 111165556+(flanking_kb*1000)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | female eGFRcrea |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 209,516

file_path <- file.path(project_dir, "results_regenie", 
    "eGFR_female_HRC_imputed_regenie_allchrs.regenie")
fread(file_path, nrows = 6)

egfr_f <- fread(file_path) %>%
    filter(., CHROM == 13 & GENPOS >= min_bp & GENPOS <= max_bp) %>%
    dplyr::mutate(varbeta = SE^2,
                  P = 10^-LOG10P,
                  MAF = case_when(
                    A1FREQ >= 0.5 ~ 1-A1FREQ,
                    A1FREQ < 0.5 ~ A1FREQ
                  )) %>%
    dplyr::select(
      CHR = CHROM,
      BP = GENPOS,
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      SNP = ID,
      beta = BETA,
      varbeta,
      pvalues = P,
      MAF,
      N)

fwrite(egfr_f, file = here(project_dir, "data_for_coloc", stringr::str_c("eGFRcrea_females", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | male eGFRcrea |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 178,711

file_path <- file.path(project_dir, "results_regenie", 
    "eGFR_male_HRC_imputed_regenie_allchrs.regenie")
fread(file_path, nrows = 6)

egfr_m <- fread(file_path) %>%
    filter(., CHROM == 13 & GENPOS >= min_bp & GENPOS <= max_bp) %>%
    dplyr::mutate(varbeta = SE^2,
                  P = 10^-LOG10P,
                  MAF = case_when(
                    A1FREQ >= 0.5 ~ 1-A1FREQ,
                    A1FREQ < 0.5 ~ A1FREQ
                  )) %>%
    dplyr::select(
      CHR = CHROM,
      BP = GENPOS,
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      SNP = ID,
      beta = BETA,
      varbeta,
      pvalues = P,
      MAF,
      N)

fwrite(egfr_m, file = here(project_dir, "data_for_coloc", stringr::str_c("eGFRcrea_males", out_ext)), sep = "\t")


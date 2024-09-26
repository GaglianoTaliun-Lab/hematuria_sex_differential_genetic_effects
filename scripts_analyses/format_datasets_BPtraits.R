# process GWAS summary statistics (BP traits from Yang et al. 2024) colocalization of COL4A1/2 region on chr 13

# Summary statistics (quantitative traits) must have the following columns:
# snp (chr:bp_ref_alt), beta, varbeta, pvalues, MAF, N

# Region is limited to chr 13 bp: 110,801,310 - 250kb … 111,165,556 + 250kb (GRCh37)

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
#                                               | female DBP |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 174,664

file_path <- file.path(project_dir, "datasets", 
    "GCST90301699.h.DBP_females.tsv.gz")
fread(file_path, nrows = 6)

DBP_f <- fread(file_path) %>%
    separate(data = ., col = pos.hg19, into = c("CHR", "BP"), sep = ":") %>%
    filter(., CHR == 13 & BP >= min_bp & BP <= max_bp) %>%
    colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
    dplyr::mutate(varbeta = standard_error^2,
                  MAF = case_when(
                    effect_allele_frequency >= 0.5 ~ 1-effect_allele_frequency,
                    effect_allele_frequency < 0.5 ~ effect_allele_frequency
                  ),
                  N = 174664) %>%
    dplyr::select(
      CHR,
      BP,
      effect_allele,
      other_allele,
      SNP,
      beta,
      varbeta,
      pvalues = p_value,
      MAF,
      N)

DBP_f <- DBP_f %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

fwrite(DBP_f, file = here(project_dir, "data_for_coloc",stringr::str_c("DBP_females", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | female PP |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 174,664

file_path <- file.path(project_dir, "datasets", 
    "GCST90301702.h.PP_females.tsv.gz")
fread(file_path, nrows = 6)

PP_f <- fread(file_path) %>%
    separate(., pos.hg19, c("CHR", "BP"), sep = ":") %>%
    filter(., CHR == 13 & BP >= min_bp & BP <= max_bp) %>%
    colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
    dplyr::mutate(varbeta = standard_error^2,
                  MAF = case_when(
                    effect_allele_frequency >= 0.5 ~ 1-effect_allele_frequency,
                    effect_allele_frequency < 0.5 ~ effect_allele_frequency
                  ),
                  N = 174664) %>%
    dplyr::select(
      CHR,
      BP,
      effect_allele,
      other_allele,
      SNP,
      beta,
      varbeta,
      pvalues = p_value,
      MAF,
      N)

PP_f <- PP_f %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

fwrite(PP_f, file = here(project_dir, "data_for_coloc",stringr::str_c("PP_females", out_ext)), sep = "\t")


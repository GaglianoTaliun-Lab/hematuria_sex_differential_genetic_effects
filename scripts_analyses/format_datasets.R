# process GWAS summary statistics (hematuria and CAD - SC, F, M) colocalization of COL4A1/2 region on chr 13

# Summary statistics must have the following columns:
# snp (chr:bp_ref_alt), beta, varbeta, pvalues, MAF, N, s (ncase/N)

# Region is limited to chr 13 bp: 110,801,310 - 250kb … 111,165,556 + 250kb

# summary statistics are all in GRCh37 already
# effect alleles:
#   hematuria sex-combined = alt
#   hematuria sex-stratified = Allele2 (alt)
#   CAD = reference_allele -> check consistency

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
#                                               | hematuria sex-combined |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 16,409 cases and 379,936 controls

file_path <- file.path(project_dir, "datasets", 
    "hematuria_chr13_HRC_imputed_GRCh37_sexcombined.txt.gz")
fread(file_path, nrows = 6)

hematuria_sc <- fread(file_path) %>%
    filter(., pos >= min_bp & pos <= max_bp) %>%
    mutate(CHR = chrom, BP = pos) %>%
    colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
    dplyr::mutate(snp = stringr::str_c(CHR,":", BP,"_",ref,"_",alt),
                  varbeta = sebeta^2,
                  MAF = case_when(
                    af >= 0.5 ~ 1-af,
                    af < 0.5 ~ af
                  ),
                  N = 16409 + 379936,
                  s = 16409/(16409 + 379936)) %>%
    dplyr::select(
      CHR,
      BP,
      SNP,
      beta,
      varbeta,
      pvalues = pval,
      MAF,
      N,
      s)

hematuria_sc <- hematuria_sc %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

fwrite(hematuria_sc, file = here(project_dir, "data_for_coloc",stringr::str_c("hematuria_sexcombined", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria male-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 9,774 cases and 177,666 controls

file_path <- file.path(project_dir, "datasets", 
    "male-X593.ukb_chr13_v3.SAIGE.txt.gz")
fread(file_path, nrows = 6)

hematuria_m <- fread(file_path) %>%
    filter(., POS >= min_bp & POS <= max_bp) %>%
    filter(., AC_Allele2 > 20 & AC_Allele2 < 374860) %>%
    filter(., imputationInfo > 0.4) %>%
    dplyr::mutate(varbeta = SE^2,
                  MAF = case_when(
                    AF_Allele2 >= 0.5 ~ 1-AF_Allele2,
                    AF_Allele2 < 0.5 ~ AF_Allele2
                  ),
                  N = 9774 + 177666,
                  s = 9774/(9774 + 177666)) %>%
    dplyr::select(
      CHR,
      BP = POS,
      SNP = rsid,
      beta = BETA,
      varbeta,
      pvalues = p.value,
      MAF,
      N,
      s)

fwrite(hematuria_m, file = here(project_dir, "data_for_coloc",stringr::str_c("hematuria_males", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria female-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 7,061 cases and 213,033 controls

file_path <- file.path(project_dir, "datasets", 
    "female-X593.ukb_chr13_v3.SAIGE.txt.gz")
fread(file_path, nrows = 6)

hematuria_f <- fread(file_path) %>%
    filter(., POS >= min_bp & POS <= max_bp) %>%
    filter(., AC_Allele2 > 20 & AC_Allele2 < 440168) %>% # MAC > 20
    filter(., imputationInfo > 0.4) %>%
    dplyr::mutate(varbeta = SE^2,
                  MAF = case_when(
                    AF_Allele2 >= 0.5 ~ 1-AF_Allele2,
                    AF_Allele2 < 0.5 ~ AF_Allele2
                  ),
                  N = 7061 + 213033,
                  s = 7061/(7061 + 213033)) %>%
    dplyr::select(
      CHR,
      BP = POS,
      SNP = rsid,
      other_allele = Allele1,
      effect_allele = Allele2,
      beta = BETA,
      varbeta,
      pvalues = p.value,
      MAF,
      N,
      s)

fwrite(hematuria_f, file = here(project_dir, "data_for_coloc",stringr::str_c("hematuria_females", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CAD (read general file) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(home_dir, "Aragam_2022_CARDIoGRAM_CAD_GWAS", 
    "GWAMA.CAD.INTERMEDIATE.SENSITIVITY.UKBB.PLUS.EXTRA.STUDIES.QCed.txt.gz")
fread(file_path, nrows = 6)

cad <- fread(file_path)

# ~~~~~~~~~~~~
# check effect allele consistency among datasets (CAD and hematuria)
# ~~~~~~~~~~~~

hematuria_alleles <- fread(here(project_dir, "datasets", 
    "female-X593.ukb_chr13_v3.SAIGE.txt.gz")) %>%
  dplyr::select(rsid, 
                Allele1, 
                Allele2)

names(cad)

cad <- left_join(cad, hematuria_alleles, by = c('rsid_ukb'='rsid')) %>%
  mutate(effect_dir = case_when(
    Allele2 == reference_allele ~ 1,
    Allele1 == reference_allele ~ -1
  )) %>%
  dplyr::select(-Allele1, -Allele2) %>%
    mutate(beta = beta * effect_dir,
           male_beta = male_beta * effect_dir,
           female_beta = female_beta * effect_dir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CAD sex-combined |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size sex-combined: 120,788 cases and 856,183 controls

cad_sc <- cad %>%
    filter(., CHR == 13 & BP >= min_bp & BP <= max_bp) %>%
    dplyr::mutate(varbeta = se^2,
                  snp = stringr::str_c(CHR,":",BP,"_", reference_allele, "_", other_allele),
                  MAF = case_when(
                    eaf >= 0.5 ~ 1-eaf,
                    eaf < 0.5 ~ eaf
                  ),
                  N = 120788 + 856183,
                  s = 120788/(120788 + 856183)) %>%
    dplyr::select(
      CHR,
      BP,
      SNP = rsid_ukb,
      beta,
      varbeta,
      pvalues = p_value,
      MAF,
      N,
      s)

fwrite(cad_sc, file = here(project_dir, "data_for_coloc",stringr::str_c("cad_sexcombined", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CAD male-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size males: 82,060 cases and 440,505 controls

cad_m <- cad %>%
    filter(., CHR == 13 & BP >= min_bp & BP <= max_bp) %>%
    dplyr::mutate(varbeta = male_se^2,
                  snp = stringr::str_c(CHR,":",BP,"_", reference_allele, "_", other_allele),
                  MAF = case_when(
                    male_eaf >= 0.5 ~ 1-male_eaf,
                    male_eaf < 0.5 ~ male_eaf
                  ),
                  N = 82060 + 440505,
                  s = 82060/(82060 + 440505)) %>%
    dplyr::select(
      CHR,
      BP,
      SNP = rsid_ukb,
      beta = male_beta,
      varbeta,
      pvalues = male_p_value,
      MAF,
      N,
      s)

fwrite(cad_m, file = here(project_dir, "data_for_coloc",stringr::str_c("cad_males", out_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | CAD female-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size females: 38,728 cases and 415,678 controls

cad_f <- cad %>%
    filter(., CHR == 13 & BP >= min_bp & BP <= max_bp) %>%
    dplyr::mutate(varbeta = female_se^2,
                  snp = stringr::str_c(CHR,":",BP,"_", reference_allele, "_", other_allele),
                  MAF = case_when(
                    female_eaf >= 0.5 ~ 1-female_eaf,
                    female_eaf < 0.5 ~ female_eaf
                  ),
                  N = 38728 + 415678,
                  s = 38728/(38728 + 415678)) %>%
    dplyr::select(
      CHR,
      BP,
      SNP = rsid_ukb,
      beta = female_beta,
      varbeta,
      pvalues = female_p_value,
      MAF,
      N,
      s)

fwrite(cad_f, file = here(project_dir, "data_for_coloc",stringr::str_c("cad_females", out_ext)), sep = "\t")

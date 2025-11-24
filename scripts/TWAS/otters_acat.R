# Combine OTTERS TWAS p-values into one ACAT p-value

# Load packages -----------------------------------------------------------

library(dplyr)
library(here)
library(stringr)
library(tidyr)
library(data.table)
library(ACAT)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/hematuria_sexspecific"

args <- commandArgs(TRUE)
# 1) chromosome
# 2) sex
# 3) tissue
chr <- as.character(args[1])
sex <- as.character(args[2])
tissue <- as.character(args[3])

# TWAS output directory
twas_dir=here(project_dir, "otter_twas", stringr::str_c(tissue, "_", sex), stringr::str_c("chr", chr))

# Functions ---------------------------------------------------------------

get_pvalues = function(method){
  
  twas_file <- file.path(twas_dir, stringr::str_c(method, ".txt"))
  twas <- read.table(twas_file, header = T, fill = TRUE, na.strings = "nan")
  
  out <- twas %>%
    dplyr::select(CHROM, GeneStart, GeneEnd, TargetID, FUSION_PVAL)
  
  return(out)
  
}

ACAT_with_NA = function(p_vec){
  
  p_vec = p_vec[is.na(p_vec) == F]
  
  return(ACAT(p_vec))
  
}

# Main ---------------------------------------------------------------

lassosum <- get_pvalues("lassosum") %>% dplyr::rename(FUSION_PVAL_lassosum = FUSION_PVAL)
P0_001 <- get_pvalues("P0.001") %>% dplyr::rename(FUSION_PVAL_P0_001 = FUSION_PVAL)
P0_05 <- get_pvalues("P0.05") %>% dplyr::rename(FUSION_PVAL_P0_05 = FUSION_PVAL)
PRScs <- get_pvalues("PRScs") %>% dplyr::rename(FUSION_PVAL_PRScs = FUSION_PVAL)
SDPR <- get_pvalues("SDPR") %>% dplyr::rename(FUSION_PVAL_SDPR = FUSION_PVAL)

twas_all <- dplyr::full_join(lassosum, P0_001, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID")) %>%
    full_join(., P0_05, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID")) %>%
    full_join(., PRScs, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID")) %>%
    full_join(., SDPR, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID")) %>%
    rowwise() %>%
    mutate(otters_pval = ACAT_with_NA(c(FUSION_PVAL_P0_001, FUSION_PVAL_P0_05, FUSION_PVAL_PRScs, FUSION_PVAL_SDPR)))

write.table(twas_all, here(project_dir, "otter_acat", stringr::str_c(tissue, "_", sex), stringr::str_c("chr", chr, ".txt")), sep = "\t", row.names = F, col.names = F, quote = F)
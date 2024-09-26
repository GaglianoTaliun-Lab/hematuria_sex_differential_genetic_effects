# Description: subset of the script 'get_sample_overlaps.R' from Regina,
# to process after running 'ldsc.sh' across all pairs of traits and obtain a correlation matrix
# and the sample overlaps for hyprcoloc.

# Load packages -----------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readr)

# Set arguments -----------------------------------------------------------

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

# phenotypes to include in the matrix
phenotypes <- read.table(here(project_dir, "data_for_ldsc", "gwas_sample_sizes.tsv"), sep = "\t", header = F) %>% .[,1] %>% as.array() %>% sort()

# read and write input/output directories:
gwas_dir <- here(project_dir, "data_for_ldsc")

log_dir <- here(project_dir,"results_ldsc")
out_dir <- here(project_dir,"results_ldsc")

# number of lines to skip in the ldsc output:
n_skip=60
# number of lines to read from the output:
n_read=1

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    log_dir,
    pattern = "_rg.log",
    full.names = T
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){
  
  files[[i]] <-
    read_table(
      file = file_paths[i],
      skip = n_skip, # Number of lines to skip in log file
      n_max = n_read # Number of lines to read
    )
  
}

# extract rg values in matrix
all_rg <-
  files %>%
  rbindlist(., fill = TRUE) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/")),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz") %>%
      stringr::str_remove(str_c(gwas_dir,"/"))
  ) %>%
  select(
    p1,
    p2,
    rg,
    se,
    z,
    p,
    h2_obs,
    h2_obs_se,
    h2_int,
    h2_int_se,
    gcov_int,
    gcov_int_se
  )

###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

n <- length(phenotypes)
covar_matrix <- matrix(NA,n,n)

rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

for (k in 1:length(phenotypes)){
  
  for(i in phenotypes) {
    for(j in phenotypes) {

      gcov_int <- try(dplyr::filter(all_rg, p1 == i, p2 == j) %>% .[["gcov_int"]])

        if (inherits(gcov_int, "try-error")) {
            next
        }

      covar_matrix[i,j] <-
        gcov_int    

    }
  }
  
  # Sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
  if (!all(t(covar_matrix)==covar_matrix)) {
    covar_matrix[lower.tri(covar_matrix)] <- t(covar_matrix)[lower.tri(covar_matrix)]
  }
  
  # Standardise the matrix
  covar_matrix <-
    covar_matrix %>%
    cov2cor() %>%
    round(digits = 5)
  
}

# Save data ---------------------------------------------------------------

write.table(
  all_rg,
  file = here(out_dir, stringr::str_c("ldsc_correlations.txt")),
  sep = "\t",
  row.names = F,
  quote = F
)

write.table(
  covar_matrix,
  file = here(out_dir, stringr::str_c("sample_overlap.txt")),
  quote = F,
  row.names = phenotypes,
  col.names = NA,
  sep = "\t"
)

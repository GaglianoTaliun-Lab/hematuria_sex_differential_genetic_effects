library(RSQLite)
library(here)
library(tidyverse)
library(stringr)

project_dir="/Users/frida/Documents/research-projects/col4a2_hematuria/build_sex_stratified_PrediXcan_weights"

for (sex in c("males", "females", "all")) {
  
  filename <- here(project_dir, stringr::str_c("gtex_v8_imputed_b37_eur_Whole_Blood_", sex, "_v2.db"))
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver,
                  dbname = filename)
  
  dbListTables(db)
  modelsummary_info <- dbReadTable(db,"model_summaries") %>%
    dplyr::rename(genename = gene_name, 
                  n.snps.in.model = n_snps_in_model, 
                  pred.perf.R2 = rho_avg_squared,
                  pred.perf.pval = zscore_pval) %>%
    mutate(pred.perf.qval = NA)
  weights2 <- dbReadTable(db, "weights") %>%
    dplyr::rename(weight = beta, ref_allele = ref, eff_allele = alt)
  dbRemoveTable(db, "weights2")
  dbRemoveTable(db, "model_summaries")
  dbWriteTable(conn = db, name = "weights", weights2, header = TRUE)
  dbWriteTable(conn = db, name = "extra", modelsummary_info, header = TRUE)
  dbDisconnect(db)
}


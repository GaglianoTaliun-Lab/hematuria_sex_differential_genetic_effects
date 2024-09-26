# Map gene names from metaXcan output to chr and position.

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ---------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

# Main --------------------------------------

### read reference from Homo sapiens obtained from: https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
### this file includes gene symbols and IDs and chromosome, but not start/end position of the gene
gene_symbol_reference <- fread(here(project_dir, "reference_data", "Homo_sapiens.gene_info")) %>%
  filter(., chromosome %in% (1:22)) %>%
  filter(., type_of_gene == "protein-coding")

### read reference from Homo sapiens obtained from: https://ftp.ncbi.nih.gov/gene/DATA/gene_neighbors.gz
### this file includes gene IDs, chromosome, start/end positions, but not gene symbols.
gene_ID_reference <- fread(here(project_dir, "reference_data", "gene_neighbors")) %>%
  filter(., chromosome %in% (1:22))

### read output from metaXcan:
males <- read.table(here(project_dir, "spredixcan_output", "hematuria_males_gtex_v8_WholeBlood.csv"), sep = ",", header = T)
females <- read.table(here(project_dir, "spredixcan_output", "hematuria_females_gtex_v8_WholeBlood.csv"), sep = ",", header = T)

### join gene_symbol table with metaXcan output to include chr and geneID:
genes_mapped_males <- left_join(males, gene_symbol_reference, by = c("gene_name" = "Symbol"))
genes_mapped_females <- left_join(females, gene_symbol_reference, by = c("gene_name" = "Symbol"))

### join gene_ID table with genes_mapped list to include start/end positions of gene:
genes_mapped2_males <- left_join(genes_mapped_males, gene_ID_reference, by = c("GeneID", "chromosome"))
genes_mapped2_females <- left_join(genes_mapped_females, gene_ID_reference, by = c("GeneID", "chromosome"))

genes_mapped2_males %>%
  filter(., !is.na(start_position)) %>%
  dplyr::select(
    gene,
    gene_name,
    zscore,
    effect_size,
    pvalue,
    var_g,
    pred_perf_r2,
    pred_perf_pval,
    pred_perf_qval,
    n_snps_used,
    n_snps_in_cov,
    n_snps_in_model,
    GeneID,
    chromosome,
    map_location,
    start_position,
    end_position) %>%
    distinct(., gene_name, .keep_all = TRUE) %>%
  write.table(., here(project_dir, "spredixcan_output","hematuria_males_gtex_v8_WholeBlood_mapped.tsv"), sep = "\t", row.names = F, quote = F)

  genes_mapped2_females %>%
  filter(., !is.na(start_position)) %>%
  dplyr::select(
    gene,
    gene_name,
    zscore,
    effect_size,
    pvalue,
    var_g,
    pred_perf_r2,
    pred_perf_pval,
    pred_perf_qval,
    n_snps_used,
    n_snps_in_cov,
    n_snps_in_model,
    GeneID,
    chromosome,
    map_location,
    start_position,
    end_position) %>%
    distinct(., gene_name, .keep_all = TRUE) %>%
  write.table(., here(project_dir, "spredixcan_output","hematuria_females_gtex_v8_WholeBlood_mapped.tsv"), sep = "\t", row.names = F, quote = F)
# format GWAS summary statistics (sex-specific hematuria) for GCTA-COJO

# Summary statistics must have the following columns (in that order, column name doesn't matter):
# SNP EA OA Freq_EA BETA SE pvalue N

# summary statistics are all in GRCh37 already (bfile HRC imputed is GRCh37 too)
# effect alleles:
#   hematuria sex-stratified = Allele2 (alt)

# Packages --------------------------------------------------------------------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir="/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"
out_ext=".txt"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria female-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 7,061 cases and 213,033 controls

file_path <- file.path(project_dir, "datasets", 
    "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")

hematuria_f <- fread(file_path)

# SNP EA OA Freq_EA BETA SE pvalue N
hematuria_f %>%
    filter(., !is.na(SNP)) %>%
    filter(., SNP != "") %>%
    dplyr::select(
        SNP,
        effect_allele,
        other_allele,
        EAF,
        BETA,
        SE,
        P,
        N
    ) %>% write.table(., here(project_dir, "data_for_gcta", "hematuria_females.txt"), sep = "\t", row.names = F, quote = F)

hematuria_f_list <- hematuria_f %>%
    filter(., P < 5e-08) %>% 
    filter(., CHR != 6) %>% # remove HLA region
    group_split(CHR)

# obtain top SNP of each chromosome
hematuria_top_df <- lapply(hematuria_f_list, function(x) arrange(x, P) %>% .[1,]) %>% rbindlist()

top_females <- hematuria_top_df %>%
    mutate(sex = "females") %>%
    rownames_to_column(., "locus") %>%
    dplyr::select(locus, CHR, BP, SNP, sex)
write.table(top_females, here(project_dir, "data_for_gcta", "top_signals_hematuria_females.txt"), sep = "\t", row.names = F, quote = F, col.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | hematuria male-specific |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample size: 9,774 cases and 177,666 controls

file_path <- file.path(project_dir, "datasets", 
    "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz")

hematuria_m <- fread(file_path) 

# SNP EA OA Freq_EA BETA SE pvalue N
hematuria_m %>%
    filter(., !is.na(SNP)) %>%
    filter(., SNP != "") %>%
    dplyr::select(
        SNP,
        effect_allele,
        other_allele,
        EAF,
        BETA,
        SE,
        P,
        N
    ) %>% write.table(., here(project_dir, "data_for_gcta", "hematuria_males.txt"), sep = "\t", row.names = F, quote = F)

hematuria_m_list <- hematuria_m %>%
    filter(., P < 5e-08) %>% 
    filter(., CHR != 6) %>% # remove HLA region)
    group_split(CHR)

hematuria_top_df <- lapply(hematuria_m_list, function(x) arrange(x, P) %>% .[1,]) %>%
    rbindlist()

top_males <- hematuria_top_df %>%
    mutate(sex = "males") %>%
    rownames_to_column(., "locus") %>%
    dplyr::select(locus, CHR, BP, SNP, sex)
write.table(top_males, here(project_dir, "data_for_gcta", "top_signals_hematuria_males.txt"), sep = "\t", row.names = F, quote = F, col.names = F)

rbind(top_females, top_males) %>%
    write.table(., here(project_dir, "data_for_gcta", "top_signals_both_sexes.txt"), sep = "\t", row.names = F, quote = F, col.names = F)


##### testing for many independent signals in a chromosome:
#for (i in 1:nrow(hematuria_top_df)) {
#    chrom=hematuria_top_df$CHR[i]
#    BP_from=hematuria_top_df$BP[i]-500000
#    BP_to=hematuria_top_df$BP[i]+500000

    
#    hematuria_tmp <- filter(hematuria_f, CHR == chrom & BP <= BP_from & BP >= BP_to) %>%
#        arrange(., BP) %>% filter(., P < 5e-08)
#}

# list of all top signals per 500kb
#hematuria_signals <- lapply(hematuria_f_list, function(x) 
#                            arrange(x, BP, P) %>%
#                            mutate(diff = BP - lag(BP)) %>%
#                            filter(., diff > 500000) %>%
#                            dplyr::select(-diff)) %>%
#                            rbindlist()
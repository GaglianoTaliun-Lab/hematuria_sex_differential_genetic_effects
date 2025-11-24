library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(forcats)
library(vctrs)
library(locuszoomr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

# https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

# Read files -----------------------------------------------------

# p-values from female-only hematuria
hematuria_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_females.txt"), sep = "\t", header = T)

# p-values from male-only hematuria
hematuria_m <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_males.txt"), sep = "\t", header = T)

# p-values from female-only eGFRcrea
egfr_f <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_female_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
  mutate(., SNP = ID, P = 10^-LOG10P)

# p-values from male-only eGFRcrea
egfr_m <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_male_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
  mutate(., SNP = ID, P = 10^-LOG10P)

# LD from plink
ld <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix_chr13_COL4A2.phased.vcor1"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix_chr13_COL4A2.phased.vcor1.vars"), sep = "\t", header = F) 
colnames(ld) <- ld_snps$V1

###########################################################################################
# Plots with locuszoomr for hematuria -----------------------------------------------------
###########################################################################################

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

##### ------------- Hematuria females:

hematuria_f_ld <- left_join(hematuria_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs7323228,
         ld = case_when(
           LD_rs7323228 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs7323228 != 0 ~ LD_rs7323228
         )) %>%
  dplyr::filter(., !is.na(LD_rs7323228))

loc_f <- locus(data = hematuria_f_ld, gene = "COL4A2", flank = 1e2, LD = "r2",
             ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

# loc_f <- link_LD(loc_f, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_f <- link_recomb(loc_f)

summary(loc_f)
females_plot <- locus_ggplot(loc_f, labels = "index", italics = TRUE,
                             filter_gene_name = c("COL4A1", "COL4A2", "COL4A2-AS2"))

##### -------------- Hematuria males:

hematuria_m_ld <- left_join(hematuria_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs7323228,
         ld = case_when(
           LD_rs7323228 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs7323228 != 0 ~ LD_rs7323228
         )) %>%
  dplyr::filter(., !is.na(LD_rs7323228))

loc_m <- locus(data = hematuria_m_ld, gene = "COL4A2", flank = 1e2, LD = "r2", index_snp = "rs7323228",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

# loc_m <- link_LD(loc_m, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_m <- link_recomb(loc_m)

# summary(loc_m)
males_plot <- locus_ggplot(loc_m, labels = "index", 
                           ylim = c(0,max(-log10(hematuria_f$pvalues))),
                           nudge_y = 1.5, italics = TRUE, filter_gene_name = c("COL4A1", "COL4A2", "COL4A2-AS2"))

# Save plot ----------------------------------------------------------------

figure <- grid.arrange(females_plot, males_plot,
                       nrow = 1, widths = c(3.2, 3.2))

ggsave(here(project_dir, "regional_plots", "COL4A2_hematuria_sex_specific_locuszoomr_UKB_LD.jpg"),
       plot = figure,
       width = 14, height = 6, dpi = 300)

###########################################################################################
# Plots with locuszoomr for eGFRcrea -----------------------------------------------------
###########################################################################################

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs11838637 = ld$rs11838637) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs11838637 = round(LD_rs11838637^2,2))

##### ------------- eGFR females:

egfr_f_ld <- left_join(egfr_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs11838637,
         ld = case_when(
           LD_rs11838637 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs11838637 != 0 ~ LD_rs11838637
         )) %>%
  dplyr::filter(., !is.na(LD_rs11838637))

loc_f <- locus(data = egfr_f_ld, gene = "COL4A2", flank = 1e2, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHROM", pos = "GENPOS", labs = "SNP", p = "P")

# loc_f <- link_LD(loc_f, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_f <- link_recomb(loc_f)

summary(loc_f)
females_plot <- locus_ggplot(loc_f, labels = "index", italics = TRUE,
                             filter_gene_name = c("COL4A1", "COL4A2", "COL4A2-AS2"))

##### -------------- eGFR males:

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs9521719 = ld$rs9521719) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs9521719 = round(LD_rs9521719^2,2))

egfr_m_ld <- left_join(egfr_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs9521719,
         ld = case_when(
           LD_rs9521719 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs9521719 != 0 ~ LD_rs9521719
         )) %>%
  dplyr::filter(., !is.na(LD_rs9521719))

loc_m <- locus(data = egfr_m_ld, gene = "COL4A2", flank = 1e2, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHROM", pos = "GENPOS", labs = "SNP", p = "P")

# loc_m$index_snp <- "rs530725799"

# loc_m <- link_LD(loc_m, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_m <- link_recomb(loc_m)

# summary(loc_m)
males_plot <- locus_ggplot(loc_m, labels = "index", 
                           ylim = c(0,max(-log10(egfr_f_ld$P))), 
                           italics = TRUE, filter_gene_name = c("COL4A1", "COL4A2", "COL4A2-AS2"))

# Save plot ----------------------------------------------------------------

figure <- grid.arrange(females_plot, males_plot,
                       nrow = 1, widths = c(3.2, 3.2))

ggsave(here(project_dir, "regional_plots", "COL4A2_eGFRcrea_sex_specific_locuszoomr_UKB_LD.jpg"),
       plot = figure,
       width = 14, height = 6, dpi = 300)


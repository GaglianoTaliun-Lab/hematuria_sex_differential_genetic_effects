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

# ~~~~~~ FEMALE TRAITS

# 1) hematuria (lead SNP = rs7323228)
hematuria_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_females.txt"), sep = "\t", header = T)

# 2) eGFRcrea (lead SNP = rs11838637)
egfr_f <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_female_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
  mutate(., SNP = ID, P = 10^-LOG10P)

# 3) DBP (lead SNP = rs55940034)
dbp_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_DBP_females.txt"), sep = "\t", header = T)

# 4) PP (lead SNP = rs11838776)
pp_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_PP_females.txt"), sep = "\t", header = T)

# 5) CAD
cad_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_cad_females.txt"), sep = "\t", header = T)

# ~~~~~~ MALE TRAITS
# 1) hematuria (lead SNP = rs9515162)
# hematuria_m <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_males.txt"), sep = "\t", header = T)

# 2) eGFRcrea (lead SNP = rs530725799)
# egfr_m <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_male_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
#   mutate(., SNP = ID, P = 10^-LOG10P)

# LD from plink
ld <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1.vars"), sep = "\t", header = F) 
colnames(ld) <- ld_snps$V1

###########################################################################################
# Plots with locuszoomr for females -----------------------------------------------------
###########################################################################################

# 1) hematuria (X593)
lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

hematuria_f_ld <- left_join(hematuria_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs7323228,
         ld = case_when(
           LD_rs7323228 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs7323228 != 0 ~ LD_rs7323228
         )) %>%
  dplyr::filter(., !is.na(LD_rs7323228))

loc_x593 <- locus(data = hematuria_f_ld, gene = "COL4A2", flank = 1e5, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_x593 <- link_recomb(loc_x593)

summary(loc_x593)
X593_plot <- gg_scatter(loc_x593, labels = "index", xticks = FALSE, legend_pos = "topright",
                        ylim = c(0,max(-log10(dbp_f$pvalues))))

# 2) eGFRcrea
lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs11838637 = ld$rs11838637) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs11838637 = round(LD_rs11838637^2,2))

egfr_f_ld <- left_join(egfr_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs11838637,
         ld = case_when(
           LD_rs11838637 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs11838637 != 0 ~ LD_rs11838637
         )) %>%
  dplyr::filter(., !is.na(LD_rs11838637))

loc_egfr <- locus(data = egfr_f_ld, gene = "COL4A2", flank = 1e5, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHROM", pos = "GENPOS", labs = "SNP", p = "P")

loc_egfr <- link_recomb(loc_egfr)

summary(loc_egfr)
eGFR_plot <- gg_scatter(loc_egfr, labels = "index", xticks = FALSE, legend_pos = "topright",
                        ylim = c(0,max(-log10(dbp_f$pvalues))))

# 3) DBP
lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

dbp_f_ld <- left_join(dbp_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs55940034,
         ld = case_when(
           LD_rs55940034 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs55940034 != 0 ~ LD_rs55940034
         )) %>%
  dplyr::filter(., !is.na(LD_rs55940034))

loc_dbp <- locus(data = dbp_f_ld, gene = "COL4A2", flank = 1e5, LD = "r2",
                  ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_dbp <- link_recomb(loc_dbp)

summary(loc_dbp)
DBP_plot <- gg_scatter(loc_dbp, labels = "index", xticks = FALSE, legend_pos = "topright",
                       ylim = c(0,max(-log10(dbp_f$pvalues))))

# 4) PP
lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs11838776 = ld$rs11838776) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs11838776 = round(LD_rs11838776^2,2))

pp_f_ld <- left_join(pp_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs11838776,
         ld = case_when(
           LD_rs11838776 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs11838776 != 0 ~ LD_rs11838776
         )) %>%
  dplyr::filter(., !is.na(LD_rs11838776))

loc_pp <- locus(data = pp_f_ld, gene = "COL4A2", flank = 1e5, LD = "r2",
                 ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_pp <- link_recomb(loc_pp)

summary(loc_pp)
PP_plot <- gg_scatter(loc_pp, labels = "index", xticks = FALSE, legend_pos = "topright",
                      ylim = c(0,max(-log10(dbp_f$pvalues))))

# 5) CAD
lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs11838776 = ld$rs11838776) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs11838776 = round(LD_rs11838776^2,2))

cad_f_ld <- left_join(cad_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs11838776,
         ld = case_when(
           LD_rs11838776 == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD_rs11838776 != 0 ~ LD_rs11838776
         )) %>%
  dplyr::filter(., !is.na(LD_rs11838776))

loc_cad <- locus(data = cad_f_ld, gene = "COL4A2", flank = 1e5, LD = "r2",
                ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_cad <- link_recomb(loc_cad)

summary(loc_cad)
CAD_plot <- locus_ggplot(loc_cad, labels = "index", italics = TRUE,
                        filter_gene_name = c("COL4A1", "COL4A2", "COL4A2-AS2"), legend_pos = "topright",
                        ylim = c(0,max(-log10(dbp_f$pvalues))))

# Save plot ----------------------------------------------------------------

figure <- grid.arrange(X593_plot, eGFR_plot, DBP_plot, PP_plot, CAD_plot, 
                       ncol = 1, heights = c(5, 5, 5, 5, 7))

ggsave(here(project_dir, "regional_plots", "coloc_female_specific_locuszoomr_UKB_LD.jpg"),
       plot = figure,
       width = 7, height = 18, dpi = 300)


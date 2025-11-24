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

# ~~~~~~ MALE TRAITS
# 1) hematuria (lead SNP = rs9515162)
hematuria_m <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_males.txt"), sep = "\t", header = T)

# 2) eGFRcrea (lead SNP = rs9521719)
egfr_m <- read.table(here(project_dir, "regional_plots_data", "eGFRcrea_males.txt"), sep = "\t", header = T)

# 3) DBP (lead SNP = rs693119)
dbp_m <- read.table(here(project_dir, "regional_plots_data", "DBP_males.txt"), sep = "\t", header = T)

# 4) PP (lead SNP = rs150847806)
pp_m <- read.table(here(project_dir, "regional_plots_data", "PP_males.txt"), sep = "\t", header = T)

# 5) CAD (lead SNP = rs9515203)
cad_m <- read.table(here(project_dir, "regional_plots_data", "cad_males.txt"), sep = "\t", header = T)

# LD from plink
ld <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1.vars"), sep = "\t", header = F) 
colnames(ld) <- ld_snps$V1

###########################################################################################
# Plots with locuszoomr for males -----------------------------------------------------
###########################################################################################

# 1) hematuria (X593)
lead_snp <- data.frame(SNP = ld_snps$V1, LD = ld$rs9515162) %>%  # filter to retain only LD for index SNP
  mutate(LD = round(LD^2,2))

hematuria_m_ld <- left_join(hematuria_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD,
         ld = case_when(
           LD == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD != 0 ~ LD
         )) %>%
  dplyr::filter(., !is.na(LD))

loc_x593 <- locus(data = hematuria_m_ld, seqname = 13, xrange = c(110701370,111265312), LD = "r2",
                  ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_x593 <- link_recomb(loc_x593)

summary(loc_x593)
X593_plot <- gg_scatter(loc_x593, labels = "index", xticks = FALSE, legend_pos = "topright",
                        ylim = c(0,max(-log10(cad_m$pvalues))))

# 2) eGFRcrea
lead_snp <- data.frame(SNP = ld_snps$V1, LD = ld$rs9521719) %>%  # filter to retain only LD for index SNP
  mutate(LD = round(LD^2,2))

egfr_m_ld <- left_join(egfr_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD,
         ld = case_when(
           LD == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD != 0 ~ LD
         )) %>%
  dplyr::filter(., !is.na(LD))

loc_egfr <- locus(data = egfr_m_ld, seqname = 13, xrange = c(110701370,111265312), LD = "r2",
                  ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_egfr <- link_recomb(loc_egfr)

summary(loc_egfr)
egfr_plot <- gg_scatter(loc_egfr, labels = "index", xticks = FALSE, legend_pos = "topright",
                        ylim = c(0,max(-log10(cad_m$pvalues))))

# 3) DBP
lead_snp <- data.frame(SNP = ld_snps$V1, LD = ld$rs693119) %>%  # filter to retain only LD for index SNP
  mutate(LD = round(LD^2,2))

dbp_m_ld <- left_join(dbp_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD,
         ld = case_when(
           LD == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD != 0 ~ LD
         )) %>%
  dplyr::filter(., !is.na(LD))

loc_dbp <- locus(data = dbp_m_ld, seqname = 13, xrange = c(110701370,111265312), LD = "r2",
                  ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_dbp <- link_recomb(loc_dbp)

summary(loc_dbp)
dbp_plot <- gg_scatter(loc_dbp, labels = "index", xticks = FALSE, legend_pos = "topright",
                       ylim = c(0,max(-log10(cad_m$pvalues))))

# 4) PP
lead_snp <- data.frame(SNP = ld_snps$V1, LD = ld$rs150847806) %>%  # filter to retain only LD for index SNP
  mutate(LD = round(LD^2,2))

pp_m_ld <- left_join(pp_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD,
         ld = case_when(
           LD == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD != 0 ~ LD
         )) %>%
  dplyr::filter(., !is.na(LD))

loc_pp <- locus(data = pp_m_ld, seqname = 13, xrange = c(110701370,111265312), LD = "r2",
                 ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_pp <- link_recomb(loc_pp)

summary(loc_pp)
pp_plot <- gg_scatter(loc_pp, labels = "index", xticks = FALSE, legend_pos = "topright",
                      ylim = c(0,max(-log10(cad_m$pvalues))))

# 5) CAD
lead_snp <- data.frame(SNP = ld_snps$V1, LD = ld$rs9515203) %>%  # filter to retain only LD for index SNP
  mutate(LD = round(LD^2,2))

cad_m_ld <- left_join(cad_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD,
         ld = case_when(
           LD == 0 ~ 0.0000000001, # if it is zero then it appears as NA in the plot
           LD != 0 ~ LD
         )) %>%
  dplyr::filter(., !is.na(LD))

loc_cad <- locus(data = cad_m_ld, seqname = 13, xrange = c(110701370,111265312), LD = "r2",
                ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "BP", labs = "SNP", p = "pvalues")

loc_cad <- link_recomb(loc_cad)

summary(loc_cad)
cad_plot <- locus_ggplot(loc_cad, labels = "index", legend_pos = "topright", italics = TRUE,
                         ylim = c(0,max(-log10(cad_m$pvalues))))

# Save plot ----------------------------------------------------------------

figure <- grid.arrange(X593_plot, egfr_plot, dbp_plot, pp_plot, cad_plot, 
                       ncol = 1, heights = c(5, 5, 5, 5, 7))

ggsave(here(project_dir, "regional_plots", "coloc_male_specific_locuszoomr_UKB_LD.jpg"),
       plot = figure,
       width = 7, height = 18, dpi = 300)


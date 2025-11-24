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
hematuria_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_chr19_females.txt"), sep = " ", header = T) %>%
  mutate(., SNP = rsid)

# p-values from male-only hematuria
hematuria_m <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_chr19_males.txt"), sep = " ", header = T) %>%
  mutate(., SNP = rsid)

# LD from plink
ld <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix_chr19.phased.vcor1"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix_chr19.phased.vcor1.vars"), sep = "\t", header = F) 
colnames(ld) <- ld_snps$V1

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs56254331 = ld$rs56254331) %>%
  mutate(LD_rs56254331 = round(LD_rs56254331^2,2))

# Plots with locuszoomr -----------------------------------------------------

##### ------------- Females:

hematuria_f_ld <- left_join(hematuria_f, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs56254331,
         ld = case_when(
           LD_rs56254331 == 0 ~ 0.0000000001,
           LD_rs56254331 != 0 ~ LD_rs56254331
         )) %>%
  dplyr::filter(., !is.na(LD_rs56254331))

loc_f <- locus(data = hematuria_f_ld, gene = "TGFB1", flank = 1e4, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "POS", labs = "SNP", p = "p.value")

# loc_f <- link_LD(loc_f, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_f <- link_recomb(loc_f)

summary(loc_f)
females_plot <- locus_ggplot(loc_f, labels = "index", italics = TRUE, ylim = c(0,max(-log10(hematuria_m$p.value))))

##### -------------- Males:

hematuria_m_ld <- left_join(hematuria_m, lead_snp, by = "SNP") %>%
  mutate(ld = LD_rs56254331,
         ld = case_when(
           LD_rs56254331 == 0 ~ 0.0000000001,
           LD_rs56254331 != 0 ~ LD_rs56254331
         )) %>%
  dplyr::filter(., !is.na(LD_rs56254331))

loc_m <- locus(data = hematuria_m_ld, gene = "TGFB1", flank = 1e4, LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75", chrom = "CHR", pos = "POS", labs = "SNP", p = "p.value")

loc_m$index_snp <- "rs56254331"

# loc_m <- link_LD(loc_m, token = "bf1f585049ed", pop = c("GBR", "CEU", "TSI", "IBS"), genome_build = 'grch37')
loc_m <- link_recomb(loc_m)

summary(loc_m)
males_plot <- locus_ggplot(loc_m, labels = "index", italics = TRUE)

# Save plot ----------------------------------------------------------------

figure <- grid.arrange(females_plot, males_plot,
                       nrow = 1, widths = c(3.2, 3.2))

ggsave(here(project_dir, "regional_plots", "TGFB1_hematuria_sex_specific_locuszoomr_UKB_LD.jpg"),
       plot = figure,
       width = 14, height = 6, dpi = 300)


# Description: compare GWAS on the COL4A2 region between SAIGE and REGENIE

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_regenie <- read.table(here(project_dir, "regenie", "results", "X593_TOPMED_regenie_chr13_COL4A2_with_rsids.regenie"), sep = "\t", header = T)
female_regenie <- read.table(here(project_dir, "regenie", "results", "X593_TOPMED_female_regenie_chr13_COL4A2_with_rsids.regenie"), sep = "\t", header = T)
male_regenie <- read.table(here(project_dir, "regenie", "results", "X593_TOPMED_male_regenie_chr13_COL4A2_with_rsids.regenie"), sep = "\t", header = T) 

all_saige <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_sexcombined.txt"), sep = "\t", header = T) %>%
  mutate(se = sqrt(varbeta))
female_saige <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_females.txt"), sep = "\t", header = T) %>%
  mutate(se = sqrt(varbeta))
male_saige <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_males.txt"), sep = "\t", header = T)  %>%
  mutate(se = sqrt(varbeta))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# merge per sex
all <- inner_join(all_regenie, all_saige, by = "SNP", suffix = c(".regenie", ".saige")) %>%
  mutate(BETA = BETA*-1)
female <- inner_join(female_regenie, female_saige, by = "SNP", suffix = c(".regenie", ".saige")) %>%
  mutate(BETA = BETA*-1)
male <- inner_join(male_regenie, male_saige, by = "SNP", suffix = c(".regenie", ".saige")) %>%
  mutate(BETA = BETA*-1)

data_list <- list(all, female, male)

top_list <- lapply(data_list, function(x) filter(x, P <= 5e-08 | pvalues <= 5e-08))

for(sex in c("all", "females", "males")) {
  
  if(sex == "all") {
    list_el = 1
  } else if (sex == "females") {
    list_el = 2
  } else if (sex == "males") {
    list_el = 3
  }

  # pvalues
  ggplot(data_list[[list_el]], aes(x = -log10(P), y = -log10(pvalues))) +
    geom_point(size = 2, alpha = 6/10) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    labs( x = "-log10(P) REGENIE", y = "-log10(P) SAIGE") +
    theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top_list[[list_el]], colour = "red", size = 2, alpha = 7/10) +
    geom_hline(yintercept=-log10(5e-08), colour = "grey", linetype="dashed") +
    geom_vline(xintercept=-log10(5e-08), colour = "grey", linetype="dashed")
  ggsave(here(project_dir, "comparison_HRC_TOPMed", stringr::str_c("pvalues_chr13_COL4A2_", sex, ".png")), width = 15, height = 15, units = "cm")
  
  # betas
  ggplot(data_list[[list_el]], aes(x = BETA, y = beta)) +
    geom_point(size = 2, alpha = 6/10) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    labs( x = "effect size REGENIE", y = "effect size SAIGE") +
    theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top_list[[list_el]], colour = "red", size = 2, alpha = 7/10) +
    geom_hline(yintercept=0, colour = "grey", linetype="dashed") +
    geom_vline(xintercept=0, colour = "grey", linetype="dashed")
  ggsave(here(project_dir, "comparison_HRC_TOPMed", stringr::str_c("betas_chr13_COL4A2_", sex, ".png")), width = 15, height = 15, units = "cm")

}

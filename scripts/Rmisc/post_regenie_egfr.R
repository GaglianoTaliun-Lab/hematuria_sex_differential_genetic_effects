# Description: script to wrangle regenie eGFR association results (chr 13 COL4A2 with TOPMed in UKB RAP)

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

males <- read.table(here(project_dir, "regenie", "results", "logeGFR2009_TOPMED_male_regenie_chr13_COL4A2_logeGFR2009.regenie"), sep = " ", header = T) %>%
  mutate(CHR = CHROM, BP = GENPOS) %>%
  colochelpR::convert_loc_to_rs(., dbSNP = dbsnp_144) %>%
  mutate(P = 10^-LOG10P,
         SNP = case_when(
           is.na(SNP) ~ stringr::str_c(CHR,":", BP),
           !is.na(SNP) ~ SNP
         ),
         CHR_BP = stringr::str_c("chr", CHR, "_", BP)) %>%
  dplyr::select(-CHROM, -GENPOS, -ID, -EXTRA)
write.table(males, here(project_dir, "regenie", "results", "logeGFR2009_TOPMED_male_regenie_chr13_COL4A2_with_rsids.regenie"), sep = "\t", row.names = F, quote = F)


females <- read.table(here(project_dir, "regenie", "results", "logeGFR2009_TOPMED_female_regenie_chr13_COL4A2_logeGFR2009.regenie"), sep = " ", header = T) %>%
  mutate(CHR = CHROM, BP = GENPOS) %>%
  colochelpR::convert_loc_to_rs(., dbSNP = dbsnp_144) %>%
  mutate(P = 10^-LOG10P,
         SNP = case_when(
           is.na(SNP) ~ stringr::str_c(CHR,":", BP),
           !is.na(SNP) ~ SNP
         ),
         CHR_BP = stringr::str_c("chr", CHR, "_", BP)) %>%
  dplyr::select(-CHROM, -GENPOS, -ID, -EXTRA)
write.table(females, here(project_dir, "regenie", "results", "logeGFR2009_TOPMED_female_regenie_chr13_COL4A2_with_rsids.regenie"), sep = "\t", row.names = F, quote = F)

chr13 <- list(females = females, males = males)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regional plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load LD:
# LD from plink with TOPMed imputed:
ld <- read.table(here(project_dir, "regional_plots_data", "COL4A2_LD_TOPMed_imputed.ld"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "COL4A2_LD_TOPMed_imputed.bim"), sep = "\t", header = F)
colnames(ld_snps) <- c("CHR", "rsid", "cm", "BP", "A1", "A2")
ld_snps <- ld_snps %>% 
  mutate(CHR_BP = stringr::str_c("chr", CHR, "_", BP)) %>%
  dplyr::select(CHR, BP, CHR_BP)


colnames(ld) <- ld_snps$CHR_BP

for (sex in c("females", "males")) {
  
  if(sex == "females") {
    list_el = 1
  } else if (sex == "males") {
    list_el = 2
  }
  
  # 1) choose rs7323228 as index SNP:
  lead_snp <- data.frame(CHR_BP = ld_snps$CHR_BP, LD_rs7323228 = ld$chr13_110412070) %>%  # filter to retain only LD for index SNP
    mutate(LD_rs7323228 = round(LD_rs7323228^2,2))
  
  p1 <- chr13[[list_el]] %>%
    left_join(., lead_snp, by = "CHR_BP") %>%
    filter(., !is.na(LD_rs7323228)) %>%
    mutate(log10p = -log10(P))
  
  rsq=expression(paste(r^2))
  
  top <- p1 %>%
    filter(., SNP == "rs7323228")
  
  p1 <- p1 %>%
    mutate(is_top = case_when(
      SNP == top$SNP[1] ~ "annotate",
      SNP != top$SNP[1] ~ "none"
    ))
  
  # plot the p-value results (Y axis; colored by LD).
  ggplot(p1, aes_string(x = "BP", colour = "LD_rs7323228", y = "log10p")) +
    scale_y_continuous() +  
    geom_point(colour = "black", size = 2, alpha = 7/10) +
    geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
    #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
    scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    guides(fill="none") + labs( x = "Chromosome 13 (GRCh38)", y = "-log10(p-value)") +
    theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
    geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
    geom_text_repel(data=p1[p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
  ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_logeGFR2009_", sex,"_TOPMed_imputed_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")
  
  # 2) choose rs11838637 as index SNP:
  lead_snp <- data.frame(CHR_BP = ld_snps$CHR_BP, LD_rs11838637 = ld$chr13_110377576) %>%  # filter to retain only LD for index SNP
    mutate(LD_rs11838637 = round(LD_rs11838637^2,2))
  
  p2 <- chr13[[list_el]] %>%
    left_join(., lead_snp, by = "CHR_BP") %>%
    filter(., !is.na(LD_rs11838637)) %>%
    mutate(log10p = -log10(P))
  
  rsq=expression(paste(r^2))
  
  top <- p2 %>%
    filter(., SNP == "rs11838637")
  
  p2 <- p2 %>%
    mutate(is_top = case_when(
      SNP == top$SNP[1] ~ "annotate",
      SNP != top$SNP[1] ~ "none"
    ))
  
  # plot the p-value results (Y axis; colored by LD).
  ggplot(p2, aes_string(x = "BP", colour = "LD_rs11838637", y = "log10p")) +
    scale_y_continuous() +  
    geom_point(colour = "black", size = 2, alpha = 7/10) +
    geom_point(aes(colour = LD_rs11838637), size = 1.5, shape = 19, alpha = 7/10) +
    #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
    scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    guides(fill="none") + labs( x = "Chromosome 13 (GRCh38)", y = "-log10(p-value)") +
    theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
    geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
    geom_text_repel(data=p2[p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
  ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_logeGFR2009_", sex, "_TOPMed_imputed_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")
  
  # 3) choose rs55940034 as index SNP:
  lead_snp <- data.frame(CHR_BP = ld_snps$CHR_BP, LD_rs55940034 = ld$chr13_110390962) %>%  # filter to retain only LD for index SNP
    mutate(LD_rs55940034 = round(LD_rs55940034^2,2))
  
  p2 <- chr13[[list_el]] %>%
    left_join(., lead_snp, by = "CHR_BP") %>%
    filter(., !is.na(LD_rs55940034)) %>%
    mutate(log10p = -log10(P))
  
  rsq=expression(paste(r^2))
  
  top <- p2 %>%
    filter(., SNP == "rs55940034")
  
  p2 <- p2 %>%
    mutate(is_top = case_when(
      SNP == top$SNP[1] ~ "annotate",
      SNP != top$SNP[1] ~ "none"
    ))
  
  # plot the p-value results (Y axis; colored by LD).
  ggplot(p2, aes_string(x = "BP", colour = "LD_rs55940034", y = "log10p")) +
    scale_y_continuous() +  
    geom_point(colour = "black", size = 2, alpha = 7/10) +
    geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
    #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
    scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    guides(fill="none") + labs( x = "Chromosome 13 (GRCh38)", y = "-log10(p-value)") +
    theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
    geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
    geom_text_repel(data=p2[p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
  ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_logeGFR2009_", sex, "_TOPMed_imputed_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")
  
  
  # 4) choose rs9521719 as index SNP:
  lead_snp <- data.frame(CHR_BP = ld_snps$CHR_BP, LD_rs9521719 = ld$chr13_110365437) %>%  # filter to retain only LD for index SNP
    mutate(LD_rs9521719 = round(LD_rs9521719^2,2))
  
  p2 <- chr13[[list_el]] %>%
    left_join(., lead_snp, by = "CHR_BP") %>%
    filter(., !is.na(LD_rs9521719)) %>%
    mutate(log10p = -log10(P))
  
  rsq=expression(paste(r^2))
  
  top <- p2 %>%
    filter(., SNP == "rs9521719")
  
  p2 <- p2 %>%
    mutate(is_top = case_when(
      SNP == top$SNP[1] ~ "annotate",
      SNP != top$SNP[1] ~ "none"
    ))
  
  # plot the p-value results (Y axis; colored by LD).
  ggplot(p2, aes_string(x = "BP", colour = "LD_rs9521719", y = "log10p")) +
    scale_y_continuous() +  
    geom_point(colour = "black", size = 2, alpha = 7/10) +
    geom_point(aes(colour = LD_rs9521719), size = 1.5, shape = 19, alpha = 7/10) +
    #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
    scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
    theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
    guides(fill="none") + labs( x = "Chromosome 13 (GRCh38)", y = "-log10(p-value)") +
    theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
    theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
    geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
    geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
    geom_text_repel(data=p2[p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
  ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_logeGFR2009_", sex, "_TOPMed_imputed_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

}


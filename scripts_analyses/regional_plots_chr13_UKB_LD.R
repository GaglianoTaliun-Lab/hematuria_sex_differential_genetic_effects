# Description: get regional plots with UKB HRC-imputed LD

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(LDlinkR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"
my_token = "bf1f585049ed"  # LDlink token

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# p-values from female-only hematuria
hematuria_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_hematuria_females.txt"), sep = "\t", header = T)

# LD from plink
ld <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1"), sep = "\t", header = F)
ld_snps <- read.table(here(project_dir, "regional_plots_data", "LD_r_matrix.phased.vcor1.vars"), sep = "\t", header = F) 
colnames(ld) <- ld_snps$V1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 1: regional plot of hematuria female-specific chr13 COL4A2 region
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

hematuria_f_p1 <- hematuria_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs7323228)) %>%
  mutate(log10p = -log10(pvalues))
  
rsq=expression(paste(r^2))

top <- hematuria_f_p1 %>%
  filter(., SNP == "rs7323228")

hematuria_f_p1 <- hematuria_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD r^2).
ggplot(hematuria_f_p1, aes_string(x = "BP", colour = "LD_rs7323228", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
    axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=hematuria_f_p1[hematuria_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_hematuria_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 2: regional plot of female-specific chr13 COL4A2 region, with rs55940034 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# no need to load hematuria sumstats again
# no need to load LD again

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

hematuria_f_p2 <- hematuria_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs55940034)) %>%
  mutate(log10p = -log10(pvalues))

rsq=expression(paste(r^2))

top <- hematuria_f_p2 %>%
  filter(., SNP == "rs55940034")

hematuria_f_p2 <- hematuria_f_p2 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(hematuria_f_p2, aes_string(x = "BP", colour = "LD_rs55940034", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=hematuria_f_p2[hematuria_f_p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_hematuria_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 3: regional plot of female-specific DBP chr13 COL4A2 region, with rs55940034 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region

# load DBP regional data
dbp_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_DBP_females.txt"), sep = "\t", header = T) %>%
  filter(., BP >= min(hematuria_f$BP) & BP <= max(hematuria_f$BP))

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

dbp_f_p1 <- dbp_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs55940034)) %>%
  mutate(log10p = -log10(pvalues))

rsq=expression(paste(r^2))

top <- dbp_f_p1 %>%
  filter(., SNP == "rs55940034")

dbp_f_p1 <- dbp_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(dbp_f_p1, aes_string(x = "BP", colour = "LD_rs55940034", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=dbp_f_p1[dbp_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_DBP_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 4: regional plot of female-specific DBP chr13 COL4A2 region, with rs7323228 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region
# no need to reload DBP data

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

dbp_f_p2 <- dbp_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs7323228)) %>%
  mutate(log10p = -log10(pvalues))

rsq=expression(paste(r^2))

top <- dbp_f_p2 %>%
  filter(., SNP == "rs7323228")

dbp_f_p2 <- dbp_f_p2 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(dbp_f_p2, aes_string(x = "BP", colour = "LD_rs7323228", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=dbp_f_p2[dbp_f_p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_DBP_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 5: regional plot of female-specific PP chr13 COL4A2 region, with rs55940034 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region

# load PP regional data
pp_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_PP_females.txt"), sep = "\t", header = T) %>%
  filter(., BP >= min(hematuria_f$BP) & BP <= max(hematuria_f$BP))

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

pp_f_p1 <- pp_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs55940034)) %>%
  mutate(log10p = -log10(pvalues))

rsq=expression(paste(r^2))

top <- pp_f_p1 %>%
  filter(., SNP == "rs55940034")

pp_f_p1 <- pp_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(pp_f_p1, aes_string(x = "BP", colour = "LD_rs55940034", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=pp_f_p1[pp_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_PP_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 6: regional plot of female-specific PP chr13 COL4A2 region, with rs7323228 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region
# no need to reload PP data

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

pp_f_p2 <- pp_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs7323228)) %>%
  mutate(log10p = -log10(pvalues))

rsq=expression(paste(r^2))

top <- pp_f_p2 %>%
  filter(., SNP == "rs7323228")

pp_f_p2 <- pp_f_p2 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(pp_f_p2, aes_string(x = "BP", colour = "LD_rs7323228", y = "log10p")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=pp_f_p2[pp_f_p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_PP_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 7: regional plot of female-specific eGFR chr13 COL4A2 region, with rs11838637 (lead SNP) as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region

# load eGFR regional data
egfr_f <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_female_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
  filter(., GENPOS >= min(hematuria_f$BP) & GENPOS <= max(hematuria_f$BP)) %>%
  mutate(., SNP = ID, P = 10^-LOG10P)

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs11838637 = ld$rs11838637) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs11838637 = round(LD_rs11838637^2,2))

egfr_f_p1 <- egfr_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs11838637))

rsq=expression(paste(r^2))

top <- egfr_f_p1 %>%
  filter(., SNP == "rs11838637")

egfr_f_p1 <- egfr_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(egfr_f_p1, aes_string(x = "GENPOS", colour = "LD_rs11838637", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs11838637), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=egfr_f_p1[egfr_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_eGFRcrea_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 8: regional plot of female-specific eGFR chr13 COL4A2 region, with rs7323228 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region
# no need to reload eGFR sumstats

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

egfr_f_p1 <- egfr_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs7323228))

rsq=expression(paste(r^2))

top <- egfr_f_p1 %>%
  filter(., SNP == "rs7323228")

egfr_f_p1 <- egfr_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(egfr_f_p1, aes_string(x = "GENPOS", colour = "LD_rs7323228", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=egfr_f_p1[egfr_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_eGFRcrea_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 8: regional plot of female-specific eGFR chr13 COL4A2 region, with rs55940034 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region
# no need to reload eGFR sumstats

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

egfr_f_p1 <- egfr_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs55940034))

rsq=expression(paste(r^2))

top <- egfr_f_p1 %>%
  filter(., SNP == "rs55940034")

egfr_f_p1 <- egfr_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(egfr_f_p1, aes_string(x = "GENPOS", colour = "LD_rs55940034", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=egfr_f_p1[egfr_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_eGFRcrea_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 9: regional plot of male-specific eGFR chr13 COL4A2 region, with rs7323228 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region

# load eGFR regional data
egfr_m <- read.table(here(project_dir, "regenie", "results", "step2_logeGFR2009_HRC_male_regenie_chr13_logeGFR2009.regenie"), sep = " ", header = T) %>%
  filter(., GENPOS >= min(hematuria_f$BP) & GENPOS <= max(hematuria_f$BP)) %>%
  mutate(., SNP = ID, P = 10^-LOG10P)

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs7323228 = ld$rs7323228) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs7323228 = round(LD_rs7323228^2,2))

egfr_m_p1 <- egfr_m %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs7323228))

rsq=expression(paste(r^2))

top <- egfr_m_p1 %>%
  filter(., SNP == "rs7323228")

egfr_m_p1 <- egfr_m_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(egfr_m_p1, aes_string(x = "GENPOS", colour = "LD_rs7323228", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs7323228), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=egfr_m_p1[egfr_m_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_eGFRcrea_males_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 10: regional plot of male-specific eGFR chr13 COL4A2 region, with rs9521719 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region
# no need to reload eGFR sumstats

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs9521719 = ld$rs9521719) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs9521719 = round(LD_rs9521719^2,2))

egfr_m_p2 <- egfr_m %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs9521719))

rsq=expression(paste(r^2))

top <- egfr_m_p2 %>%
  filter(., SNP == "rs9521719")

egfr_m_p2 <- egfr_m_p2 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(egfr_m_p2, aes_string(x = "GENPOS", colour = "LD_rs9521719", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs9521719), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=egfr_m_p2[egfr_m_p2$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_eGFRcrea_males_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot 11: regional plot of female-specific CAD chr13 COL4A2 region, with rs55940034 as index SNP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use same UKB LD
# use hematuria_f coordinates to restrict COL4A2 region

# load CAD regional data
cad_f <- read.table(here(project_dir, "regional_plots_data", "regional_plot_data_cad_females.txt"), sep = "\t", header = T) %>%
  filter(., BP >= min(hematuria_f$BP) & BP <= max(hematuria_f$BP)) %>%
  mutate(., P = pvalues, LOG10P = -log10(pvalues))

lead_snp <- data.frame(SNP = ld_snps$V1, LD_rs55940034 = ld$rs55940034) %>%  # filter to retain only LD for index SNP
  mutate(LD_rs55940034 = round(LD_rs55940034^2,2))

cad_f_p1 <- cad_f %>%
  left_join(., lead_snp, by = "SNP") %>%
  filter(., !is.na(LD_rs55940034))

rsq=expression(paste(r^2))

top <- cad_f_p1 %>%
  filter(., SNP == "rs55940034")

cad_f_p1 <- cad_f_p1 %>%
  mutate(is_top = case_when(
    SNP == top$SNP[1] ~ "annotate",
    SNP != top$SNP[1] ~ "none"
  ))

# plot the p-value results (Y axis; colored by LD).
ggplot(cad_f_p1, aes_string(x = "BP", colour = "LD_rs55940034", y = "LOG10P")) +
  scale_y_continuous() +  
  geom_point(colour = "black", size = 2, alpha = 7/10) +
  geom_point(aes(colour = LD_rs55940034), size = 1.5, shape = 19, alpha = 7/10) +
  #scale_colour_gradient(low = "blue", high = "red", name="LD (r2)", breaks=c(0.25, 0.5, 0.75, 1)) +
  scale_colour_stepsn(colours=c("blue2","deepskyblue","limegreen","orange2","#CC0000"), name=rsq, breaks=c(0.2, 0.4, 0.6, 0.8), labels=c("0.2","0.4","0.6","0.8")) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  guides(fill="none") + labs( x = "Chromosome 13 (GRCh37)", y = "-log10(p-value)") +
  theme(axis.text.x = element_text(face="bold", size=8), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12)) +
  geom_point(data=top, colour = "black", size = 4.4, shape = 18) + geom_point(data=top, colour="#FFFF00", shape = 18, size = 4) +
  geom_hline(yintercept=-log10(5e-08), colour = "black", linetype="dashed") +
  geom_label_repel(data=cad_f_p1[cad_f_p1$is_top=="annotate",], aes(label=as.factor(SNP)), size=3, color = "black", alpha = 0.9)
ggsave(here(project_dir, "regional_plots", stringr::str_c("regional_plot_CAD_females_with_UKB_LD_COL4A2_", top$SNP[1],".png")), width = 15, height = 15, units = "cm")






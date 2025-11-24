# Create Miami plots of TWAS with OTTERS results
# SNP, CHR, POS, PVALUE (in that order, header name does not matter)

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(R.utils)
library(gprofiler2)
library(ggrepel)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
source(here(project_dir, "scripts", "ggmirror_FLD.R"))

# Main --------------------------------------

for (tissue in c("Glom", "Tube")) {
  
  # Read summary statistics for females
  females <- fread(here(project_dir, "TWAS", "OTTERS", stringr::str_c("otters_final_pvals_", tissue, "_FEMALE.txt"))) %>%
    mutate(SNP = as.character(TargetID), BP = as.numeric(GeneStart), P = as.numeric(otters_pval_gc_excl_LASSO)) %>%
    dplyr::select(
      SNP,
      CHR = CHROM,
      BP,
      P)
  
  p_threshold_females = 0.05/nrow(females)*2
  sign_F <- females %>% filter(P < 5e-08)
  
  sign_gconvert_F <- gconvert(query = sign_F$SNP, organism = "hsapiens", target = "UCSC", mthreshold = 1) 
  
  ucsc_all_F <- sign_gconvert_F %>% 
    dplyr::select(., input, name) %>%
    right_join(., females, by = c("input" = "SNP")) %>%
    mutate(SNP = case_when(
      is.na(name) ~ input,
      name == "None" ~ input,
      TRUE ~ name
    )) %>%
    dplyr::select(
      SNP,
      CHR,
      BP,
      P
    )
  
  sign_F_plot <- sign_gconvert_F %>% dplyr::select(., name)
  
  # Read summary statistics for males
  males <- fread(here(project_dir, "TWAS", "OTTERS", stringr::str_c("otters_final_pvals_", tissue, "_MALE.txt"))) %>%
    mutate(SNP = as.character(TargetID), BP = as.numeric(GeneStart), P = as.numeric(otters_pval_gc_excl_LASSO)) %>%
    dplyr::select(
      SNP,
      CHR = CHROM,
      BP,
      P)
  
  p_threshold_males = 0.05/nrow(males)*2
  sign_M <- males %>% filter(P < 5e-08)
  
  sign_gconvert_M <- gconvert(query = sign_M$SNP, organism = "hsapiens", target = "UCSC", mthreshold = 1) 
  
  ucsc_all_M <- sign_gconvert_M %>% 
    dplyr::select(., input, name) %>%
    right_join(., males, by = c("input" = "SNP")) %>%
    mutate(SNP = case_when(
      is.na(name) ~ input,
      name == "None" ~ input,
      TRUE ~ name
    )) %>%
    dplyr::select(
      SNP,
      CHR,
      BP,
      P
    )
  
  sign_M_plot <- sign_gconvert_M %>% dplyr::select(., name)
  
  # Read summary statistics for sex-combined
  sc <- fread(here(project_dir, "TWAS", "OTTERS", stringr::str_c("otters_final_pvals_", tissue, "_sexcombined.txt"))) %>%
    mutate(SNP = as.character(TargetID), BP = as.numeric(GeneStart), P = as.numeric(otters_pval_gc_excl_LASSO)) %>%
    dplyr::select(
      SNP,
      CHR = CHROM,
      BP,
      P)
  
  p_threshold_sc = 0.05/nrow(sc)*2
  sign_SC <- sc %>% filter(P < 5e-08)
  
  sign_gconvert_SC <- gconvert(query = sign_SC$SNP, organism = "hsapiens", target = "UCSC", mthreshold = 1) 
  
  ucsc_all_SC <- sign_gconvert_SC %>% 
    dplyr::select(., input, name) %>%
    right_join(., sc, by = c("input" = "SNP")) %>%
    mutate(SNP = case_when(
      is.na(name) ~ input,
      name == "None" ~ input,
      TRUE ~ name
    )) %>%
    dplyr::select(
      SNP,
      CHR,
      BP,
      P
    )
  
  sign_SC_plot <- sign_gconvert_SC %>% dplyr::select(., name)
  
  # Read summary statistics for sex-combined METAL
  scm <- fread(here(project_dir, "TWAS", "OTTERS", stringr::str_c("otters_final_pvals_", tissue, "_sexcombined_METAL.txt"))) %>%
    mutate(SNP = as.character(TargetID), BP = as.numeric(GeneStart), P = as.numeric(otters_pval_gc_excl_LASSO)) %>%
    dplyr::select(
    SNP,
    CHR = CHROM,
    BP,
    P)
  
  p_threshold_scm = 0.05/nrow(scm)*2
  sign_SCM <- scm %>% filter(P < 5e-08)
  
  sign_gconvert_SCM <- gconvert(query = sign_SCM$SNP, organism = "hsapiens", target = "UCSC", mthreshold = 1) 
  
  ucsc_all_SCM <- sign_gconvert_SCM %>% 
    dplyr::select(., input, name) %>%
    right_join(., scm, by = c("input" = "SNP")) %>%
    mutate(SNP = case_when(
      is.na(name) ~ input,
      name == "None" ~ input,
      TRUE ~ name
    )) %>%
    dplyr::select(
      SNP,
      CHR,
      BP,
      P
    )
  
  sign_SCM_plot <- sign_gconvert_SCM %>% dplyr::select(., name)
  
  if (tissue == "Glom") {
    sign_color = "orange"
  } else if (tissue == "Tube") {
    sign_color = "purple"
  }

  # array of significant genes
  sign_twas_F_M <- c(sign_F_plot$name, sign_M_plot$name)
  sign_twas_F_SC <- c(sign_F_plot$name, sign_SC_plot$name)
  sign_twas_M_SC <- c(sign_M_plot$name, sign_SC_plot$name)
  sign_twas_SC_SCM <- c(sign_SC_plot$name, sign_SCM_plot$name)
  
  # Plot F-M
  p <- gmirror_v2(top=ucsc_all_F, bottom=ucsc_all_M, tline=p_threshold_females, bline=p_threshold_males,
                  toptitle="Females", bottomtitle="Males", 
                  highlight_p = c(p_threshold_females,p_threshold_males), highlighter=sign_color, annotate_snp = sign_twas_F_M, freey = FALSE)
  
  # Save plot
  ggsave(filename = here(project_dir, "TWAS", "OTTERS", "miami_plots", stringr::str_c("sex_specific_TWAS_NETPUNE_", tissue, "_gc_excl_LASSOSUM.png")),
         plot = p, units="in", height=7, width=12, dpi=300)
  
  # Plot F-SC
  p <- gmirror_v2(top=ucsc_all_F, bottom=ucsc_all_SC, tline=p_threshold_females, bline=p_threshold_sc,
                  toptitle="Females", bottomtitle="Sex-combined", 
                  highlight_p = c(p_threshold_females,p_threshold_sc), highlighter=sign_color, annotate_snp = sign_twas_F_SC, freey = FALSE)
  
  # Save plot
  ggsave(filename = here(project_dir, "TWAS", "OTTERS", "miami_plots", stringr::str_c("females_sex-combined_TWAS_NETPUNE_", tissue, "_gc_excl_LASSOSUM.png")),
         plot = p, units="in", height=7, width=12, dpi=300)
  
  # Plot M-SC
  p <- gmirror_v2(top=ucsc_all_M, bottom=ucsc_all_SC, tline=p_threshold_males, bline=p_threshold_sc,
                  toptitle="Males", bottomtitle="Sex-combined", 
                  highlight_p = c(p_threshold_males,p_threshold_sc), highlighter=sign_color, annotate_snp = sign_twas_M_SC, freey = FALSE)
  
  # Save plot
  ggsave(filename = here(project_dir, "TWAS", "OTTERS", "miami_plots", stringr::str_c("males_sex-combined_TWAS_NETPUNE_", tissue, "_gc_excl_LASSOSUM.png")),
         plot = p, units="in", height=7, width=12, dpi=300)
  
  # Plot SC-SCM
  p <- gmirror_v2(top=ucsc_all_SC, bottom=ucsc_all_SCM, tline=p_threshold_sc, bline=p_threshold_scm,
                  toptitle="Sex-combined", bottomtitle="Sex-combined METAL", 
                  highlight_p = c(p_threshold_sc,p_threshold_scm), highlighter=sign_color, annotate_snp = sign_twas_SC_SCM, freey = FALSE)
  
  # Save plot
  ggsave(filename = here(project_dir, "TWAS", "OTTERS","miami_plots", stringr::str_c("sex-combined_sex-combined-METAL_TWAS_NETPUNE_", tissue, "_gc_excl_LASSOSUM.png")),
         plot = p, units="in", height=7, width=12, dpi=300)
  
}

#---------------------------- Scatter plot for looking at between sex significant signals

both <- inner_join(females, males, by = c("SNP", "CHR", "BP"), suffix = c("_females", "_males")) %>%
  mutate(significance = case_when(
    P_females < p_threshold_females & P_males < p_threshold_males ~ "both",
    P_females >= p_threshold_females & P_males < p_threshold_males ~ "males",
    P_females < p_threshold_females & P_males >= p_threshold_males ~ "females",
    P_females >= p_threshold_females & P_males >= p_threshold_males ~ "none",
  ))

sign_both <- filter(both, significance %in% c("both", "females", "males"))

sign_both_plot <- gconvert(query = sign_both$SNP, organism = "hsapiens", target = "UCSC", mthreshold = 1) %>%
  dplyr::select(SNP = input, name) %>%
  left_join(., both, by = "SNP") %>%
  mutate(gene_name = case_when(
    is.na(name) ~ SNP,
    name == "None" ~ SNP,
    TRUE ~ name
  ))

ggplot(both, aes(x = -log10(P_females), y = -log10(P_males))) +
  geom_point(size = 2, alpha = 8/10, color = "black", shape = 21, aes(fill = significance)) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  labs( x = "-log10 females p-value", y = "-log10 males p-value") +
  geom_text_repel(data = sign_both_plot, aes(label = gene_name)) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12))
ggsave(here(project_dir, "TWAS", "OTTERS", "miami_plots", stringr::str_c("scatter_plot_comparison_females_males_", tissue, "_TWAS.png")), width = 20, height = 15, units = "cm")

# --------- extra to compare sex-combined with METAL and NEPTUNE eQTLs.

# scm_sign <- filter(scm, P < p_threshold_scm) %>% dplyr::select(SNP)
# sc_sign <- filter(sc, P < p_threshold_sc) %>% dplyr::select(SNP)
# 
# scm_sign_df <- filter(scm, SNP %in% c(scm_sign$SNP, sc_sign$SNP)) %>% dplyr::rename(P_SCM = P)
# sc_sign_df <- filter(sc, SNP %in% c(scm_sign$SNP, sc_sign$SNP)) %>% dplyr::rename(P_SC = P)
# 
# comp <- inner_join(scm_sign_df, sc_sign_df, by = c("SNP", "CHR", "BP"))
# 
# plot_sign <- ggplot(comp, aes(x = -log10(P_SC), y = -log10(P_SCM))) +
#   geom_point(size = 2, alpha = 5/10) +
#   theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
#   labs( x = "-log10 p-value (sex-combined NEPTUNE)", y = "-log10 p-value (sex-combined METAL)") +
#   theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
#         axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
#   theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=12))
# ggsave(here(project_dir, "comparison_P_significant_sexcombined_eqtls.png"), width = 15, height = 15, units = "cm", plot = plot_sign)

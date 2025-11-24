# analyse associations grouped by menopause status (and male comparisons with age stratification)

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(forcats)
library(vctrs)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"
cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#D55E00", "#009E73", "#0072B2", "#F0E442")

# Main --------------------------------------

top_variants = c("2:227052367", "2:227309802", "13:110415646", "16:19021066", "19:41319286", "19:41320115")

results_females <- setNames(list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_females_menopause*", full.names = T) %>%
                              stringr::str_subset(string = ., pattern = ".log", negate = TRUE),
                            nm = list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_females_menopause*", full.names = F) %>%
                              stringr::str_subset(string = ., pattern = ".log", negate = TRUE) %>%
                              stringr::str_extract(., "G[:digit:]"))

results_males <- setNames(list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_TOPMed_male_G*", full.names = T) %>%
                              stringr::str_subset(string = ., pattern = ".log", negate = TRUE),
                            nm = list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_TOPMed_male_G*", full.names = F) %>%
                              stringr::str_subset(string = ., pattern = ".log", negate = TRUE) %>%
                              stringr::str_extract(., "G[:digit:]"))
  
results_variants <- lapply(results_females, function(x) fread(x) %>%
                               mutate(CHR_POS = stringr::str_c(CHROM, ":", GENPOS)) %>%
                               filter(CHR_POS %in% top_variants) %>%
                               mutate(rsid = 
                                        case_when(
                                          CHR_POS == "2:227052367" ~ "rs35138315",
                                          CHR_POS == "2:227309802" ~ "rs4305262",
                                          CHR_POS == "13:110415646" ~ "rs1927355",
                                          CHR_POS == "16:19021066" ~ "rs12919754",
                                          CHR_POS == "19:41319286" ~ "rs73045269",
                                          CHR_POS == "19:41320115" ~ "rs56254331"
                                        )))
  
results_df <- rbindlist(results_variants, idcol = TRUE) %>%
  mutate(chr_rsid = stringr::str_c("chr", CHROM, ": ", rsid),
         menopause_status = case_when(
           .id == "G1" ~ 0,
           .id == "G2" ~ 1
         ),
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr16: rs12919754", "chr19: rs56254331", "chr19: rs73045269"))) %>%
  dplyr::select(-N_CASES, -N_CONTROLS)

# add results non-stratified
NS_results_list <- list.files(here(project_dir, "regenie", "results"), pattern = "X593_TOPMED_female_regenie_*", full.names = T) %>%
  stringr::str_subset(string = ., pattern = ".log", negate = TRUE)

NS_results_variants <- lapply(NS_results_list, function(x) fread(x) %>%
                                mutate(CHR_POS = stringr::str_c(CHROM, ":", GENPOS)) %>%
                                filter(CHR_POS %in% top_variants) %>%
                                filter(CHR_POS != "19:41320115") %>%
                                mutate(rsid = 
                                         case_when(
                                           CHR_POS == "2:227052367" ~ "rs35138315",
                                           CHR_POS == "2:227309802" ~ "rs4305262",
                                           CHR_POS == "13:110415646" ~ "rs1927355",
                                           CHR_POS == "16:19021066" ~ "rs12919754",
                                           CHR_POS == "19:41319286" ~ "rs73045269",
                                           CHR_POS == "19:41320115" ~ "rs56254331"
                                         )))

NS_results_df <- rbindlist(NS_results_variants) %>%
  mutate(chr_rsid = stringr::str_c("chr", CHROM, ": ", rsid),
         menopause_status = "0+1",
         .id = "G0",
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr16: rs12919754", "chr19: rs56254331", "chr19: rs73045269"))
  )

all <- rbind(results_df, NS_results_df) %>%
  mutate(ord_menopause_status = factor(menopause_status, levels = c("1", "0", "0+1"))) %>%
  dplyr::select(ord_chr_rsid, BETA, SE, ord_menopause_status)

# add info (non stratified) for rs4305262:
miss_info <- data.frame(ord_chr_rsid = "chr2: rs4305262", BETA = -0.103344, SE = 0.018326, ord_menopause_status = "0+1")

all <- rbind(all, miss_info)

# write.table(all, here(project_dir, "age_sex_analyses", "menopause_status_beta_TOPMed_variants_table.txt"), sep = "\t", row.names = F, quote = F)

# test for differential effects between menopause groups:
menopause_0 <- all %>% filter(ord_menopause_status == 0) %>%
  dplyr::select(ord_chr_rsid, BETA_menopause0 = BETA, SE_menopause0 = SE)
menopause_1 <- all %>% filter(ord_menopause_status == 1) %>%
  dplyr::select(ord_chr_rsid, BETA_menopause1 = BETA, SE_menopause1 = SE)

DE <- inner_join(menopause_0, menopause_1, by = "ord_chr_rsid") %>%
  mutate(zscore = (BETA_menopause0 - BETA_menopause1)/(sqrt(SE_menopause0^2 + SE_menopause1^2)),
         zscore_pval = 2*pnorm(q=abs(zscore), lower.tail = FALSE),
         fdr = p.adjust(zscore_pval, method = "BH"),
         significance = case_when(zscore_pval < 0.05 ~ "*",
                                  zscore_pval >= 0.05 ~ ""))

# plot 

all <- all %>%
  left_join(., DE, by = "ord_chr_rsid") %>%
  mutate(dummy_sign = case_when(
    significance == "*" & ord_menopause_status == "0+1" ~ "*",
    significance == "*" & ord_menopause_status != "0+1" ~ "",
    significance == "" ~ ""
  )) %>%
  dplyr::select(ord_chr_rsid, BETA, SE, ord_menopause_status, dummy_sign)

sign <- all %>%
  filter(., dummy_sign == "*")

plot_f <- ggplot(all, aes(x=BETA, y=fct_rev(ord_chr_rsid), color = as.factor(ord_menopause_status))) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") + 
  # Add dot plot and error bars
  geom_errorbar(aes(xmin = BETA-SE, xmax = BETA+SE), width = 0.25, size = 1, position=position_dodge(width=0.5)) +
  geom_point(size = 3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=cbbPalette, breaks=c("0+1", "0", "1"), labels = c("no stratification", "pre-menopause", "post-menopause")) +
  labs(x="Effect size", y = "", title = "A") +
  # scale_x_continuous(breaks=seq(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI)),), limits=c(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI))), expand=c(0,0) ) +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, color = "black"), axis.text.y = element_text(size=14, color = "black"),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15, angle = )) +
  theme(plot.title = element_text(hjust =0),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  guides(color = guide_legend(title = "Menopause status"))

ggsave(here(project_dir, "age_sex_analyses", "menopause_status_beta_TOPMed_variants_updated.jpg"),
       plot = plot_f,
       width = 9, height = 6, dpi = 300)

## -------------------------------- compare with male results (using age groups)

results_males <- setNames(list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_TOPMed_male_G*", full.names = T) %>%
                            stringr::str_subset(string = ., pattern = ".log", negate = TRUE),
                          nm = list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_TOPMed_male_G*", full.names = F) %>%
                            stringr::str_subset(string = ., pattern = ".log", negate = TRUE) %>%
                            stringr::str_extract(., "G[:digit:]"))

results_variants <- lapply(results_males, function(x) fread(x) %>%
                             mutate(CHR_POS = stringr::str_c(CHROM, ":", GENPOS)) %>%
                             filter(CHR_POS %in% top_variants) %>%
                             mutate(rsid = 
                                      case_when(
                                        CHR_POS == "2:227052367" ~ "rs35138315",
                                        CHR_POS == "2:227309802" ~ "rs4305262",
                                        CHR_POS == "13:110415646" ~ "rs1927355",
                                        CHR_POS == "19:41319286" ~ "rs73045269",
                                        CHR_POS == "19:41320115" ~ "rs56254331"
                                      )))

results_df <- rbindlist(results_variants, idcol = TRUE) %>%
  mutate(chr_rsid = stringr::str_c("chr", CHROM, ": ", rsid),
         agegroup_status = case_when(
           .id == "G1" ~ 1,
           .id == "G2" ~ 2
         ),
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr19: rs56254331", "chr19: rs73045269")))

# add results non-stratified
NS_results_list <- list.files(here(project_dir, "regenie", "results"), pattern = "X593_TOPMED_male_regenie_*", full.names = T) %>%
  stringr::str_subset(string = ., pattern = ".log", negate = TRUE)

NS_results_variants <- lapply(NS_results_list, function(x) fread(x) %>%
                                mutate(CHR_POS = stringr::str_c(CHROM, ":", GENPOS)) %>%
                                filter(CHR_POS %in% top_variants) %>%
                                filter(CHR_POS != "19:41319286") %>%
                                mutate(rsid = 
                                         case_when(
                                           CHR_POS == "2:227052367" ~ "rs35138315",
                                           CHR_POS == "2:227309802" ~ "rs4305262",
                                           CHR_POS == "13:110415646" ~ "rs1927355",
                                           CHR_POS == "19:41319286" ~ "rs73045269",
                                           CHR_POS == "19:41320115" ~ "rs56254331"
                                         )))

NS_results_df <- rbindlist(NS_results_variants) %>%
  mutate(chr_rsid = stringr::str_c("chr", CHROM, ": ", rsid),
         agegroup_status = "1+2",
         .id = "G0",
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr19: rs56254331", "chr19: rs73045269"))
  )

all <- rbind(results_df, NS_results_df) %>%
  mutate(ord_agegroup_status = factor(agegroup_status, levels = c("2", "1", "1+2"))) %>%
  dplyr::select(ord_chr_rsid, BETA, SE, ord_agegroup_status)

# add info (non stratified) for rs4305262 (2:227309802):
miss_info <- data.frame(ord_chr_rsid = "chr2: rs4305262", BETA = -0.028, SE = 0.01582, ord_agegroup_status = "1+2")

all <- rbind(all, miss_info)

# write.table(all, here(project_dir, "age_sex_analyses", "male_age_groups_status_beta_TOPMed_variants_table.txt"), sep = "\t", row.names = F, quote = F)

# test for differential effects between menopause groups:
agegroup_1 <- all %>% filter(ord_agegroup_status == 1) %>%
  dplyr::select(ord_chr_rsid, BETA_agegroup1 = BETA, SE_agegroup1 = SE)
agegroup_2 <- all %>% filter(ord_agegroup_status == 2) %>%
  dplyr::select(ord_chr_rsid, BETA_agegroup2 = BETA, SE_agegroup2 = SE)

DE <- inner_join(agegroup_1, agegroup_2, by = "ord_chr_rsid") %>%
  mutate(zscore = (BETA_agegroup1 - BETA_agegroup2)/(sqrt(SE_agegroup1^2 + SE_agegroup2^2)),
         zscore_pval = 2*pnorm(q=abs(zscore), lower.tail = FALSE),
         fdr = p.adjust(zscore_pval, method = "BH"),
         significance = case_when(zscore_pval < 0.05 ~ "*",
                                  zscore_pval >= 0.05 ~ ""))

all_out <- all %>%
  left_join(., DE, by = "ord_chr_rsid") %>%
  mutate(dummy_sign = case_when(
    significance == "*" & ord_agegroup_status == "1+2" ~ "*",
    significance == "*" & ord_agegroup_status != "1+2" ~ "",
    significance == "" ~ ""
  )) %>%
  dplyr::select(ord_chr_rsid, BETA, SE, ord_agegroup_status, dummy_sign)

sign <- all_out %>%
  filter(., dummy_sign == "*")

# plot 
plot_m <- ggplot(all_out, aes(x=BETA, y=fct_rev(ord_chr_rsid), color = as.factor(ord_agegroup_status))) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") + 
  # Add dot plot and error bars
  geom_errorbar(aes(xmin = BETA-SE, xmax = BETA+SE), width = 0.25, size = 1, position=position_dodge(width=0.5)) +
  geom_point(size = 3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=cbbPalette, breaks = c("1+2", "1", "2"), labels = c("no stratification", "< 50 years of age", "\u2265 50 years of age")) +
  labs(x="Effect size", y = "", title = "B") +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, color = "black"), axis.text.y = element_text(size=14, color = "black"),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15, angle = )) +
  theme(plot.title = element_text(hjust =0),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  guides(color = guide_legend(title = "Male age groups")) +
  geom_text(data=sign, colour = "black", size = 12, aes(y = ord_chr_rsid, x = BETA+SE+0.5, label = dummy_sign))

ggsave(here(project_dir, "age_sex_analyses", "male_age_groups_status_beta_TOPMed_variants_updated.jpg"),
       plot = plot_m,
       width = 9, height = 6, dpi = 300)

figure_all <- grid.arrange(plot_f, plot_m, nrow = 1, widths = c(6, 6))

ggsave(here(project_dir, "age_sex_analyses", "female_menopause_male_age_groups_beta_TOPMed_variants_updated.jpg"),
       plot = figure_all,
       width = 15, height = 6, dpi = 300)

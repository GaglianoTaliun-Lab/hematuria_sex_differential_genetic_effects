# analyse associations grouped by menopause status

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

top_variants = c("2:227052367", "2:227309802", "13:110415646", "19:41319286", "19:41320115")

results_list <- setNames(list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_females_menopause*", full.names = T) %>%
  stringr::str_subset(string = ., pattern = ".log", negate = "TRUE"),
  nm = list.files(here(project_dir, "age_sex_analyses", "TOPMed_regenie_results"), pattern = "step2_ukb_X593_females_menopause*", full.names = F) %>%
    stringr::str_subset(string = ., pattern = ".log", negate = "TRUE") %>%
    stringr::str_extract(., "G[:digit:]"))

results_variants <- lapply(results_list, function(x) fread(x) %>%
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
         menopause_status = case_when(
           .id == "G1" ~ 0,
           .id == "G2" ~ 1
         ),
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr19: rs56254331", "chr19: rs73045269")))

# add results non-stratified
NS_results_list <- list.files(here(project_dir, "regenie", "results"), pattern = "X593_TOPMED_female_regenie_*", full.names = T) %>%
  stringr::str_subset(string = ., pattern = ".log", negate = "TRUE")

NS_results_variants <- lapply(NS_results_list, function(x) fread(x) %>%
                                mutate(CHR_POS = stringr::str_c(CHROM, ":", GENPOS)) %>%
                                filter(CHR_POS %in% top_variants) %>%
                                filter(CHR_POS != "19:41320115") %>%
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
         menopause_status = "0+1",
         .id = "G0",
         ord_chr_rsid = factor(chr_rsid, levels = c("chr2: rs35138315", "chr2: rs4305262","chr13: rs1927355", "chr19: rs56254331", "chr19: rs73045269"))
  )

all <- rbind(results_df, NS_results_df) %>%
  mutate(ord_menopause_status = factor(menopause_status, levels = c("0+1", 0, 1))) %>%
  dplyr::select(ord_chr_rsid, BETA, SE, ord_menopause_status)

# add info (non stratified) for rs4305262:
miss_info <- data.frame(ord_chr_rsid = "chr2: rs4305262", BETA = -0.103344, SE = 0.018326, ord_menopause_status = "0+1")

all <- rbind(all, miss_info)

# plot 
ggplot(all, aes(x = ord_chr_rsid, y = BETA, color = as.factor(ord_menopause_status))) +
  geom_point(position=position_dodge(width=0.3)) +
  geom_errorbar(aes(ymin=BETA-SE, ymax=BETA+SE), width = 0.1, position=position_dodge(width=0.3)) +
  scale_color_manual(values=cbbPalette) +
  labs(x="", y = "Effect size", title = "Association of lead variants") +
  theme(plot.title = element_text(hjust =0.5),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.length=unit(0.3,"cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
  ) + theme_classic(base_size=14) +
  theme(axis.text.x = element_text(face="bold", size=12, angle = 90, vjust = 0.5), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14)) +
  guides(color = guide_legend(title = "Menopause status"))

ggsave(here(project_dir, "age_sex_analyses", "menopause_status_beta_TOPMed_variants.jpg"), width=10, height=7)


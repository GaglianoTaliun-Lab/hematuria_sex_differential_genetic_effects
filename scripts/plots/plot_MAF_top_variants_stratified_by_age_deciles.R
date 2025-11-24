# create figures for MAF by age deciles and controls
# for cases, the age used was the smallest age at diagnosis (i.e., oldest time stamp)
# for rs35138315 the freqs were calculated in DNA Nexus with UKB TOPMed imputed variants.
# for other variants, freqs were calculated using the UKB HRC-imputed variants.

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

# Read tables and format for plots --------------------------------------

males <- setNames(object = list.files(here(project_dir, "age_sex_analyses", "MAF_age_decile_groups_results"), pattern = "*.afreq", full.names = T) %>%
                    stringr::str_subset(., "_male.afreq"), 
                      nm = list.files(here(project_dir, "age_sex_analyses", "MAF_age_decile_groups_results"), pattern = "*.afreq", full.names = F) %>%
                        stringr::str_subset(., "_male.afreq") %>%
                        stringr::str_remove(., "_male.afreq") %>%
                        stringr::str_remove(., "freq_") %>%
                        stringr::str_remove(., "age_"))

df_males <- lapply(males, function(x) fread(x)) %>%
  rbindlist(., idcol = TRUE) %>%
  separate(.id, c("group", "variant_id"), remove = TRUE) %>% 
  mutate(freq_SE = sqrt((ALT_FREQS*(1-ALT_FREQS))/OBS_CT),
         sex = "Males",
         group = case_when(
           group == "all" ~ "Overall",
           group == "cases" ~ "Cases (all)",
           group == "controls" ~ "Controls",
           group == "decile1" ~ "Decile 1",
           group == "decile2" ~ "Decile 2",
           group == "decile3" ~ "Decile 3",
           group == "decile4" ~ "Decile 4",
           group == "decile5" ~ "Decile 5",
           group == "decile6" ~ "Decile 6",
           group == "decile7" ~ "Decile 7",
           group == "decile8" ~ "Decile 8",
           group == "decile9" ~ "Decile 9",
           group == "decile10" ~ "Decile 10",
         ),
         ord_group = factor(group, levels = c("Overall", "Controls", "Cases (all)", "Decile 1", "Decile 2", "Decile 3", "Decile 4", "Decile 5", 
                                              "Decile 6","Decile 7", "Decile 8", "Decile 9", "Decile 10"))) %>%
  filter(., group != "Overall") %>%
  dplyr::select(group = ord_group,
                sex,
                variant_id,
                freq = ALT_FREQS,
                freq_SE,
                N_alleles = OBS_CT)
  

females <- setNames(object = list.files(here(project_dir, "age_sex_analyses", "MAF_age_decile_groups_results"), pattern = "*.afreq", full.names = T) %>%
                               stringr::str_subset(., "_female.afreq"), 
                         nm = list.files(here(project_dir, "age_sex_analyses", "MAF_age_decile_groups_results"), pattern = "*.afreq", full.names = F) %>%
                           stringr::str_subset(., "_female.afreq") %>%
                           stringr::str_remove(., "_female.afreq") %>%
                           stringr::str_remove(., "freq_") %>%
                           stringr::str_remove(., "age_"))

df_females <- lapply(females, function(x) fread(x)) %>%
  rbindlist(., idcol = TRUE) %>%
  separate(.id, c("group", "variant_id"), remove = TRUE) %>% 
  mutate(freq_SE = sqrt((ALT_FREQS*(1-ALT_FREQS))/OBS_CT),
         sex = "Females",
         group = case_when(
           group == "all" ~ "Overall",
           group == "cases" ~ "Cases (all)",
           group == "controls" ~ "Controls",
           group == "decile1" ~ "Decile 1",
           group == "decile2" ~ "Decile 2",
           group == "decile3" ~ "Decile 3",
           group == "decile4" ~ "Decile 4",
           group == "decile5" ~ "Decile 5",
           group == "decile6" ~ "Decile 6",
           group == "decile7" ~ "Decile 7",
           group == "decile8" ~ "Decile 8",
           group == "decile9" ~ "Decile 9",
           group == "decile10" ~ "Decile 10",
         ),
         ord_group = factor(group, levels = c("Overall", "Controls", "Cases (all)", "Decile 1", "Decile 2", "Decile 3", "Decile 4", "Decile 5", 
                                              "Decile 6","Decile 7", "Decile 8", "Decile 9", "Decile 10"))) %>%
  filter(., group != "Overall") %>%
  dplyr::select(group = ord_group,
                sex,
                variant_id,
                freq = ALT_FREQS,
                freq_SE,
                N_alleles = OBS_CT)

list_by_variant <- rbind(df_males, df_females) %>%
  dplyr::group_split(., variant_id)

# get age range for each decile -------------------------

inds_deciles <- setNames(object = list.files(here(project_dir, "age_sex_analyses", "age_decile_groups"), pattern = "ukb_inds_IDs_age*", full.names = TRUE),
                         nm = list.files(here(project_dir, "age_sex_analyses", "age_decile_groups"), pattern = "ukb_inds_IDs_age*", full.names = F) %>%
                                           stringr::str_remove(., "ukb_inds_IDs_age_") %>%
                                           stringr::str_remove(., ".txt"))

age_diagnosis <- read.table(here(project_dir, "age_sex_analyses", "age_at_diagnosis_cases_X593.txt"), sep = "\t", header = TRUE) %>%
  dplyr::select(FID, IID, min_age_at_diagnosis)

age_ranges_df <- data.frame("sex" = c("female", "male"), "decile1" = c(NA, NA), "decile2"= c(NA, NA), "decile3"= c(NA, NA), "decile4"= c(NA, NA),
                            "decile5"= c(NA, NA), "decile6"= c(NA, NA), "decile7"= c(NA, NA), "decile8"= c(NA, NA), "decile9"= c(NA, NA), "decile10"= c(NA, NA)) %>%
  tibble::column_to_rownames(., var = "sex")

for (i in 1:length(inds_deciles)) {
  
  ids_decile <- read.table(inds_deciles[i], sep = "\t", header = TRUE)
  decile_sex = names(inds_deciles[i]) %>% str_split(., "_", 2)
  
  age_group_filtered <- filter(age_diagnosis, FID %in% ids_decile$FID)
  min_age = min(age_group_filtered$min_age_at_diagnosis)
  max_age = max(age_group_filtered$min_age_at_diagnosis)
  age_range = stringr::str_c("[", min_age, " - ", max_age, "]")
  
  age_ranges_df[decile_sex[[1]][2], decile_sex[[1]][1]] <- age_range
  
}

colnames(age_ranges_df) <- c("Decile 1", "Decile 2", "Decile 3", "Decile 4", "Decile 5", 
                              "Decile 6","Decile 7", "Decile 8", "Decile 9", "Decile 10")

# dictionary to include gene name:
gene_names = data.frame(variant_id = c("rs1927355", "rs35138315", "rs4305262", "rs56254331", "rs12919754", "rs111479366", "rs755713332"), gene = c("COL4A2", "COL4A4", "COL4A3", "TGFB1", "TMC7", "TMC7", "ZDHHC15"))

# Wrangle and plot --------------------------------------

for (i in 1:length(list_by_variant)) {
  
  df_plot = list_by_variant[[i]] %>%
    left_join(., gene_names, by = "variant_id")
  rsid=df_plot$variant_id[1]
  gene=df_plot$gene[1]
  
  # table for N alleles:
  mytable <- dplyr::select(df_plot, group, sex, N_alleles) %>%
    tidyr::spread(., group, N_alleles) %>%
    tableGrob(., rows=NULL)
  
  # table for age ranges:
  mytable2 <- age_ranges_df %>% 
    tibble::rownames_to_column(., "sex") %>%
    filter(., sex == "male") %>%
    tibble::column_to_rownames(., "sex") %>%
    tableGrob(., rows=NULL)
  
  title1 <- textGrob("Age ranges",gp=gpar(fontsize=10, fontface = "bold"))
  title2 <- textGrob("N alleles",gp=gpar(fontsize=10, fontface = "bold"))
  
  plot1 <- ggplot(df_plot, aes(x = as.factor(group), y = freq, colour = as.factor(sex))) +
    geom_point(position=position_dodge(width=0.3), size = 3, aes(shape = as.factor(sex))) +
    geom_errorbar(aes(ymin=freq-freq_SE, ymax=freq+freq_SE), width = 0.3, size = 0.9, position=position_dodge(width=0.3)) +
    labs(x="Hematuria age deciles and controls", y = stringr::str_c("EAF ", rsid, " (", gene, ")"), title = "") +
    theme(plot.title = element_text(hjust =0.5),
          axis.line.x = element_line(linewidth = 0.6),
          axis.ticks.length=unit(0.3,"cm"),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
    ) + theme_classic(base_size=14) +
    theme(axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, color = "black"), axis.text.y = element_text(size=14, color = "black"),
          axis.title.x = element_text(face="bold", size=16), axis.title.y = element_text(face="bold", size=16)) +
    theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) +
    scale_shape_manual(values = c(16,15)) +
    guides(shape = "none") +
    scale_colour_manual(values = c("Females" = "#DDAA33", "Males" = "#BB5566"),
                        name = "Sex",
                        guide = guide_legend(override.aes = list(shape = c(16,15), size = 4, stroke = 1.5, colour = c("#DDAA33", "#BB5566"))))
  

  figure <- grid.arrange(plot1, title1, mytable2, title2, mytable,
               nrow = 5, heights = c(5, 0.3, 0.4, 0.2, 0.6), clip = FALSE)
  
  ggsave(here(project_dir, "age_sex_analyses", "MAF_age_decile_groups_results", "plots", stringr::str_c("plot_MAF_by_age_deciles_sex_", rsid, ".png")),
         plot = figure,
         width = 14, height = 10, dpi = 300)
  
}

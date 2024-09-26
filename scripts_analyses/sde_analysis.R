# Description: look at SDE results

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sde <- read.table(here(project_dir, "SDE", "SDE_Z_scores_significant_1e-05.tsv"), sep = "\t", header = T)

# Create a list for SNPnexus:
data.frame(dbsnp = "dbsnp", rsid = sde$SNP) %>%
  write.table(., here(project_dir, "SDE", "input_for_SNPNexus.txt"), sep = "\t", row.names = F, col.names = F, quote = F)

# Load annotation of SNPs with SNPNexus:
snp_annot <- read.table(here(project_dir, "SDE", "snp_annotation_snpnexus.txt"), sep = "\t", header = T) %>%
  mutate(
    gene_annotation = 
           case_when(
             !is.na(Overlapped.Gene) ~ Overlapped.Gene,
             is.na(Overlapped.Gene) & Distance.to.Nearest.Upstream.Gene < Distance.to.Nearest.Downstream.Gene ~ Nearest.Upstream.Gene,
             is.na(Overlapped.Gene) & Distance.to.Nearest.Upstream.Gene > Distance.to.Nearest.Downstream.Gene ~ Nearest.Downstream.Gene,
             TRUE ~ stringr::str_c("Upstream: ", Nearest.Upstream.Gene, "; Downstream: ", Nearest.Downstream.Gene))
             ) %>%
  dplyr::select(SNP = Variation.ID,
                gene_annotation)

sde <- sde %>%
  left_join(., snp_annot, by = "SNP")

sde %>% 
  filter(., z_score_pval <= 1e-06) %>%
  arrange(., CHR, BP) %>%
  write.table(., here(project_dir, "SDE", "SDE_significant_signals_1e-06.txt"), sep = "\t", row.names = F, quote = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Volcano plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# p-value threshold for plotting:
thres=0.05

volcano_z <- sde %>%
  mutate(males = BETA_M/SE_M,
         females = BETA_F/SE_F) %>%
  gather(., "sex", "Zscore", c("males", "females")) %>%
  dplyr::select(SNP, sex, Zscore)

volcano_p <- sde %>%
  mutate(males = -log10(P_M),
         females = -log10(P_F)) %>%
  gather(., "sex", "pvalue", c("males", "females")) %>%
  dplyr::select(SNP, sex, pvalue)

sde_volcano <- inner_join(volcano_z, volcano_p, by = c("SNP", "sex")) %>%
  filter(., pvalue >= -log10(thres))

ggplot(data = sde_volcano, aes(x = Zscore, y = pvalue)) +
  geom_point(aes(colour = sex)) + geom_line(aes(group=SNP), color="grey", 
                                            size=0.1, alpha=0.8) +
  labs( x = "Z-score", y = "-log10(p-value)") +
  geom_hline(yintercept = -log10(5e-08), linetype = 2, color = "blueviolet") +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank())
ggsave(here(project_dir, "SDE", stringr::str_c("volcano_plot_significant_", thres, ".png")))

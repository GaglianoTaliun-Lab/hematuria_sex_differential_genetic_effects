# plot partitioned h2

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

# Read table --------------------------------------
h2 <- read.table(here(project_dir, "bolt-lmm", "imputed_pruned_1Mb", "summary_partitioned_h2.txt"), sep = "\t", header = T)
 
# Plot --------------------------------------
ggplot(h2, aes(x = as.factor(locus), y = h2, color = sex)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=h2-se,ymax=h2+se), width = 0.3)  +
  labs(x="Locus in chromosome", y = "SNP h2") +
  theme(legend.title = element_blank()) +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=14))
ggsave(here(project_dir, "bolt-lmm", "imputed_pruned_1Mb", "partitioned_h2_X593.png"))

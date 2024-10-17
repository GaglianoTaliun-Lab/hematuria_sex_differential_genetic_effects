# Description: run hyprcoloc for BP traits, CAD and hematuria (females) for the COL4A1/2 region

# Packages -------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(hyprcoloc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/col4a2_hematuria"

# COL4A2 GRCh37 coordinates: 13:110,959,631-111,165,556
col4a2_start=110959631
col4a2_end=111165556

# COL4A1 GRCh37 coordinates: 13: 110,801,310-110,959,504
# col4a1_start=110801310
# col4a1_end=110959504

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read datasets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# restrict to COL4A2 which is the region of interest in hematuria:

cad <- fread(here(project_dir, "data_for_coloc", "cad_females.txt")) %>%
    mutate(se = sqrt(varbeta)) %>%
    filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
    filter(., !is.na(beta)) %>%
    filter(., !is.na(se)) %>%
    distinct(., SNP, .keep_all = TRUE)

dbp <- fread(here(project_dir, "data_for_coloc", "DBP_females.txt")) %>%
    mutate(se = sqrt(varbeta)) %>%
    filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
    filter(., !is.na(beta)) %>%
    filter(., !is.na(se)) %>%
    distinct(., SNP, .keep_all = TRUE)

pp <- fread(here(project_dir, "data_for_coloc", "PP_females.txt")) %>%
    mutate(se = sqrt(varbeta)) %>%
    filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
    filter(., !is.na(beta)) %>%
    filter(., !is.na(se)) %>%
    distinct(., SNP, .keep_all = TRUE)

hematuria <- fread(here(project_dir, "data_for_coloc", "hematuria_females.txt")) %>%
    mutate(se = sqrt(varbeta)) %>%
    filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
    filter(., !is.na(beta)) %>%
    filter(., !is.na(se)) %>%
    distinct(., SNP, .keep_all = TRUE)

# Read LD matrix:
ld <- fread(here(project_dir, "data_for_coloc", "LD_COL4A2_r_matrix_hyprcoloc.phased.vcor1"))
ld_vars <- read.table(here(project_dir, "data_for_coloc", "LD_COL4A2_r_matrix_hyprcoloc.phased.vcor1.vars"), col.names = F, sep = "\t")
colnames(ld) <- ld_vars$V1
rownames(ld) <- ld_vars$V1

ld <- as.matrix(ld)

# Trait correlations (https://stackoverflow.com/questions/57903948/creating-a-correlation-matrix-from-a-data-frame-in-r):
cor_traits <- fread(here(project_dir, "results_ldsc", "ldsc_correlations.txt")) %>%
    dplyr::select(p1,p2,rg)

data1 <- data.frame(p1 = cor_traits$p1, p2 = cor_traits$p2, rg = cor_traits$rg)  
df <- rbind(cor_traits, data1)
cor_matrix <- as.data.frame.matrix(xtabs(rg ~ ., df))
diag(cor_matrix) <- 1

cor_matrix <- as.matrix(cor_matrix)

# Sample overlaps:
overlap_matrix <- read.table(here(project_dir, "results_ldsc", "sample_overlap.txt"), sep = "\t", header = T, row.names = 1) %>%
    as.matrix()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data for hyprcoloc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# obtain common SNPs to test
common_snps <- inner_join(cad, dbp, by = "SNP") %>%
    inner_join(., pp, by = "SNP") %>%
    inner_join(., hematuria, by = "SNP") %>%
    dplyr::select(SNP)

write.table(common_snps$SNP, here(project_dir, "data_for_coloc", "hyprcoloc_SNPs_for_LDmatrix.txt"), sep = "\t", row.names = F, quote = F, col.names = F)

cad <- cad %>%
    filter(., SNP %in% common_snps$SNP)

dbp <- dbp %>%
    filter(., SNP %in% common_snps$SNP)

pp <- pp %>%
    filter(., SNP %in% common_snps$SNP)

hematuria <- hematuria %>%
    filter(., SNP %in% common_snps$SNP)

betas <- as.matrix(data.frame(CAD = cad$beta, DBP = dbp$beta, hematuria = hematuria$beta, PP = pp$beta, row.names = common_snps$SNP))
ses <- as.matrix(data.frame(CAD = cad$se, DBP = dbp$se, hematuria = hematuria$se, PP = pp$se, row.names = common_snps$SNP))
traits=c("CAD", "DBP", "hematuria", "PP")
rsid=common_snps$SNP

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run hyprcoloc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prior.1 = the ratio of the prior probability of a hypothesis in which at least one trait has a causal variant in the region by the null hypothesis that none have a causal variant
# prior.c = the probability of a snp being associated with a second trait given its association with one trait
    # if 0.05, then there is 1 in 20 chances; if 0.01, then there is 1 in 100 chances, etc. 

hyprcoloc_res <- hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid, 
                        binary.outcomes = c(1,0,1,0), ld.matrix = ld, trait.cor = cor_matrix, sample.overlap = overlap_matrix, uniform.priors = FALSE,
                        prior.1 = 1e-04, reg.steps = 2)

write.table(hyprcoloc_res$results, here(project_dir, "results_coloc", "hyprcoloc_results_p1_1e-04.txt"), sep = "\t", row.names = F, quote = F)

hyprcoloc_res_strictp1 <- hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid, 
                        binary.outcomes = c(1,0,1,0), ld.matrix = ld, trait.cor = cor_matrix, sample.overlap = overlap_matrix, uniform.priors = FALSE,
                        prior.1 = 1e-04, reg.steps = 2, prior.c = 0.02)

write.table(hyprcoloc_res_strictp1$results, here(project_dir, "results_coloc", "hyprcoloc_results_p1_1e-05_pc_0.02.txt"), sep = "\t", row.names = F, quote = F)

# sensitivity plot
pdf(here(project_dir, "results_coloc", "hyprcoloc_sensitivity_plot.pdf"))
sensitivity.plot(betas, ses, trait.names = traits, snp.id = rsid, 
                binary.outcomes = c(1,0,1,0), ld.matrix = ld, trait.cor = cor_matrix, sample.overlap = overlap_matrix, uniform.priors = FALSE,
                reg.thresh = c(0.6, 0.7, 0.8, 0.9), align.thresh = c(0.6, 0.7, 0.8, 0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE)
dev.off()
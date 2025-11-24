# Description: run hyprcoloc for BP traits, CAD and hematuria (females and males) for the COL4A1/2 region

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

project_dir <- "/home/fridald4/links/projects/def-gsarah/fridald4/col4a2_hematuria"

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

for (sex in c("females", "males")) {

    cad <- fread(here(project_dir, "data_for_coloc", stringr::str_c("cad_", sex, ".txt"))) %>%
        mutate(se = sqrt(varbeta)) %>%
        filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
        filter(., !is.na(beta)) %>%
        filter(., !is.na(se)) %>%
        distinct(., SNP, .keep_all = TRUE)

    dbp <- fread(here(project_dir, "data_for_coloc", stringr::str_c("DBP_", sex, ".txt"))) %>%
        mutate(se = sqrt(varbeta)) %>%
        filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
        filter(., !is.na(beta)) %>%
        filter(., !is.na(se)) %>%
        distinct(., SNP, .keep_all = TRUE)

    pp <- fread(here(project_dir, "data_for_coloc", stringr::str_c("PP_", sex, ".txt"))) %>%
        mutate(se = sqrt(varbeta)) %>%
        filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
        filter(., !is.na(beta)) %>%
        filter(., !is.na(se)) %>%
        distinct(., SNP, .keep_all = TRUE)

    hematuria <- fread(here(project_dir, "data_for_coloc", stringr::str_c("hematuria_", sex, ".txt"))) %>%
        mutate(se = sqrt(varbeta)) %>%
        filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
        filter(., !is.na(beta)) %>%
        filter(., !is.na(se)) %>%
        distinct(., SNP, .keep_all = TRUE)

    egfr <- fread(here(project_dir, "data_for_coloc", stringr::str_c("eGFRcrea_", sex, ".txt"))) %>%
        mutate(se = sqrt(varbeta)) %>%
        filter(., BP >= col4a2_start & BP <= col4a2_end) %>%
        filter(., !is.na(beta)) %>%
        filter(., !is.na(se)) %>%
        distinct(., SNP, .keep_all = TRUE)

    # Read LD matrix:
    # ld <- fread(here(project_dir, "data_for_coloc", "LD_COL4A2_r_matrix_hyprcoloc.phased.vcor1"))
    # ld_vars <- read.table(here(project_dir, "data_for_coloc", "LD_COL4A2_r_matrix_hyprcoloc.phased.vcor1.vars"), col.names = F, sep = "\t")
    # colnames(ld) <- ld_vars$V1
    # rownames(ld) <- ld_vars$V1

    # ld <- as.matrix(ld)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Prepare data for hyprcoloc
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # obtain common SNPs to test
    common_snps <- inner_join(cad, dbp, by = "SNP") %>%
        inner_join(., pp, by = "SNP") %>%
        inner_join(., hematuria, by = "SNP") %>%
        inner_join(., egfr, by = "SNP") %>%
        dplyr::select(SNP)

    # write.table(common_snps$SNP, here(project_dir, "data_for_coloc", stringr::str_c("hyprcoloc_SNPs_for_LDmatrix_", sex, ".txt")), sep = "\t", row.names = F, quote = F, col.names = F)

    cad <- cad %>%
        filter(., SNP %in% common_snps$SNP)

    dbp <- dbp %>%
        filter(., SNP %in% common_snps$SNP)

    pp <- pp %>%
        filter(., SNP %in% common_snps$SNP)

    hematuria <- hematuria %>%
        filter(., SNP %in% common_snps$SNP)

    egfr <- egfr %>%
        filter(., SNP %in% common_snps$SNP)

    if (sex == "females") {
        betas <- as.matrix(data.frame(CAD_females = cad$beta, DBP_females = dbp$beta, hematuria_females = hematuria$beta, PP_females = pp$beta, eGFRcrea_females = egfr$beta, row.names = common_snps$SNP))
        ses <- as.matrix(data.frame(CAD_females = cad$se, DBP_females = dbp$se, hematuria_females = hematuria$se, PP_females = pp$se, eGFRcrea_females = egfr$se, row.names = common_snps$SNP))
        traits=c("CAD_females", "DBP_females", "hematuria_females", "PP_females", "eGFRcrea_females")
        traits_loop=c("CAD_females", "DBP_females", "PP_females", "eGFRcrea_females")

    } else if (sex == "males") {
        betas <- as.matrix(data.frame(CAD_males = cad$beta, DBP_males = dbp$beta, hematuria_males = hematuria$beta, PP_males = pp$beta, eGFRcrea_males = egfr$beta, row.names = common_snps$SNP))
        ses <- as.matrix(data.frame(CAD_males = cad$se, DBP_males = dbp$se, hematuria_males = hematuria$se, PP_males = pp$se, eGFRcrea_males = egfr$se, row.names = common_snps$SNP))
        traits=c("CAD_males", "DBP_males", "hematuria_males", "PP_males", "eGFRcrea_males")
        traits_loop=c("CAD_males", "DBP_males", "PP_males", "eGFRcrea_males")
    }   
    rsid=common_snps$SNP

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run hyprcoloc
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # prior.1 = the ratio of the prior probability of a hypothesis in which at least one trait has a causal variant in the region by the null hypothesis that none have a causal variant
    # prior.c = the probability of a snp being associated with a second trait given its association with one trait
        # if 0.05, then there is 1 in 20 chances; if 0.01, then there is 1 in 100 chances, etc. 

    # use argument "trait.subset = c("traitX", "traitY")" to test only a subset of traits.

    prior.options <- c(1e-4, 1e-5, 1e-6)

    for (i in prior.options) {
        hyprcoloc_res <- hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid, 
                                binary.outcomes = c(1,0,1,0,0), uniform.priors = TRUE,
                                prior.1 = i, reg.steps = 2)

        write.table(hyprcoloc_res$results, here(project_dir, "results_coloc", stringr::str_c("hyprcoloc_results_", sex, "_p1_", i, ".txt")), sep = "\t", row.names = F, quote = F)
        
    }

    priorc.options <- c(0.05, 0.02, 0.01)

    for (i in priorc.options) {
        hyprcoloc_res_strictp1 <- hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid, 
                                binary.outcomes = c(1,0,1,0,0), uniform.priors = FALSE,
                                prior.1 = 1e-04, reg.steps = 2, prior.c = i)

        write.table(hyprcoloc_res_strictp1$results, here(project_dir, "results_coloc", stringr::str_c("hyprcoloc_results_", sex, "_p1_1e-04_pc_", i, ".txt")), sep = "\t", row.names = F, quote = F)
        
    }

    trait1=stringr::str_c("hematuria_", sex)

    for (sub_trait in traits_loop) {

        if (sub_trait == stringr::str_c("CAD_", sex)) {
            sub_trait_bin = 1
        } else {
            sub_trait_bin = 0
        }

        hyprcoloc_subset <- hyprcoloc(betas, ses, trait.names = traits, snp.id = rsid, 
                                trait.subset = c(trait1, sub_trait),
                                binary.outcomes = c(1, sub_trait_bin), uniform.priors = TRUE,
                                prior.1 = 1e-04, reg.steps = 2)

        write.table(hyprcoloc_subset$results, here(project_dir, "results_coloc", stringr::str_c("hyprcoloc_results_pairwise_", sex, "_", trait1, "_", sub_trait, "_p1_1e-04.txt")), sep = "\t", row.names = F, quote = F)
        
    }

    # sensitivity plot
    # pdf(here(project_dir, "results_coloc", stringr::str_c("hyprcoloc_sensitivity_plot", sex, ".pdf")))
    # sensitivity.plot(betas, ses, trait.names = traits, snp.id = rsid, 
    #                 binary.outcomes = c(1,0,1,0), ld.matrix = ld, trait.cor = cor_matrix, sample.overlap = overlap_matrix, uniform.priors = FALSE,
    #                 reg.thresh = c(0.6, 0.7, 0.8, 0.9), align.thresh = c(0.6, 0.7, 0.8, 0.9), prior.c = c(0.02, 0.01, 0.005), equal.thresholds = FALSE)
    # dev.off()
}

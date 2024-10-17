# **Title:** Investigation of genome-wide sex-differential effects in hematuria

**Summary:** We performed sex-specific genome-wide association analyses in the white British subset of the UK Biobank, and we report replication results in Million Veteran Program and Geisinger MyCode. We assessed sex-differing genetic effects for hematuria using three computational approaches: 1) sex by genotype interaction analysis of lead variants, 2) case-only genome-wide interaction (COGWI) analysis, and 3) genome-wide sexually dimorphic effects. Finally, we explored the pleiotropic COL4A2 region across cardiovascular and kidney traits through Bayesian colocalization analysis.

----------------------------------------------------------------------------
Code used in this paper is available within the /scripts directory, and below are the scripts used on each analysis, following the manuscript subheadings.

----------------------------------------------------------------------------

## **_Sex-specific GWAS of hematuria_**

**1. saige:** format_sex_hematuria_GWAS.R

**2. ldsc:** tsv_mungesumstats.R; munge_sumstats_ldsc.sh; ldsc.sh; corr_matrix.R

**3. bolt:** plink_bfile_bolt.sh; bolt_reml.sh

**4. regenie:** regenie_step1_X593.sexsp.sh; regenie_step1_X593_cond.sh; regenie_step2_RAP.sh; regenie_step2_RAP_cond.sh

**5. Rmisc:** analyse_regenie_conditional_results.R; get_conditional_snps_regenie_TOPMed.R; comparison_HRC_TOPMed.R

**6.** HLA association

**7. plots:** partitioned_h2.R; miami_hematuria_bysex.R

## **_Hematuria prevalence by sex_**

**1. Rmisc:** female_male_age_deciles_prev.R; age_stratification_X593_phenotype.R; analyse_menopause_status_associations.R

**2. regenie:** regenie_step1_X593_groups.sh; regenie_step2_age_groups_RAP.sh

## _**Genome-wide gene burden analysis of LoF variants in UKB WES data**_

**1. regenie:** regenie_step2_sexstratified_genome_wide_burden.sh; regenie_step2_sexstratified_genome_wide_burden_lovo.sh

**2. plots:** miami_plot_burden.R

## _**Assessment of sex-differential effects in hematuria**_

**1. Rmisc:** power_analysis_interaction.R; SDE.R; sde_analysis.R

**2. plots:** manhattan_plot_SDE.R; manhaattan_plot_sex_X593caseonly.R

**3. regenie:** regenie_step2_interaction.sh; regenie_step1_sex_C593_caseonly.sh; regenie_step1_check_assumption_COGWI.sh; regenie_step2_sex_X593_caseonly.sh; regenie_step2_check_assumption_COGWI.sh; regenie_step2_sex_X593_caseonly_chrX.sh

## _**Colocalization with cardiovascular and kidney related traits in the COL4A2 region**_

**1. regenie:** regenie_step1_egfr.sh; regenie_step2_egfr.sh; regenie_step2_egfr_chrX.sh

**2. hyprcoloc:** format_datasets.R; format_datasets_BPtraits.R; format_datasets_egfr_coloc.R; hyprcoloc.R

**3. Rmisc:** egfr_phenotype.R; post_regenie_egfr.R; regional_plots_chr13_UKB_LD.R

**4. plink:** plink_LD_UKB_TOPMed_RAP.sh; plink_regional_ld.sh

**5. plots:** miami_eGFR_bysex.R

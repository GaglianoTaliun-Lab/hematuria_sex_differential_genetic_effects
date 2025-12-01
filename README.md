# **Title:** Investigation of genome-wide sex-differential effects in hematuria

**Summary:** 
Sex differences play a significant role in kidney disease. We systematically assessed potential sex differential genetic effects for hematuria in a large population-based biobank followed by replication in independent cohorts (Million Veteran Program and Geisinger MyCode).  We identified significant sex interaction effects of genetic variants in two autosomal type IV collagen genes (COL4A3 and COL4A2). Our results indicate that the association between autosomal genetic variation and hematuria can differ by sex. Finally, we explored the pleiotropic COL4A2 region across cardiovascular and kidney traits through Bayesian colocalization analysis. Ultimately, these findings should be considered in the development of personalized therapies tailored by sex. 

----------------------------------------------------------------------------

Sex-stratified GWAS summary statistics of hematuria (HRC imputed) are available here: [Male-specific](https://my.locuszoom.org/gwas/494413/?token=39c44f7deef9427db2624e74d5d869fd) and [Female-specific](https://my.locuszoom.org/gwas/496796/?token=c508ca8a88d4412ba282cc466cab5dac)

Code used in this paper is available in this repository within the /scripts directory, and below are the scripts used on each analysis, following the manuscript subheadings.

----------------------------------------------------------------------------

## **_Sex-specific GWAS of hematuria_**

**1. saige:** format_sex_hematuria_GWAS.R ____

**2. ldsc:** tsv_mungesumstats.R; munge_sumstats_ldsc.sh; ldsc.sh; corr_matrix.R

**3. bolt:** plink_bfile_bolt.sh; bolt_reml.sh

**4. regenie:** regenie_step1_X593.sexsp.sh; regenie_step1_X593_cond.sh; regenie_step2_RAP.sh; regenie_step2_RAP_cond.sh

**5. Rmisc:** analyse_regenie_conditional_results.R; get_conditional_snps_regenie_TOPMed.R; comparison_HRC_TOPMed.R

**6.** HLA association _____

**7. plots:** partitioned_h2.R; miami_hematuria_bysex.R; locuszoomr_plots_chr13.R; locuszoomr_plots_chr19.R

## **_Hematuria prevalence by sex_**

**1. Rmisc:** female_male_age_deciles_prev.R; age_stratification_X593_phenotype.R; menopause_stratification_X593.R; analyse_menopause_status_associations.R; hematuria_ICDcodes_timestamps.R; hematuria_ICScodes.ipynb

**2. regenie:** regenie_step1_X593_groups.sh; regenie_step2_age_groups_RAP.sh

## _**Assessment of sex-differential effects in hematuria**_

**1. Rmisc:** power_analysis_interaction.R; SDE.R; sde_analysis.R

**2. plots:** manhattan_plot_SDE.R; manhaattan_plot_sex_X593caseonly.R

**3. regenie:** regenie_step2_interaction_RAP.sh; regenie_step1_sex_C593_caseonly.sh; regenie_step1_check_assumption_COGWI.sh; regenie_step2_sex_X593_caseonly.sh; regenie_step2_check_assumption_COGWI.sh; regenie_step2_sex_X593_caseonly_chrX.sh; regenie_step2_X593_WGS_RAP.sh

## _**Assocation of sex hormone levels with hematuria in the UKB**_

**1. Rmisc:** estradiol_table_updated.R; estradiol_analysis_for_hematuria.R; testosterone_associations.R

## _**Colocalization with cardiovascular and kidney related traits in the COL4A2 region**_

**1. regenie:** regenie_step1_egfr.sh; regenie_step2_egfr.sh; regenie_step2_egfr_chrX.sh

**2. hyprcoloc:** format_datasets.R; format_datasets_BPtraits.R; format_datasets.R; hyprcoloc.R; sharepro.sh

**3. Rmisc:** egfr_phenotype.R; post_regenie_egfr.R

**4. plink:** plink_LD_UKB_TOPMed_RAP.sh; plink_regional_ld.sh

**5. plots:** miami_eGFR_bysex.R; locuszoomr_plots_coloc_traits_females.R; locuszoomr_plots_coloc_traits_males.R

## _**Genome-wide gene burden analysis of LoF variants in UKB WES data**_

**1. regenie:** regenie_step2_sexstratified_genome_wide_burden.sh; regenie_step2_sexstratified_genome_wide_burden_lovo.sh

**2. plots:** miami_plot_burden.R

## _**Sex-specific TWAS of hematuria**_

**1. TWAS:** otters_stage1.sh; otters_stage2.sh; otters_acat.R; plink_otters_input.sh; coloc_eQTLs.R
   
**2. plots:** otters_miami_plots.R

## _**Sex-specific PWAS of hematuria**_

**1. PWAS:** bliss.sh
   
**2. plots:** miami_bliss_bysex.R

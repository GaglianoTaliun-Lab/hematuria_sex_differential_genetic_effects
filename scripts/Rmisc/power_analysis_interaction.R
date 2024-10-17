library(dplyr)
library(here)
library(tidyr)
library(stringr)
library(data.table)
library(powerGWASinteraction)
# https://rdrr.io/cran/powerGWASinteraction/man/powerGE.html

project_dir <- "/Users/frida/Documents/research-projects/col4a2_hematuria"

#  examples: -----------------------------------------------------------------------------
mod1 <- list(prev=0.01,pGene=0.2,pEnv=0.2,beta.LOR=log(c(1.0,1.2,1.4)),orGE=1.2,nSNP=10^6)
results <- powerGE(n=20000, model=mod1,alpha1=.01)
print(results)

mod2 <- list(prev=0.01,pGene=0.2,pEnv=0.2,beta.LOR=log(c(1.0,1.0,1.4)),orGE=1,nSNP=10^6)
results <- powerGE(power=0.8, model=mod2,alpha1=.01)
print(results)
# ----------------------------------------------------------------------------------------

# sample size
n=407472

females=7061+213033
males=9774+177666
all=females+males
prop_females=(females)/all

# pGene
raw_gen <- read.table(here(project_dir, "interaction", "raw_genotypes_chr13.raw"), sep = "\t", header = T) %>%
  mutate(rs7323228_HET_round = round(rs7323228_HET))
pGene_a=dplyr::count(raw_gen, rs7323228_HET_round) %>% .[1,2]/nrow(raw_gen)
pGene_b=dplyr::count(raw_gen, rs7323228_HET_round) %>% .[2,2]/nrow(raw_gen)

orGE=exp(-0.524157)

# OR of sex in hematuria
hem <- read.table("/Users/frida/Documents/research-projects/hematuria_project/UKB-WES-burden/files_for_RAP_regenie/case.control-wb.withcovariates-4regenie.txt", sep = "\t", header = T)
glm_model <- glm(X593 ~ f.22001.0.0,family=binomial(link='logit'),data=hem)
summary(glm_model)

# pGene = probability that a binary SNP is 1
# pEnv = frequency of the binary environmental variable
# orGE = odds ratio between the binary SNP and binary environmental variable (beta of ADD-INT_VAR in regenie)
# beta.LOR = vector (length=3) with log(odds ratio) of 1) genetic, 2) environmental (beta of regression hematuria ~ sex) and 3) interaction effect.
model1=list(prev=0.03, pGene= pGene_a, pEnv= prop_females, orGE= orGE, beta.LOR= c(-0.1, 0.51, 0.07), nSNP=10^6)
model2=list(prev=0.03, pGene= pGene_a, pEnv= prop_females, orGE= orGE, beta.LOR=log(c(1.0,1.2,1.4)), nSNP=10^6)


# ratio of cases/controls
caco=c(16833/390639, 60000/390639, 120000/390639, 240000/390639)

alpha=0.05

# alpha if screening was performed before first stage (if not, then alpha1=1)
alpha1=1

results = array()

for (i in 1:length(caco)) {
  results[i] <- powerGE(n = n, model = model1, caco = caco[i], alpha = alpha, alpha1 = 1) %>% .$power %>% .[1,1]
}


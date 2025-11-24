

###################################################################################################
# Load Packages -----------------------------
###################################################################################################

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)
# library(VGAM)
library(AER)

###################################################################################################
# Arguments ---------------------------------
###################################################################################################

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

###################################################################################################
# Read table ---------------------------------
###################################################################################################

estradiol <- read.table(here(project_dir, "estradiol", "final_estradiol_table_clean_updated.txt"), sep = "\t", header = T)

estradiol_tobit <- estradiol %>%
  mutate(f.30800.0.0_tobit = case_when(
    f.30806.bin == 1 ~ f.30800.0.0,
    f.30806.bin == 0 ~ 174,
    is.na(f.30806.bin) ~ f.30806.bin
  ))

# stratify by menopause status:
B_menopause_estradiol <- filter(estradiol, menopause_status_estradiol == 0)
B_menopause_X593 <- filter(estradiol, menopause_status_X593 == 0)

A_menopause_estradiol <- filter(estradiol, menopause_status_estradiol == 1)
A_menopause_X593 <- filter(estradiol, menopause_status_X593 == 1)

# for tobit:
B_menopause_estradiol_t <- filter(estradiol_tobit, menopause_status_estradiol == 0)
B_menopause_X593_t <- filter(estradiol_tobit, menopause_status_X593 == 0)

A_menopause_estradiol_t <- filter(estradiol_tobit, menopause_status_estradiol == 1)
A_menopause_X593_t <- filter(estradiol_tobit, menopause_status_X593 == 1)

# create empty dataframes for outputs:

output_X593 <- data.frame(matrix(nrow = 5, ncol = 8))
colnames(output_X593) <- c("OR_quant", "lowCI_quant", "highCI_quant", "pval_quant", "OR_bin", "lowCI_quant", "highCI_quant", "pval_bin")
rownames(output_X593) <- c("all", "pre_estradiol", "post_estradiol", "pre_X593", "post_X593")

###################################################################################################
# Main ---------------------------------
###################################################################################################

# ~~~ effect of binary estradiol:
# model the effect of binary estradiol on hematuria:
model <- glm(X593 ~ f.30806.bin + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + 
               f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + 
               f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
             data = estradiol, family=binomial(link='logit'))

summary(model)

beta <- coef(model)["f.30806.bin"]
SE <- summary(model)$coefficients[,2]["f.30806.bin"]
output_X593[1,5] <- round(exp(beta), 5)
output_X593[1,6] <- round(exp(beta-1.96*SE),5)
output_X593[1,7] <- round(exp(beta+1.96*SE),5)
output_X593[1,8] <- round(coef(summary(model))[, "Pr(>|z|)"]["f.30806.bin"], 3)

# ~~~ stratification by menopause status

##### ~~~ effect of binary estradiol:

# on hematuria
### before menopause
model1 <- glm(X593 ~ f.30806.bin + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = B_menopause_estradiol, family=binomial(link='logit'))
model2 <- glm(X593 ~ f.30806.bin + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = B_menopause_X593, family=binomial(link='logit'))

summary(model1)
summary(model2)

beta <- coef(model1)["f.30806.bin"]
SE <- summary(model1)$coefficients[,2]["f.30806.bin"]
output_X593[2,5] <- round(exp(beta), 5)
output_X593[2,6] <- round(exp(beta-1.96*SE),5)
output_X593[2,7] <- round(exp(beta+1.96*SE),5)
output_X593[2,8] <- round(coef(summary(model1))[, "Pr(>|z|)"]["f.30806.bin"], 3)

beta <- coef(model2)["f.30806.bin"]
SE <- summary(model2)$coefficients[,2]["f.30806.bin"]
output_X593[4,5] <- round(exp(beta), 5)
output_X593[4,6] <- round(exp(beta-1.96*SE),5)
output_X593[4,7] <- round(exp(beta+1.96*SE),5)
output_X593[4,8] <- round(coef(summary(model2))[, "Pr(>|z|)"]["f.30806.bin"], 3)

### after menopause
model1 <- glm(X593 ~ f.30806.bin + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = A_menopause_estradiol, family=binomial(link='logit'))
model2 <- glm(X593 ~ f.30806.bin + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = A_menopause_X593, family=binomial(link='logit'))

summary(model1)
summary(model2)

beta <- coef(model1)["f.30806.bin"]
SE <- summary(model1)$coefficients[,2]["f.30806.bin"]
output_X593[3,5] <- round(exp(beta), 5)
output_X593[3,6] <- round(exp(beta-1.96*SE),5)
output_X593[3,7] <- round(exp(beta+1.96*SE),5)
output_X593[3,8] <- round(coef(summary(model1))[, "Pr(>|z|)"]["f.30806.bin"], 3)

beta <- coef(model2)["f.30806.bin"]
SE <- summary(model2)$coefficients[,2]["f.30806.bin"]
output_X593[5,5] <- round(exp(beta), 5)
output_X593[5,6] <- round(exp(beta-1.96*SE),5)
output_X593[5,7] <- round(exp(beta+1.96*SE),5)
output_X593[5,8] <- round(coef(summary(model2))[, "Pr(>|z|)"]["f.30806.bin"], 3)


# ~~~~ effect of quantitative estradiol:

# model the effect of quantitative estradiol on hematuria:
model <- glm(X593 ~ f.30800.0.0 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
             data = estradiol, family=binomial(link='logit'))
summary(model)

tobit_model <- tobit(formula = f.30800.0.0_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                     left = 174, right = Inf,
                    data = estradiol_tobit)

summary(tobit_model)

beta <- coef(tobit_model)["X593"]
SE <- summary(tobit_model)$coefficients[,2]["X593"]
output_X593[1,1] <- round(exp(beta), 3)
output_X593[1,2] <- round(exp(beta-1.96*SE),2)
output_X593[1,3] <- round(exp(beta+1.96*SE),2)
output_X593[1,4] <- round(coef(summary(tobit_model))[, "Pr(>|z|)"]["X593"], 3)

###### ~~~ effect of quantitative estradiol:

# on hematuria
### before menopause
model1 <- glm(X593 ~ f.30800.0.0 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = B_menopause_estradiol, family=binomial(link='logit'))
model2 <- glm(X593 ~ f.30800.0.0 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = B_menopause_X593, family=binomial(link='logit'))

summary(model1)
summary(model2)

tobit_model1 <- tobit(formula = f.30800.0.0_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                     data = B_menopause_estradiol_t, left = 174, right = Inf)
tobit_model2 <- tobit(formula = f.30800.0.0_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                      data = B_menopause_X593_t, left = 174, right = Inf)

summary(tobit_model1)
summary(tobit_model2) # significant

beta <- coef(tobit_model1)["X593"]
SE <- summary(tobit_model1)$coefficients[,2]["X593"]
output_X593[2,1] <- round(exp(beta), 5)
output_X593[2,2] <- round(exp(beta-1.96*SE),5)
output_X593[2,3] <- round(exp(beta+1.96*SE),5)
output_X593[2,4] <- round(coef(summary(tobit_model1))[, "Pr(>|z|)"]["X593"], 3)

beta <- coef(tobit_model2)["X593"]
SE <- summary(tobit_model2)$coefficients[,2]["X593"]
output_X593[4,1] <- round(exp(beta), 5)
output_X593[4,2] <- round(exp(beta-1.96*SE),5)
output_X593[4,3] <- round(exp(beta+1.96*SE),5)
output_X593[4,4] <- round(coef(summary(tobit_model2))[, "Pr(>|z|)"]["X593"], 3)

### after menopause
model1 <- glm(X593 ~ f.30800.0.0 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = A_menopause_estradiol, family=binomial(link='logit'))
model2 <- glm(X593 ~ f.30800.0.0 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
              data = A_menopause_X593, family=binomial(link='logit'))

summary(model1)
summary(model2)

tobit_model1 <- tobit(formula = f.30800.0.0_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                      data = A_menopause_estradiol_t, left = 174, right = Inf)
tobit_model2 <- tobit(formula = f.30800.0.0_tobit ~ X593 + f.34.0.0 + f.22009.0.1 + f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + f.22009.0.9 + f.22009.0.10,
                      data = A_menopause_X593_t, left = 174, right = Inf)

summary(tobit_model1)
summary(tobit_model2)

beta <- coef(tobit_model1)["X593"]
SE <- summary(tobit_model1)$coefficients[,2]["X593"]
output_X593[3,1] <- round(exp(beta), 5)
output_X593[3,2] <- round(exp(beta-1.96*SE),5)
output_X593[3,3] <- round(exp(beta+1.96*SE),5)
output_X593[3,4] <- round(coef(summary(tobit_model1))[, "Pr(>|z|)"]["X593"], 3)

beta <- coef(tobit_model2)["X593"]
SE <- summary(tobit_model2)$coefficients[,2]["X593"]
output_X593[5,1] <- round(exp(beta), 5)
output_X593[5,2] <- round(exp(beta-1.96*SE),5)
output_X593[5,3] <- round(exp(beta+1.96*SE),5)
output_X593[5,4] <- round(coef(summary(tobit_model2))[, "Pr(>|z|)"]["X593"], 3)


write.table(output_X593, here(project_dir, "estradiol", "observational_results_hematuria_tobit_FLD.txt"), row.names = T, quote = F, sep = "\t")


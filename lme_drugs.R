# Ger Loughnane 10/11/2017
# Quick script for running linear mixed effects models on Drugs data

library(lme4)

drug_tab <- read.csv("C:/Users/loughnge/Dropbox/Drugs/Drug Analysis/Analysis_12_12_16/drug_lme.csv", 
                     header = FALSE, stringsAsFactors = FALSE)

# Run LME model
summary(lme_mod <-lmer(V3 ~ 1 + V4 + V2 + V4:V2 +
                         (V4 + V2 + V4:V2 | V1) + 
                         (1|V2), data = drug_tab, REML=FALSE, na.action = na.omit))

# extract coefficients
coefs <- data.frame(coef(summary(lme_mod)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
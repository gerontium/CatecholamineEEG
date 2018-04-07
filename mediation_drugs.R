# Ger Loughnane 10/11/2017
# Quick script for running linear mixed effects mediation on Drugs data

library(mediation)
library(lme4)

drug_tab <- read.csv("C:/Users/loughnge/Dropbox/Research/Drugs/Drug Analysis/lme_checking/drug_lme_p2.csv", 
                     header = FALSE, stringsAsFactors = FALSE)
colnames(drug_tab) <- c('Subject','Drug','Amp','RT')


# fit total effect, X -> Y
print('Total Effect')
total.fit <-lmer(RT ~ Drug + (1|Subject), 
               data = drug_tab, REML=FALSE, na.action = na.omit)
anova(total.fit)
# extract coefficients
coefs <- data.frame(coef(summary(total.fit)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


# fit mediation effect, X -> M
print('A path')
med.fit <-lmer(Amp ~ Drug + (1|Subject), 
                data = drug_tab, REML=FALSE, na.action = na.omit)
anova(med.fit)
# extract coefficients
coefs <- data.frame(coef(summary(med.fit)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


# fit outcome effect, M -> Y
print('B path')
out.fit <-lmer(RT ~ Drug + Amp + (1|Subject), 
                data = drug_tab, REML=FALSE, na.action = na.omit)
anova(out.fit)
# extract coefficients
coefs <- data.frame(coef(summary(out.fit)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

# These fitted objects can then be fed into the mediate function in the usual manner.
print('Mediation')
medmod <- mediate(med.fit, out.fit, treat = "Drug", mediator = "Amp")

# Show results of the mixed effects mediation model
summary(medmod)

# Ger Loughnane 01/11/2017
# Quick script for running repeated measures ANOVA on Drugs data, covarying for VAS measures

library(dplyr); library(plyr);
packages <- c("readxl", "tidyr", "ez", "ggplot2", "lme4")
sapply(packages, require, character.only = TRUE, quietly = TRUE)
setwd("C:/Users/loughnge/Dropbox/Research/Drugs/Drug Analysis/Analysis_12_12_16/stats_excels")

datatemp <- data.frame(read_excel('latency_slope_data_covar.xls',col_names = T),
                       stringsAsFactors = FALSE)

# change factors to factor class
datatemp$ID <- as.factor(datatemp$ID)
for (i in 2:(ncol(datatemp)-2)) {
  datatemp[i] <- lapply(datatemp[i], function(x) as.factor(x))
}

# plot
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

source("summarySE.R") 
source("summarySEwithin.R") # function to calculate within subjects Std.Er of mean 
source("normDataWithin.R")
plotdata <- summarySEwithin(datatemp, measurevar="DV", 
                            withinvars=c("Factor_1"), idvar="ID")

# error plot with indivdual data points
ggplot(datatemp, aes(x=Factor_1, y=DV)) + 
  geom_point(size=1.5, colour=cbPalette[3]) +
  geom_line(data=datatemp, aes(x=Factor_1, y=DV, group=ID),
            size=1, alpha=0.3, colour=cbPalette[3]) + 
  geom_line(data=plotdata, aes(x=Factor_1, y=DV, group=1), 
            size=1.5, colour=cbPalette[6]) +
  geom_errorbar(data=plotdata, aes(ymin=DV-ci, ymax=DV+ci), 
                colour='black', size=1, width=0.2) +
  geom_point(data=plotdata, aes(x=Factor_1, y=DV, group=1), 
             size=4, shape=21, fill=cbPalette[6]) +
  
  coord_fixed(ratio = 0.015) +
  theme_light() + 
  scale_x_discrete(name="Drug", labels=c('MPH', 'ATM', 'CIT', 'PLA')) +
  ylab("RT (ms)") +
  ggtitle("RT x Drug")


# run ANOVA
ezANOVA(
  data = datatemp
  , dv = DV
  , wid = .(ID)
  , within = .(Factor_1)
  , observed = .(Factor_1)
  , type = 1
  , detailed = TRUE
  , return_aov = TRUE
)

# pairwise.t.test(x=datatemp$DV, g=datatemp$Factor_1, 
#                 p.adjust="holm", paired = TRUE)

ptemp <- pairwise.t.test(x=datatemp$DV, g=datatemp$Factor_1, 
                         p.adjust="none", paired = TRUE)

p.adjust(ptemp$p.value[3, ], method = "holm")

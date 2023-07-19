###################################################################################
#
# PNR ground beetle data
#
# Changes in ground beetle abundance, richness, and diversity among
# disturbance treatments
# Years 2013, 2014, 2015
#
# Activity-abundance
# Species richness and species diversity
# Community-weights means for each trait and functional diversity
#
# KI Perry; 7 July 2023
#
###################################################################################

## Load necessary packages

library(lme4)
library(lmerTest)
library(blmeco)
library(car)
library(multcomp)
library(pscl)
library(vegan)
library(FD)
library(emmeans)
library(nortest)

library(tidyverse)
library(ggpubr)
library(rstatix)

################################################################################
## Start with 2013
## Species data
## Generate new variables: total abundance, species richness, and Shannon diversity

a13 <- read.csv("./PNR_Carabidae_2013.csv")
str(a13)

# change treatments and sampling interval to factors
a13$Canopy <- as.factor(a13$Canopy)
a13$Veg <- as.factor(a13$Veg)
a13$Interval <- as.factor(a13$Interval)
a13$Treatment <- as.factor(a13$Treatment)
a13$Quadrat <- as.factor(a13$Quadrat)
str(a13)

# check data
rowSums(a13[10:57])
colSums(a13[10:57])
head(a13)
tail(a13)
dim(a13)

summary(a13)

# need to remove extra control sites
a13.1 <- a13[which(a13$Disturbance == "Pre-Disturbance"),]

# need to remove columns with zero values
# removes species not collected in these sites during this year
a13.2 <- a13.1[, colSums(a13.1 !=0) > 0]

# double check
colSums(a13.2[10:50]) # good to go

# create several new variables and add to data set

# total abundance
a13.2$tabund <- rowSums(a13.2[10:50])
str(a13.2)
dotchart(a13.2$tabund, group = a13.2$Treatment, pch = 19)
hist(a13.2$tabund)

# species richness
specnumber(a13.2[10:50])
a13.2$rich <- specnumber(a13.2[10:50])
str(a13.2)
dotchart(a13.2$rich, group = a13.2$Treatment, pch = 19)

# species (alpha) diversity
diversity(a13.2[10:50], index = "shannon")
a13.2$div <- diversity(a13.2[10:50], index = "shannon")
str(a13.2)
dotchart(a13.2$div, group = a13.2$Treatment, pch = 19)

################################################################################
## Now 2014 species data

a14 <- read.csv("./PNR_Carabidae_2014.csv")
str(a14)

# change treatments and sampling interval to factors
a14$Canopy <- as.factor(a14$Canopy)
a14$Veg <- as.factor(a14$Veg)
a14$Interval <- as.factor(a14$Interval)
a14$Treatment <- as.factor(a14$Treatment)
a14$Quadrat <- as.factor(a14$Quadrat)
str(a14)

# check data
rowSums(a14[10:62], na.rm = TRUE)
colSums(a14[10:62], na.rm = TRUE)
head(a14)
tail(a14)
dim(a14)

# need to remove the first sampling interval with all the NAs
# couldn't sample in sites when treatments were being established
a14.1 <- a14[37:252,]
summary(a14.1)

# need to remove extra control sites
a14.2 <- a14.1[which(a14.1$Disturbance == "Post-Disturbance"),]

# need to remove columns with zero values
# removes species not collected in these sites during this year
a14.3 <- a14.2[, colSums(a14.2 !=0) > 0]

# double check
colSums(a14.3[10:49]) # good to go

# create several new variables and add to data set

# total abundance
a14.3$tabund <- rowSums(a14.3[10:49])
str(a14.3)
dotchart(a14.3$tabund, group = a14.3$Treatment, pch = 19)

# species richness
specnumber(a14.3[10:49])
a14.3$rich <- specnumber(a14.3[10:49])
str(a14.3)
dotchart(a14.3$rich, group = a14.3$Treatment, pch = 19)

# species diversity
diversity(a14.3[10:49], index = "shannon")
a14.3$div <- diversity(a14.3[10:49], index = "shannon")
str(a14.3)
dotchart(a14.3$div, group = a14.3$Treatment, pch = 19)

################################################################################
## Now 2015 species data

a15 <- read.csv("./PNR_Carabidae_2015.csv")
str(a15)

# change treatments and sampling interval to factors
a15$Canopy <- as.factor(a15$Canopy)
a15$Veg <- as.factor(a15$Veg)
a15$Interval <- as.factor(a15$Interval)
a15$Treatment <- as.factor(a15$Treatment)
a15$Quadrat <- as.factor(a15$Quadrat)

str(a15)

# check data
rowSums(a15[10:66])
colSums(a15[10:66])
head(a15)
tail(a15)
dim(a15)

summary(a15)

# need to remove extra control sites
a15.1 <- a15[which(a15$Disturbance == "Post-Disturbance"),]

# need to remove columns with zero values
# removes species not collected in these sites during this year
a15.2 <- a15.1[, colSums(a15.1 !=0) > 0]

# double check
colSums(a15.2[10:51]) # good to go

# create several new variables and add to data set

# total abundance
a15.2$tabund <- rowSums(a15.2[10:51])
str(a15.2)
dotchart(a15.2$tabund, group = a15.2$Treatment, pch = 19)

# species richness
specnumber(a15.2[10:51])
a15.2$rich <- specnumber(a15.2[10:51])
str(a15.2)
dotchart(a15.2$rich, group = a15.2$Treatment, pch = 19)

# species diversity
diversity(a15.2[10:51], index = "shannon")
a15.2$div <- diversity(a15.2[10:51], index = "shannon")
str(a15.2)
dotchart(a15.2$div, group = a15.2$Treatment, pch = 19)

################################################################################







################################################################################
## GLMMs for community-weighted means

# Cannot do repeated measures because CWM and functional alpha diversity
# calculations do not work if a site has zero beetles collected.
# Did not want to remove sites because then would have unbalanced design.
# Decided to pool abundances across sampling intervals for these analyses

# Start with 2013
cwm.13 <- read.csv("./CWMs_2013.csv", row.names = 1)
str(cwm.13)

# change predictors to factors
cwm.13$Treatment <- as.factor(cwm.13$Treatment)
cwm.13$Canopy <- as.factor(cwm.13$Canopy)
cwm.13$Veg <- as.factor(cwm.13$Veg)
str(cwm.13)

# GLMs for all response variables

# Dispersal capacity: Brachypterous
hist(cwm.13$dis_0)
plot(cwm.13$dis_0, pch = 19)
boxplot(cwm.13$dis_0 ~ cwm.13$Canopy)
boxplot(cwm.13$dis_0 ~ cwm.13$Veg)
boxplot(cwm.13$dis_0 ~ cwm.13$Treatment)
dotchart(cwm.13$dis_0, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(dis_0 ~ Canopy))
with(cwm.13, bartlett.test(dis_0 ~ Veg))
with(cwm.13, ad.test(dis_0))

dis0_mod <- glm(dis_0 ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(dis0_mod)
plot(dis0_mod)
qqnorm(resid(dis0_mod))
qqline(resid(dis0_mod))
densityPlot(rstudent(dis0_mod))
Anova(dis0_mod)
emmeans(dis0_mod, pairwise ~ Canopy)
emmeans(dis0_mod, pairwise ~ Veg)

# Dispersal capacity: Di-morphic
hist(cwm.13$dis_1) # not enough data to run model

# Dispersal capacity: Macropterous
hist(cwm.13$dis_2)
plot(cwm.13$dis_2, pch = 19)
boxplot(cwm.13$dis_2 ~ cwm.13$Canopy)
boxplot(cwm.13$dis_2 ~ cwm.13$Veg)
boxplot(cwm.13$dis_2 ~ cwm.13$Treatment)
dotchart(cwm.13$dis_2, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(dis_2 ~ Canopy))
with(cwm.13, bartlett.test(dis_2 ~ Veg))
with(cwm.13, ad.test(dis_2))

dis2_mod <- glm(dis_2 ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(dis2_mod)
plot(dis2_mod)
qqnorm(resid(dis2_mod))
qqline(resid(dis2_mod))
densityPlot(rstudent(dis2_mod))
Anova(dis2_mod)
emmeans(dis2_mod, pairwise ~ Canopy)
emmeans(dis2_mod, pairwise ~ Veg)

# Body length
hist(cwm.13$bl)
plot(cwm.13$bl, pch = 19)
boxplot(cwm.13$bl ~ cwm.13$Canopy)
boxplot(cwm.13$bl ~ cwm.13$Veg)
boxplot(cwm.13$bl ~ cwm.13$Treatment)
dotchart(cwm.13$bl, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(bl ~ Canopy))
with(cwm.13, bartlett.test(bl ~ Veg))
with(cwm.13, ad.test(bl))

bl_mod <- glm(bl ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(bl_mod)
plot(bl_mod)
qqnorm(resid(bl_mod))
qqline(resid(bl_mod))
densityPlot(rstudent(bl_mod))
Anova(bl_mod)
emmeans(bl_mod, pairwise ~ Canopy)
emmeans(bl_mod, pairwise ~ Veg)

# Body length range
hist(cwm.13$blr)
plot(cwm.13$blr, pch = 19)
boxplot(cwm.13$blr ~ cwm.13$Canopy)
boxplot(cwm.13$blr ~ cwm.13$Veg)
boxplot(cwm.13$blr ~ cwm.13$Treatment)
dotchart(cwm.13$blr, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(blr ~ Canopy))
with(cwm.13, bartlett.test(blr ~ Veg))
with(cwm.13, ad.test(blr))

blr_mod <- glm(blr ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(blr_mod)
plot(blr_mod)
qqnorm(resid(blr_mod))
qqline(resid(blr_mod))
densityPlot(rstudent(blr_mod))
Anova(blr_mod)
emmeans(blr_mod, pairwise ~ Canopy)
emmeans(blr_mod, pairwise ~ Veg)

# Head width
hist(cwm.13$rhw)
plot(cwm.13$rhw, pch = 19)
boxplot(cwm.13$rhw ~ cwm.13$Canopy)
boxplot(cwm.13$rhw ~ cwm.13$Veg)
boxplot(cwm.13$rhw ~ cwm.13$Treatment)
dotchart(cwm.13$rhw, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(rhw ~ Canopy))
with(cwm.13, bartlett.test(rhw ~ Veg))
with(cwm.13, ad.test(rhw))

rhw_mod <- glm(rhw ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(rhw_mod)
plot(rhw_mod)
qqnorm(resid(rhw_mod))
qqline(resid(rhw_mod))
densityPlot(rstudent(rhw_mod))
Anova(rhw_mod)
emmeans(rhw_mod, pairwise ~ Canopy)
emmeans(rhw_mod, pairwise ~ Veg)

# Mandible length
hist(cwm.13$rml)
plot(cwm.13$rml, pch = 19)
boxplot(cwm.13$rml ~ cwm.13$Canopy)
boxplot(cwm.13$rml ~ cwm.13$Veg)
boxplot(cwm.13$rml ~ cwm.13$Treatment)
dotchart(cwm.13$rml, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(rml ~ Canopy))
with(cwm.13, bartlett.test(rml ~ Veg))
with(cwm.13, ad.test(rml))

rml_mod <- glm(rml ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(rml_mod)
plot(rml_mod)
qqnorm(resid(rml_mod))
qqline(resid(rml_mod))
densityPlot(rstudent(rml_mod))
Anova(rml_mod)
emmeans(rml_mod, pairwise ~ Canopy)
emmeans(rml_mod, pairwise ~ Veg)

# Eye width
hist(cwm.13$rew)
plot(cwm.13$rew, pch = 19)
boxplot(cwm.13$rew ~ cwm.13$Canopy)
boxplot(cwm.13$rew ~ cwm.13$Veg)
boxplot(cwm.13$rew ~ cwm.13$Treatment)
dotchart(cwm.13$rew, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(rew ~ Canopy))
with(cwm.13, bartlett.test(rew ~ Veg))
with(cwm.13, ad.test(rew))

rew_mod <- glm(rew ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(rew_mod)
plot(rew_mod)
qqnorm(resid(rew_mod))
qqline(resid(rew_mod))
densityPlot(rstudent(rew_mod))
Anova(rew_mod)
emmeans(rew_mod, pairwise ~ Canopy)
emmeans(rew_mod, pairwise ~ Veg)

outlierTest(bl.mod.null)
influenceIndexPlot(bl.mod.null, vars = c("Cook"), id = list(n = 3))
# significant outlier, let's remove it and see if the model fits better
lec_1.mod.red2 <- update(lec_1.mod.red, subset = -c(1))
summary(lec_1.mod.red2)
compareCoefs(lec_1.mod.red, lec_1.mod.red2) # compares estimated coefficients and their standard errors

# Leg length
hist(cwm.13$bll)
plot(cwm.13$bll, pch = 19)
boxplot(cwm.13$bll ~ cwm.13$Canopy)
boxplot(cwm.13$bll ~ cwm.13$Veg)
boxplot(cwm.13$bll ~ cwm.13$Treatment)
dotchart(cwm.13$bll, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(bll ~ Canopy))
with(cwm.13, bartlett.test(bll ~ Veg))
with(cwm.13, ad.test(bll))

bll_mod <- glm(bll ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(bll_mod)
plot(bll_mod)
qqnorm(resid(bll_mod))
qqline(resid(bll_mod))
densityPlot(rstudent(bll_mod))
Anova(bll_mod)
emmeans(bll_mod, pairwise ~ Canopy)
emmeans(bll_mod, pairwise ~ Veg)

# Functional alpha diversity
hist(cwm.13$rao13)
plot(cwm.13$rao13, pch = 19)
boxplot(cwm.13$rao13 ~ cwm.13$Canopy)
boxplot(cwm.13$rao13 ~ cwm.13$Veg)
boxplot(cwm.13$rao13 ~ cwm.13$Treatment)
dotchart(cwm.13$rao13, group = cwm.13$Treatment, pch = 19)
with(cwm.13, bartlett.test(rao13 ~ Canopy))
with(cwm.13, bartlett.test(rao13 ~ Veg))
with(cwm.13, ad.test(rao13))

rao13_mod <- glm(rao13 ~ Canopy * Veg, family = gaussian, data = cwm.13)
summary(rao13_mod)
plot(rao13_mod)
qqnorm(resid(rao13_mod))
qqline(resid(rao13_mod))
densityPlot(rstudent(rao13_mod))
Anova(rao13_mod)
emmeans(rao13_mod, pairwise ~ Canopy)
emmeans(rao13_mod, pairwise ~ Veg)


# 2014
cwm.14 <- read.csv("./CWMs_2014.csv", row.names = 1)
str(cwm.14)

# change predictors to factors
cwm.14$Treatment <- as.factor(cwm.14$Treatment)
cwm.14$Canopy <- as.factor(cwm.14$Canopy)
cwm.14$Veg <- as.factor(cwm.14$Veg)
str(cwm.14)

# GLMs for all response variables

# Dispersal capacity: Brachypterous
hist(cwm.14$dis_0)
plot(cwm.14$dis_0, pch = 19)
boxplot(cwm.14$dis_0 ~ cwm.14$Canopy)
boxplot(cwm.14$dis_0 ~ cwm.14$Veg)
boxplot(cwm.14$dis_0 ~ cwm.14$Treatment)
dotchart(cwm.14$dis_0, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(dis_0 ~ Canopy))
with(cwm.14, bartlett.test(dis_0 ~ Veg))
with(cwm.14, ad.test(dis_0))

dis0_mod <- glm(dis_0 ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(dis0_mod)
plot(dis0_mod)
qqnorm(resid(dis0_mod))
qqline(resid(dis0_mod))
densityPlot(rstudent(dis0_mod))
Anova(dis0_mod)
emmeans(dis0_mod, pairwise ~ Canopy)
emmeans(dis0_mod, pairwise ~ Veg)

# Dispersal capacity: Di-morphic
hist(cwm.14$dis_1) # not enough data to run model

# Dispersal capacity: Macropterous
hist(cwm.14$dis_2)
plot(cwm.14$dis_2, pch = 19)
boxplot(cwm.14$dis_2 ~ cwm.14$Canopy)
boxplot(cwm.14$dis_2 ~ cwm.14$Veg)
boxplot(cwm.14$dis_2 ~ cwm.14$Treatment)
dotchart(cwm.14$dis_2, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(dis_2 ~ Canopy))
with(cwm.14, bartlett.test(dis_2 ~ Veg))
with(cwm.14, ad.test(dis_2))

dis2_mod <- glm(dis_2 ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(dis2_mod)
plot(dis2_mod)
qqnorm(resid(dis2_mod))
qqline(resid(dis2_mod))
densityPlot(rstudent(dis2_mod))
Anova(dis2_mod)
emmeans(dis2_mod, pairwise ~ Canopy)
emmeans(dis2_mod, pairwise ~ Veg)

# Body length
hist(cwm.14$bl)
plot(cwm.14$bl, pch = 19)
boxplot(cwm.14$bl ~ cwm.14$Canopy)
boxplot(cwm.14$bl ~ cwm.14$Veg)
boxplot(cwm.14$bl ~ cwm.14$Treatment)
dotchart(cwm.14$bl, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(bl ~ Canopy))
with(cwm.14, bartlett.test(bl ~ Veg))
with(cwm.14, ad.test(bl))

bl_mod <- glm(bl ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(bl_mod)
plot(bl_mod)
qqnorm(resid(bl_mod))
qqline(resid(bl_mod))
densityPlot(rstudent(bl_mod))
Anova(bl_mod)
emmeans(bl_mod, pairwise ~ Canopy)
emmeans(bl_mod, pairwise ~ Veg)

# Body length range
hist(cwm.14$blr)
plot(cwm.14$blr, pch = 19)
boxplot(cwm.14$blr ~ cwm.14$Canopy)
boxplot(cwm.14$blr ~ cwm.14$Veg)
boxplot(cwm.14$blr ~ cwm.14$Treatment)
dotchart(cwm.14$blr, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(blr ~ Canopy))
with(cwm.14, bartlett.test(blr ~ Veg))
with(cwm.14, ad.test(blr))

blr_mod <- glm(blr ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(blr_mod)
plot(blr_mod)
qqnorm(resid(blr_mod))
qqline(resid(blr_mod))
densityPlot(rstudent(blr_mod))
Anova(blr_mod)
emmeans(blr_mod, pairwise ~ Canopy)
emmeans(blr_mod, pairwise ~ Veg)

# Head width
hist(cwm.14$rhw)
plot(cwm.14$rhw, pch = 19)
boxplot(cwm.14$rhw ~ cwm.14$Canopy)
boxplot(cwm.14$rhw ~ cwm.14$Veg)
boxplot(cwm.14$rhw ~ cwm.14$Treatment)
dotchart(cwm.14$rhw, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(rhw ~ Canopy))
with(cwm.14, bartlett.test(rhw ~ Veg))
with(cwm.14, ad.test(rhw))

rhw_mod <- glm(rhw ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(rhw_mod)
plot(rhw_mod)
qqnorm(resid(rhw_mod))
qqline(resid(rhw_mod))
densityPlot(rstudent(rhw_mod))
Anova(rhw_mod)
emmeans(rhw_mod, pairwise ~ Canopy)
emmeans(rhw_mod, pairwise ~ Veg)

# Mandible length
hist(cwm.14$rml)
plot(cwm.14$rml, pch = 19)
boxplot(cwm.14$rml ~ cwm.14$Canopy)
boxplot(cwm.14$rml ~ cwm.14$Veg)
boxplot(cwm.14$rml ~ cwm.14$Treatment)
dotchart(cwm.14$rml, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(rml ~ Canopy))
with(cwm.14, bartlett.test(rml ~ Veg))
with(cwm.14, ad.test(rml))

rml_mod <- glm(rml ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(rml_mod)
plot(rml_mod)
qqnorm(resid(rml_mod))
qqline(resid(rml_mod))
densityPlot(rstudent(rml_mod))
Anova(rml_mod)
emmeans(rml_mod, pairwise ~ Canopy)
emmeans(rml_mod, pairwise ~ Veg)

# Eye width
hist(cwm.14$rew)
plot(cwm.14$rew, pch = 19)
boxplot(cwm.14$rew ~ cwm.14$Canopy)
boxplot(cwm.14$rew ~ cwm.14$Veg)
boxplot(cwm.14$rew ~ cwm.14$Treatment)
dotchart(cwm.14$rew, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(rew ~ Canopy))
with(cwm.14, bartlett.test(rew ~ Veg))
with(cwm.14, ad.test(rew))

rew_mod <- glm(rew ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(rew_mod)
plot(rew_mod)
qqnorm(resid(rew_mod))
qqline(resid(rew_mod))
densityPlot(rstudent(rew_mod))
Anova(rew_mod)
emmeans(rew_mod, pairwise ~ Canopy)
emmeans(rew_mod, pairwise ~ Veg)

# Leg length
hist(cwm.14$bll)
plot(cwm.14$bll, pch = 19)
boxplot(cwm.14$bll ~ cwm.14$Canopy)
boxplot(cwm.14$bll ~ cwm.14$Veg)
boxplot(cwm.14$bll ~ cwm.14$Treatment)
dotchart(cwm.14$bll, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(bll ~ Canopy))
with(cwm.14, bartlett.test(bll ~ Veg))
with(cwm.14, ad.test(bll))

bll_mod <- glm(bll ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(bll_mod)
plot(bll_mod)
qqnorm(resid(bll_mod))
qqline(resid(bll_mod))
densityPlot(rstudent(bll_mod))
Anova(bll_mod)
emmeans(bll_mod, pairwise ~ Canopy)
emmeans(bll_mod, pairwise ~ Veg)

# Functional alpha diversity
hist(cwm.14$rao14)
plot(cwm.14$rao14, pch = 19)
boxplot(cwm.14$rao14 ~ cwm.14$Canopy)
boxplot(cwm.14$rao14 ~ cwm.14$Veg)
boxplot(cwm.14$rao14 ~ cwm.14$Treatment)
dotchart(cwm.14$rao14, group = cwm.14$Treatment, pch = 19)
with(cwm.14, bartlett.test(rao14 ~ Canopy))
with(cwm.14, bartlett.test(rao14 ~ Veg))
with(cwm.14, ad.test(rao14))

rao14_mod <- glm(rao14 ~ Canopy * Veg, family = gaussian, data = cwm.14)
summary(rao14_mod)
plot(rao14_mod)
qqnorm(resid(rao14_mod))
qqline(resid(rao14_mod))
densityPlot(rstudent(rao14_mod))
Anova(rao14_mod)
emmeans(rao14_mod, pairwise ~ Canopy)
emmeans(rao14_mod, pairwise ~ Veg)


# 2015
cwm.15 <- read.csv("./CWMs_2015.csv", row.names = 1)
str(cwm.15)

# change predictors to factors
cwm.15$Treatment <- as.factor(cwm.15$Treatment)
cwm.15$Canopy <- as.factor(cwm.15$Canopy)
cwm.15$Veg <- as.factor(cwm.15$Veg)
str(cwm.15)

# GLMs for all response variables

# Dispersal capacity: Brachypterous
hist(cwm.15$dis_0)
plot(cwm.15$dis_0, pch = 19)
boxplot(cwm.15$dis_0 ~ cwm.15$Canopy)
boxplot(cwm.15$dis_0 ~ cwm.15$Veg)
boxplot(cwm.15$dis_0 ~ cwm.15$Treatment)
dotchart(cwm.15$dis_0, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(dis_0 ~ Canopy))
with(cwm.15, bartlett.test(dis_0 ~ Veg))
with(cwm.15, ad.test(dis_0))

dis0_mod <- glm(dis_0 ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(dis0_mod)
plot(dis0_mod)
qqnorm(resid(dis0_mod))
qqline(resid(dis0_mod))
densityPlot(rstudent(dis0_mod))
Anova(dis0_mod)
emmeans(dis0_mod, pairwise ~ Canopy)
emmeans(dis0_mod, pairwise ~ Veg)

# Dispersal capacity: Di-morphic
hist(cwm.15$dis_1) # not enough data to run model

# Dispersal capacity: Macropterous
hist(cwm.15$dis_2)
plot(cwm.15$dis_2, pch = 19)
boxplot(cwm.15$dis_2 ~ cwm.15$Canopy)
boxplot(cwm.15$dis_2 ~ cwm.15$Veg)
boxplot(cwm.15$dis_2 ~ cwm.15$Treatment)
dotchart(cwm.15$dis_2, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(dis_2 ~ Canopy))
with(cwm.15, bartlett.test(dis_2 ~ Veg))
with(cwm.15, ad.test(dis_2))

dis2_mod <- glm(dis_2 ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(dis2_mod)
plot(dis2_mod)
qqnorm(resid(dis2_mod))
qqline(resid(dis2_mod))
densityPlot(rstudent(dis2_mod))
Anova(dis2_mod)
emmeans(dis2_mod, pairwise ~ Canopy)
emmeans(dis2_mod, pairwise ~ Veg)

# Body length
hist(cwm.15$bl)
plot(cwm.15$bl, pch = 19)
boxplot(cwm.15$bl ~ cwm.15$Canopy)
boxplot(cwm.15$bl ~ cwm.15$Veg)
boxplot(cwm.15$bl ~ cwm.15$Treatment)
dotchart(cwm.15$bl, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(bl ~ Canopy))
with(cwm.15, bartlett.test(bl ~ Veg))
with(cwm.15, ad.test(bl))

bl_mod <- glm(bl ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(bl_mod)
plot(bl_mod)
qqnorm(resid(bl_mod))
qqline(resid(bl_mod))
densityPlot(rstudent(bl_mod))
Anova(bl_mod)
emmeans(bl_mod, pairwise ~ Canopy)
emmeans(bl_mod, pairwise ~ Veg)

# Body length range
hist(cwm.15$blr)
plot(cwm.15$blr, pch = 19)
boxplot(cwm.15$blr ~ cwm.15$Canopy)
boxplot(cwm.15$blr ~ cwm.15$Veg)
boxplot(cwm.15$blr ~ cwm.15$Treatment)
dotchart(cwm.15$blr, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(blr ~ Canopy))
with(cwm.15, bartlett.test(blr ~ Veg))
with(cwm.15, ad.test(blr))

blr_mod <- glm(blr ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(blr_mod)
plot(blr_mod)
qqnorm(resid(blr_mod))
qqline(resid(blr_mod))
densityPlot(rstudent(blr_mod))
Anova(blr_mod)
emmeans(blr_mod, pairwise ~ Canopy)
emmeans(blr_mod, pairwise ~ Veg)

# Head width
hist(cwm.15$rhw)
plot(cwm.15$rhw, pch = 19)
boxplot(cwm.15$rhw ~ cwm.15$Canopy)
boxplot(cwm.15$rhw ~ cwm.15$Veg)
boxplot(cwm.15$rhw ~ cwm.15$Treatment)
dotchart(cwm.15$rhw, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(rhw ~ Canopy))
with(cwm.15, bartlett.test(rhw ~ Veg))
with(cwm.15, ad.test(rhw))

rhw_mod <- glm(rhw ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(rhw_mod)
plot(rhw_mod)
qqnorm(resid(rhw_mod))
qqline(resid(rhw_mod))
densityPlot(rstudent(rhw_mod))
Anova(rhw_mod)
emmeans(rhw_mod, pairwise ~ Canopy)
emmeans(rhw_mod, pairwise ~ Veg)

# Mandible length
hist(cwm.15$rml)
plot(cwm.15$rml, pch = 19)
boxplot(cwm.15$rml ~ cwm.15$Canopy)
boxplot(cwm.15$rml ~ cwm.15$Veg)
boxplot(cwm.15$rml ~ cwm.15$Treatment)
dotchart(cwm.15$rml, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(rml ~ Canopy))
with(cwm.15, bartlett.test(rml ~ Veg))
with(cwm.15, ad.test(rml))

rml_mod <- glm(rml ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(rml_mod)
plot(rml_mod)
qqnorm(resid(rml_mod))
qqline(resid(rml_mod))
densityPlot(rstudent(rml_mod))
Anova(rml_mod)
emmeans(rml_mod, pairwise ~ Canopy)
emmeans(rml_mod, pairwise ~ Veg)

# Eye width
hist(cwm.15$rew)
plot(cwm.15$rew, pch = 19)
boxplot(cwm.15$rew ~ cwm.15$Canopy)
boxplot(cwm.15$rew ~ cwm.15$Veg)
boxplot(cwm.15$rew ~ cwm.15$Treatment)
dotchart(cwm.15$rew, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(rew ~ Canopy))
with(cwm.15, bartlett.test(rew ~ Veg))
with(cwm.15, ad.test(rew))

rew_mod <- glm(rew ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(rew_mod)
plot(rew_mod)
qqnorm(resid(rew_mod))
qqline(resid(rew_mod))
densityPlot(rstudent(rew_mod))
Anova(rew_mod)
emmeans(rew_mod, pairwise ~ Canopy)
emmeans(rew_mod, pairwise ~ Veg)

# Leg length
hist(cwm.15$bll)
plot(cwm.15$bll, pch = 19)
boxplot(cwm.15$bll ~ cwm.15$Canopy)
boxplot(cwm.15$bll ~ cwm.15$Veg)
boxplot(cwm.15$bll ~ cwm.15$Treatment)
dotchart(cwm.15$bll, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(bll ~ Canopy))
with(cwm.15, bartlett.test(bll ~ Veg))
with(cwm.15, ad.test(bll))

bll_mod <- glm(bll ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(bll_mod)
plot(bll_mod)
qqnorm(resid(bll_mod))
qqline(resid(bll_mod))
densityPlot(rstudent(bll_mod))
Anova(bll_mod)
emmeans(bll_mod, pairwise ~ Canopy)
emmeans(bll_mod, pairwise ~ Veg)

# Functional alpha diversity
hist(cwm.15$rao15)
plot(cwm.15$rao15, pch = 19)
boxplot(cwm.15$rao15 ~ cwm.15$Canopy)
boxplot(cwm.15$rao15 ~ cwm.15$Veg)
boxplot(cwm.15$rao15 ~ cwm.15$Treatment)
dotchart(cwm.15$rao15, group = cwm.15$Treatment, pch = 19)
with(cwm.15, bartlett.test(rao15 ~ Canopy))
with(cwm.15, bartlett.test(rao15 ~ Veg))
with(cwm.15, ad.test(rao15))

rao15_mod <- glm(rao15 ~ Canopy * Veg, family = gaussian, data = cwm.15)
summary(rao15_mod)
plot(rao15_mod)
qqnorm(resid(rao15_mod))
qqline(resid(rao15_mod))
densityPlot(rstudent(rao15_mod))
Anova(rao15_mod)
emmeans(rao15_mod, pairwise ~ Canopy)
emmeans(rao15_mod, pairwise ~ Veg)


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

# Cannot do repeated measures because CWM calculations do not work if a site has
# zero beetles collected. Did not want to remove sites because then would have 
# unbalanced design. Decided to pool abundances across sampling intervals for 
# CWM analyses

cwm.13 <- read.csv("./CWMs_2013.csv")
str(cwm.13)

cwm.14 <- read.csv("./CWMs_2014.csv")
str(cwm.14)

cwm.15 <- read.csv("./CWMs_2015.csv")
str(cwm.15)

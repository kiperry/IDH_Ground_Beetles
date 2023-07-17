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

################################################################################
# need to remove extra control sites

a13.1 <- a13[which(a13$Disturbance == "Pre-Disturbance"),]

# need to remove columns with zero values
# removes species not collected in these sites during this year
a13.2 <- a13.1[, colSums(a13.1 !=0) > 0]

# double check
colSums(a13.2[10:50]) # good to go

################################################################################
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

# species diversity
diversity(a13.2[10:50], index = "shannon")
a13.2$div <- diversity(a13.2[10:50], index = "shannon")
str(a13.2)
dotchart(a13.2$div, group = a13.2$Treatment, pch = 19)

################################################################################


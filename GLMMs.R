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
# load the trait data
t <- read.csv("./PNR_Carabidae_Traits.csv", row.names = 1)

names(t)
str(t)
plot(t)

# remove traits that are highly correlated
t2 <- t[,-3] # minimum body length
t2 <- t2[,-3] # maximum body length
t2 <- t2[,-7] # antennae length

str(t2)

# change dispersal ability to a factor
t2$dis <- as.factor(t2$dis)

# check traits for normality
hist(t2$bl)

hist(t2$blr)
hist(log(t2$blr))
t2$blr <- log(t2$blr + 1)

hist(t2$rhw)
hist(log(t2$rhw))
t2$rhw <- log(t2$rhw + 1)

hist(t2$rml)
hist(log(t2$rml))
t2$rml <- log(t2$rml + 1)

hist(t2$rew)
hist(log(t2$rew))
t2$rew <- log(t2$rew + 1)

hist(t2$bll)
hist(log(t2$bll))

# now we have to match up the trait matrix with every abundance matrix

# start with 2013
intersect(colnames(a13.2[10:50]), rownames(t2))

setdiff(colnames(a13.2[10:50]), rownames(t2))
a13.3 <- a13.2[,-10]
a13.3 <- a13.3[,-31]
setdiff(colnames(a13.3[10:48]), rownames(t2))

setdiff(rownames(t2), colnames(a13.3[10:48])) # need to update this now that I removed columns with zeros
t.13 <- t2[-3,]
t.13 <- t.13[-c(4:7),]
t.13 <- t.13[-5,]
t.13 <- t.13[-c(8:10),]
setdiff(rownames(t.13), colnames(a13.2[10:55]))
t.13 <- t.13[-11,]
t.13 <- t.13[-c(16:18),]
t.13 <- t.13[-22,]
t.13 <- t.13[-26,]
t.13 <- t.13[-29,]
t.13 <- t.13[-32,]
t.13 <- t.13[-35,]
t.13 <- t.13[-c(36:37),]
t.13 <- t.13[-39,]
t.13 <- t.13[-40,]
setdiff(rownames(t.13), colnames(a13.2[10:55]))

intersect(colnames(a13.2[10:55]), rownames(t.13))
rownames(t.13) == colnames(a13.2[10:55]) 

# now calculate trait variables

# community-weighed means
functcomp(t.13, as.matrix(a13.2[10:55]), CWM.type = "all")
cwm.13 <- functcomp(t.13, as.matrix(a13.2[10:55]), CWM.type = "all") # maybe this will not run with species not collected
cwm.13
write.csv(cwm.13, file = "CWMs_2013.csv")


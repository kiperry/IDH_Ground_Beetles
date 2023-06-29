###################################################################################
#
# PNR ground beetle data
#
# Changes in ground beetle species composition among disturbance treatments
# Years 2013, 2014, 2015
#
# Taxonomic & functional beta-diversity
# NMDS
# PERMANOVA
#
# KI Perry; 21 February 2023
#
###################################################################################

## Load necessary packages

library(reshape2)
library(vegan)
library(ggplot2)
library(viridis)
library(betapart)

#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

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

# need to pool species data across sampling intervals for community analyses
a13.2 <- melt(a13, id = c("Interval", "Quadrat", "Canopy", "Veg",
                         "Treatment", "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a13.2)

# double check only species in the variable column
levels(a13.2$variable)
names(a13.2)[10] <- "Species"

# now create data frame for community analyses
a13.3 <- dcast(a13.2, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same
colSums(a13.3[5:52])
colSums(a13.3[5:52]) == colSums(a13[10:57])

str(a13.3)

# need to remove extra control sites
a13.4 <- a13.3[1:24,]

# remove columns with zero values
# removes species not collected in these sites
# seven species removed
colSums(a13.3[5:52])
colSums(a13.4[5:52])

a13.4 <- a13.3[, colSums(a13.3 != 0) > 0]
colSums(a13.4[6:46]) #can use this format for community composition analyses
a13.4 <- a13.4[, colSums(a13.4 != 0) > 0]
colSums(a13.4[6:46]) #can use this format for community composition analyses

############################

## Partition beta-diversity of ground beetle communities among treatments

#change data set to presence/absence for this part of the analyses
a13.5 <- a13.4
a13.5[a13.5 > 0] <- 1

str(a13.5)
rowSums(a13.5[6:46])
colSums(a13.5[6:46])

# create beta part object for analyses
a13.core <- betapart.core(a13.5[6:46])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
a13.dist <- beta.pair(a13.core, index.family = "sorensen")
str(a13.dist)

# run NMDS models for total beta-diversity component

## Beta.sor - Total beta-diversity
nmds.a13.sor <- metaMDS(a13.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.a13.sor
stressplot(nmds.a13.sor)
goodness(nmds.a13.sor)
nmds.a13.sor$stress #stress is quality of fit
plot(nmds.a13.sor)

a13.5$Treatment <- factor(a13.5$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))
a13.5$Treatment_Code <- factor(a13.5$Treatment_Code, levels = c("A", "B", "C", "D"))


ordiplot(nmds.a13.sor, type="n", xlim = c(-0.8, 0.6), ylim = c(-0.6, 0.6))

points(nmds.a13.sor, display = "sites", pch = pchvec1[a13.5$Treatment], cex = 1.5, 
       col = colvec1[a13.5$Treatment])

ordiellipse(nmds.a13.sor, groups = a13.5$Treatment, display = "sites", draw = "lines",
            lwd = 3, conf = 0.90, col = colvec1[a13.5$Treatment])

legend("topleft", legend = c("Canopy","Understory","Canopy+Understory", "Undisturbed"), 
       pch = pchvec1[a13.5$Treatment], cex = 1.2, bty = "n", col = colvec1[a13.5$Treatment])

orditorp(nmds.a13.sor, "sites") #used to double check the legend!

# colors are messed up here, not sure why... yet.

## Test for differences in ground beetle composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(a13.dist$beta.sor ~ a13.5$Treatment, permutations = 999)
adonis2(a13.dist$beta.sor ~ a13.5$Canopy, permutations = 999)
adonis2(a13.dist$beta.sor ~ a13.5$Veg, permutations = 999)
pairwise.adonis(a13.dist$beta.sor, a13.5$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a13.beta.sor <- betadisper(a13.dist$beta.sor, a13.5$Treatment, type = c("median"))
a13.beta.sor
anova(a13.beta.sor)
plot(a13.beta.sor)
boxplot(a13.beta.sor, ylab = "Distance to median")
TukeyHSD(a13.beta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(a13.dist$beta.sim ~ a13.5$Treatment, permutations = 999)
adonis2(a13.dist$beta.sim ~ a13.5$Canopy, permutations = 999)
adonis2(a13.dist$beta.sim ~ a13.5$Veg, permutations = 999)

# nestedness
adonis2(a13.dist$beta.sne ~ a13.5$Treatment, permutations = 999)
adonis2(a13.dist$beta.sne ~ a13.5$Canopy, permutations = 999)
adonis2(a13.dist$beta.sne ~ a13.5$Veg, permutations = 999)

################################################################################
## Now 2014
## Species data
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
a14.2 <- a14[37:252,]

summary(a14.2)

# need to pool species data across sampling intervals
a14.3 <- melt(a14.2, id = c("Interval", "Quadrat", "Canopy", "Veg",
                          "Treatment", "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a14.3)

# double check only species in the variable column
levels(a14.3$variable)
names(a14.3)[10] <- "Species"

# now create data frame for community analyses
a14.4 <- dcast(a14.3, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same
colSums(a14.4[5:57])
colSums(a14.4[5:57]) == colSums(a14[37:252,10:62], na.rm = TRUE)

str(a14.4)

# need to remove extra control sites
a14.5 <- a14.4[1:24,]

# remove columns with zero values
# removes species not collected in these sites
# thirteen species removed
colSums(a14.5[5:57])

a14.5 <- a14.5[, colSums(a14.5 != 0) > 0]
colSums(a14.5[5:44])

## Compare ground beetle species richness among treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

levels(a14.5$Treatment)
rare.a14.C <- a14.5[which(a14.5$Treatment == "Control"),]
rare.a14.L <- a14.5[which(a14.5$Treatment == "Light"),]
rare.a14.V <- a14.5[which(a14.5$Treatment == "Veg"),]
rare.a14.LV <- a14.5[which(a14.5$Treatment == "Light+Veg"),]

sp.a14.C <- specaccum(rare.a14.C[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.L <- specaccum(rare.a14.L[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.V <- specaccum(rare.a14.V[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.LV <- specaccum(rare.a14.LV[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(sp.a14.C, pch = 19, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 500), ylim = c(0, 35))
plot(sp.a14.L, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 3, col = "#FDE725FF")
plot(sp.a14.V, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 3, col = "#73D055FF")
plot(sp.a14.LV, add = TRUE, pch = 9, xvar = c("individuals"), lty = 3, lwd = 3, col = "#2D708EFF")
legend("bottomright", legend = c("Canopy","Understory","Canopy+Understory", "Undisturbed"), 
       lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

levels(a14.5$Treatment)

#calculates species richness for each site
specnumber(a14.5[5:44])

#calculates species richness by treatment
specnumber(a14.5[5:44], groups = a14.5$Treatment)
str(a14.5)

bark.sp.j1 <- diversitycomp(a14.5[5:44], y = a14.5, factor1 = "Canopy", index = "jack1") #this isn't working and I'm not sure why... yet.
bark.sp.j1
bark.sp.j2 <- diversitycomp(a14.5[5:44], y = a14.5, factor1 = "Treatment", index = "jack2")
bark.sp.j2

bark.j1 <- diversityresult(a14.5[5:44], y=NULL, index = "jack1")
bark.j1
bark.j2 <- diversityresult(a14.5[5:44], y=NULL, index = "jack2")
bark.j2


#########################################################

## Partition beta-diversity of ground beetle communities among treatments

#change data set to presence/absence for this part of the analyses
a14.6 <- a14.5
a14.6[a14.6 > 0] <- 1

str(a14.6)
rowSums(a14.6[5:44])
colSums(a14.6[5:44])

# create beta part object for analyses
a14.core <- betapart.core(a14.6[5:44])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
a14.dist <- beta.pair(a14.core, index.family = "sorensen")
str(a14.dist)

# run NMDS models for total beta-diversity component

## Beta.sor - Total beta-diversity
nmds.a14.sor <- metaMDS(a14.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.a14.sor
stressplot(nmds.a14.sor)
goodness(nmds.a14.sor)
nmds.a14.sor$stress #stress is quality of fit
plot(nmds.a14.sor)

ordiplot(nmds.a14.sor, type="n", xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))
points(nmds.a14.sor, display = "sites", pch = pchvec1[a14.6$Treatment], cex = 1.5, 
       col = colvec2[a14.6$Treatment])
ordiellipse(nmds.a14.sor, groups = a14.6$Treatment, display = "sites", draw = "lines", 
            col = colvec2[a14.6$Treatment], lwd = 3, conf = 0.90)
ordihull(nmds.a14.sor, groups = a14.6$Treatment, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[a14.6$Treatment], lwd = 2.5)
legend("topleft", legend = c("Canopy","Understory","Canopy+Understory", "Undisturbed"), 
       pch = pchvec1[a14.6$Treatment], cex = 1.2, bty = "n", col = colvec2[a14.6$Treatment])
orditorp(nnmds.a14.sor, "sites") #used to double check the legend!

# colors are messed up here, not sure why... yet.

## Test for differences in ground beetle composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(a14.dist$beta.sor ~ a14.6$Treatment, permutations = 999)
adonis2(a14.dist$beta.sor ~ a14.6$Canopy, permutations = 999)
adonis2(a14.dist$beta.sor ~ a14.6$Veg, permutations = 999)
pairwise.adonis(a14.dist$beta.sor, a14.6$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a14.beta.sor <- betadisper(a14.dist$beta.sor, a14.6$Treatment, type = c("median"))
a14.beta.sor
anova(a14.beta.sor)
plot(a14.beta.sor)
boxplot(a14.beta.sor, ylab = "Distance to median")
TukeyHSD(a14.beta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(a14.dist$beta.sim ~ a14.6$Treatment, permutations = 999)
adonis2(a14.dist$beta.sim ~ a14.6$Canopy, permutations = 999)
adonis2(a14.dist$beta.sim ~ a14.6$Veg, permutations = 999)

# nestedness
adonis2(a14.dist$beta.sne ~ a14.6$Treatment, permutations = 999)
adonis2(a14.dist$beta.sne ~ a14.6$Canopy, permutations = 999)
adonis2(a14.dist$beta.sne ~ a14.6$Veg, permutations = 999)


################################################################################
## Now 2015
## Species data

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

###################################################################################

# need to pool species data across sampling intervals for community analyses
a15.2 <- melt(a15, id = c("Interval", "Quadrat", "Canopy", "Veg", "Treatment",
                          "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a15.2)

# double check only species in the variable column
levels(a15.2$variable) #yup
names(a15.2)[10] <- "Species"

# now create a data frame for community analyses
a15.3 <- dcast(a15.2, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same and we didn't mess anything up
colSums(a15.3[5:61])
colSums(a15.3[5:61]) == colSums(a15[10:66])

###################################################################################
# need to remove extra control sites
a15.3.2 <- a15.3[1:24,]
colSums(a15.3.2[5:61])

# need to remove columns with zero values
# removes species not collected in these sites during this year
a15.4 <- a15.3.2[, colSums(a15.3.2 !=0) > 0]

# removed fifteen species
colSums(a15.4[5:46])
# can use this format for community composition analyses

##################################################################################

## Partition beta-diversity of ground beetle communities among treatments

# change data set to presence/absence
a15.5 <- a15.4 # make a copy
a15.5[a15.5 > 0] <- 1

str(a15.5)
rowSums(a15.5[5:46])
colSums(a15.5[5:46])

a15.5$Treatment <- factor(a15.5$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))

## Taxonomic beta-diversity

# create a beta part object for analyses
a15.core <- betapart.core(a15.5[5:46])

# returns three dissimilarity matrices containing
# pairwise between-site values of each beta-diversity component
# total (sor), turnover (sim), and nestedness (sne)
a15.t.dist <- beta.pair(a15.core, index.family = "sorensen")
str(a15.t.dist)

# run NMDS mdoel for total beta-diversity component (sor)
nmds.a15.t.sor <- metaMDS(a15.t.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.a15.t.sor
stressplot(nmds.a15.t.sor)
goodness(nmds.a15.t.sor)
nmds.a15.t.sor$stress
plot(nmds.a15.t.sor)

# plot the NMDS model
ordiplot(nmds.a15.t.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a15.t.sor, a15.5$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("bottomleft", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!


## Test for differences in ground beetle species composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(a15.t.dist$beta.sor ~ a15.5$Treatment, permutations = 999)
adonis2(a15.t.dist$beta.sor ~ a15.5$Canopy, permutations = 999)
adonis2(a15.t.dist$beta.sor ~ a15.5$Veg, permutations = 999)
pairwise.adonis(a15.t.dist$beta.sor, a15.5$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a15.tbeta.sor <- betadisper(a15.t.dist$beta.sor, a15.5$Treatment, type = c("median"))
a15.tbeta.sor
anova(a15.tbeta.sor)
plot(a15.tbeta.sor)
boxplot(a15.tbeta.sor, ylab = "Distance to median")
TukeyHSD(a15.tbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(a15.t.dist$beta.sim ~ a15.5$Treatment, permutations = 999)
adonis2(a15.t.dist$beta.sim ~ a15.5$Canopy, permutations = 999)
adonis2(a15.t.dist$beta.sim ~ a15.5$Veg, permutations = 999)

# nestedness
adonis2(a15.t.dist$beta.sne ~ a15.5$Treatment, permutations = 999)
adonis2(a15.t.dist$beta.sne ~ a15.5$Canopy, permutations = 999) # canopy gap nested within closed canopy
pairwise.adonis(a15.t.dist$beta.sne, a15.5$Canopy)
adonis2(a15.t.dist$beta.sne ~ a15.5$Veg, permutations = 999)

##################################################################################






## code for trait data and functional beta-diversity analyses
## will get to this soon, after taxonomic beta-diversity analyses
## Trait data
t <- read.csv("./bb_rtraits_v3.csv", row.names=1)
#trim the trait dataset by identifying traits that are highly correlated
#or lack sufficient variance among species
names(t)
str(t)
plot(t)

#Double check that all species are present in both datasets
intersect(colnames(a), rownames(t3))
intersect(colnames(a), p4$tip.label)
intersect(p4$tip.label, rownames(t3))
names(a)

#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))
setdiff(rownames (t3), p4$tip.label)
setdiff(colnames(a), p4$tip.label)

#Double check all species names are in the same order
rownames(t3) == colnames(a) 
colnames(a) == p4$tip.label
rownames(t3) == p4$tip.label

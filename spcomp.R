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

library(BiodiversityR)
citation("BiodiversityR")

library(betapart)
citation("betapart")

#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
citation("pairwiseAdonis")

# set up color and point vectors for figures
colvec1 <- c("darkorange1", "darkgoldenrod1", "black", "chartreuse4")
colvec2 <- c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF")
pchvec1 <- c(16, 17, 15, 18)
pchvec2 <- c(19, 15, 4, 9)
ltyvec <- c(1, 2, 3, 4)

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

# need to pool species data across sampling intervals
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
colSums(a13.4[5:52])

a13.4 <- a13.4[, colSums(a13.4 != 0) > 0]
colSums(a13.4[5:45])

## Compare ground beetle species richness among treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

levels(a13.4$Treatment)
rare.a13.C <- a13.4[which(a13.4$Treatment == "Control"),]
rare.a13.L <- a13.4[which(a13.4$Treatment == "Light"),]
rare.a13.V <- a13.4[which(a13.4$Treatment == "Veg"),]
rare.a13.LV <- a13.4[which(a13.4$Treatment == "Light+Veg"),]

sp.a13.C <- specaccum(rare.a13.C[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.L <- specaccum(rare.a13.L[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.V <- specaccum(rare.a13.V[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.LV <- specaccum(rare.a13.LV[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(sp.a13.C, pch = 19, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 300), ylim = c(0, 35))
plot(sp.a13.L, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 3, col = "#FDE725FF")
plot(sp.a13.V, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 3, col = "#73D055FF")
plot(sp.a13.LV, add = TRUE, pch = 9, xvar = c("individuals"), lty = 3, lwd = 3, col = "#2D708EFF")
legend("bottomright", legend = c("Canopy","Understory","Canopy+Understory", "Undisturbed"), 
       lty = ltyvec, bty = "n", lwd = 3, col = colvec2)

levels(a13.4$Treatment)

#calculates species richness for each site
specnumber(a13.4[5:45])

#calculates species richness by treatment
specnumber(a13.4[5:45], groups = a13.4$Treatment)
str(a13.4)

bark.sp.j1 <- diversitycomp(a13.4[5:45], y = a13.4, factor1 = "Canopy", factor2 = "Veg", index = "jack1") #this isn't working and I'm not sure why... yet.
bark.sp.j1
bark.sp.j2 <- diversitycomp(a13.4[5:45], y = a13.4, factor1 = "Treatment", index = "jack2")
bark.sp.j2

bark.j1 <- diversityresult(a13.4[5:45], y=NULL, index = "jack1")
bark.j1
bark.j2 <- diversityresult(a13.4[5:45], y=NULL, index = "jack2")
bark.j2


############################

## Partition beta-diversity of ground beetle communities among treatments

#change data set to presence/absence for this part of the analyses
a13.5 <- a13.4
a13.5[a13.5 > 0] <- 1

str(a13.5)
rowSums(a13.5[5:45])
colSums(a13.5[5:45])

# create beta part object for analyses
a13.core <- betapart.core(a13.5[5:45])

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

ordiplot(nmds.a13.sor, type="n", xlim = c(-1, 0.8), ylim = c(-0.6, 0.6))
points(nmds.a13.sor, display = "sites", pch = pchvec1[a13.5$Treatment], cex = 1.5, 
       col = colvec2[a13.5$Treatment])
ordiellipse(nmds.a13.sor, groups = a13.5$Treatment, display = "sites", draw = "lines", 
            col = colvec2[a13.5$Treatment], lwd = 3, conf = 0.90)
ordihull(nmds.a13.sor, groups = a13.5$Treatment, display = "sites", draw = c("polygon"), col = NULL,
         border = colvec2[a13.5$Treatment], lwd = 2.5)
legend("topleft", legend = c("Canopy","Understory","Canopy+Understory", "Undisturbed"), 
       pch = pchvec1[a13.5$Treatment], cex = 1.2, bty = "n", col = colvec2[a13.5$Treatment])
orditorp(nnmds.a13.sor, "sites") #used to double check the legend!

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

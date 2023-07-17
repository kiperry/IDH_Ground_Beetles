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

###################################################################################

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

###################################################################################
# need to remove extra control sites
a13.3.2 <- a13.3[1:24,]

# need to remove columns with zero values
# removes species not collected in these sites during this year
a13.4 <- a13.3.2[, colSums(a13.3.2 !=0) > 0]

colSums(a13.4[5:45])
# can use this format for community composition analyses

##################################################################################

## Partition beta-diversity of ground beetle communities among treatments

#change data set to presence/absence
a13.5 <- a13.4 # make a copy
a13.5[a13.5 > 0] <- 1

str(a13.5)
rowSums(a13.5[6:45])
colSums(a13.5[6:45])

# order the treatments so it makes sense in the figure below
a13.5$Treatment <- factor(a13.5$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))

## Taxonomic beta-diversity

# create beta part object for analyses
a13.core <- betapart.core(a13.5[6:45])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
# total (sor), turnover (sim), and nestedness (sne)
a13.t.dist <- beta.pair(a13.core, index.family = "sorensen")
str(a13.t.dist)

# run NMDS models for total beta-diversity component (sor)
nmds.a13.t.sor <- metaMDS(a13.t.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.a13.t.sor
stressplot(nmds.a13.t.sor)
goodness(nmds.a13.t.sor)
nmds.a13.t.sor$stress #stress is quality of fit
plot(nmds.a13.t.sor)

# plot the NMDS model
ordiplot(nmds.a13.t.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.7), ylim = c(-0.6, 0.8))
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a13.t.sor, a13.5$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("topright", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!

## Test for differences in ground beetle composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(a13.t.dist$beta.sor ~ a13.5$Treatment, permutations = 999)
adonis2(a13.t.dist$beta.sor ~ a13.5$Canopy, permutations = 999)
adonis2(a13.t.dist$beta.sor ~ a13.5$Veg, permutations = 999)
pairwise.adonis(a13.t.dist$beta.sor, a13.5$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a13.tbeta.sor <- betadisper(a13.t.dist$beta.sor, a13.5$Treatment, type = c("median"))
a13.tbeta.sor
anova(a13.tbeta.sor)
plot(a13.tbeta.sor)
boxplot(a13.tbeta.sor, ylab = "Distance to median")
TukeyHSD(a13.tbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(a13.t.dist$beta.sim ~ a13.5$Treatment, permutations = 999)
adonis2(a13.t.dist$beta.sim ~ a13.5$Canopy, permutations = 999)
adonis2(a13.t.dist$beta.sim ~ a13.5$Veg, permutations = 999)

# nestedness
adonis2(a13.t.dist$beta.sne ~ a13.5$Treatment, permutations = 999)
adonis2(a13.t.dist$beta.sne ~ a13.5$Canopy, permutations = 999)
adonis2(a13.t.dist$beta.sne ~ a13.5$Veg, permutations = 999)

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
# couldn't sample in sites when treatments were being established
a14.2 <- a14[37:252,]

summary(a14.2)

###################################################################################

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

###################################################################################
# need to remove extra control sites
a14.4.2 <- a14.4[1:24,]
colSums(a14.4.2[5:57])

# need to remove columns with zero values
# removes species not collected in these sites during this year
a14.5 <- a14.4.2[, colSums(a14.4.2 !=0) > 0]

# removed thirteen species
colSums(a14.5[5:44])
# can use this format for community composition analyses

###################################################################################

## Partition beta-diversity of ground beetle communities among treatments

# change data set to presence/absence
a14.6 <- a14.5 # make a copy
a14.6[a14.6 > 0] <- 1

str(a14.6)
rowSums(a14.6[5:44])
colSums(a14.6[5:44])

# create beta part object for analyses
a14.core <- betapart.core(a14.6[5:44])

# returns three dissimilarity matrices containing 
# pairwise between-site values of each beta-diversity component
# total (sor), turnover (sim), and nestedness (sne)
a14.t.dist <- beta.pair(a14.core, index.family = "sorensen")
str(a14.t.dist)

# run NMDS models for total beta-diversity component (sor)
nmds.a14.t.sor <- metaMDS(a14.t.dist$beta.sor, trymax = 500, autotransform = TRUE)
nmds.a14.t.sor
stressplot(nmds.a14.t.sor)
goodness(nmds.a14.t.sor)
nmds.a14.t.sor$stress #stress is quality of fit
plot(nmds.a14.t.sor)

# plot the NMDS model
ordiplot(nmds.a14.t.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a14.t.sor, a14.6$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("bottomleft", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!

## Test for differences in ground beetle composition across disturbance treatments
## for total beta-diversity component


# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(a14.t.dist$beta.sor ~ a14.6$Treatment, permutations = 999)
adonis2(a14.t.dist$beta.sor ~ a14.6$Canopy, permutations = 999)
adonis2(a14.t.dist$beta.sor ~ a14.6$Veg, permutations = 999)
pairwise.adonis(a14.t.dist$beta.sor, a14.6$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a14.tbeta.sor <- betadisper(a14.t.dist$beta.sor, a14.6$Treatment, type = c("median"))
a14.tbeta.sor
anova(a14.tbeta.sor)
plot(a14.tbeta.sor)
boxplot(a14.tbeta.sor, ylab = "Distance to median")
TukeyHSD(a14.tbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(a14.t.dist$beta.sim ~ a14.6$Treatment, permutations = 999)
adonis2(a14.t.dist$beta.sim ~ a14.6$Canopy, permutations = 999)
adonis2(a14.t.dist$beta.sim ~ a14.6$Veg, permutations = 999)

# nestedness
adonis2(a14.t.dist$beta.sne ~ a14.6$Treatment, permutations = 999)
adonis2(a14.t.dist$beta.sne ~ a14.6$Canopy, permutations = 999)
adonis2(a14.t.dist$beta.sne ~ a14.6$Veg, permutations = 999)

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

## Functional beta-diversity

# load the trait data
t <- read.csv("./PNR_Carabidae_Traits.csv", row.names = 1)

names(t)
str(t)
plot(t)

library(ggplot2)
library(GGally)
cor(t, method = c("pearson"), use = "complete.obs")
cp <- ggpairs(t, upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 14))

# remove traits that are highly correlated
t2 <- t[,-3] # minimum body length
t2 <- t2[,-3] # maximum body length

t2 <- t2[,-7] # antennae length

cp <- ggpairs(t2, upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 14))

# all correlations now below 0.70

# Create pairplot - for Appendix
png("Corplot_Traits.png", width = 1800, height = 1000, pointsize = 40)
cp <- ggpairs(t, upper = list(continuous = wrap("cor", size = 5, color = "black")))
cp + theme(strip.text.x = element_text(size = 23), strip.text.y = element_text(size = 20))
dev.off()

##################################################################################

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

##################################################################################

# now we have to match up the trait matrix with every abundance matrix

# start with 2013
intersect(colnames(a13.5[5:45]), rownames(t2))

setdiff(colnames(a13.5[5:45]), rownames(t2))
a13.fun <- a13.5[,-5]
a13.fun <- a13.fun[,-26]
setdiff(colnames(a13.fun[5:43]), rownames(t2))

setdiff(rownames(t2), colnames(a13.fun[5:43]))
t.13 <- t2[-3,]
t.13 <- t.13[-c(4:7),]
t.13 <- t.13[-c(5:6),]
t.13 <- t.13[-c(7:9),]
setdiff(rownames(t.13), colnames(a13.fun[5:43]))
t.13 <- t.13[-c(10:11),]
t.13 <- t.13[-c(14:16),]
t.13 <- t.13[-20,]
t.13 <- t.13[-24,]
t.13 <- t.13[-c(26:27),]
t.13 <- t.13[-29,]
t.13 <- t.13[-30,]
t.13 <- t.13[-31,]
t.13 <- t.13[-32,]
t.13 <- t.13[-32,]
t.13 <- t.13[-c(33:35),]
t.13 <- t.13[-34,]
t.13 <- t.13[-36,]
setdiff(rownames(t.13), colnames(a13.fun[5:43]))

intersect(colnames(a13.fun[5:43]), rownames(t.13))
rownames(t.13) == colnames(a13.fun[5:43]) 

##################################################################################

library(FD)
library(picante)
library(gawdis)

# community-weighed means
cwm.13 <- functcomp(t.13, as.matrix(a13.fun[5:43]), CWM.type = "all")
cwm.13
cwm.13.full <- cbind(a13.fun[1:4], cwm.13)
write.csv(cwm.13.full, file = "CWMs_2013.csv")

# functional beta-diversity
# weight body length traits and traits on the head to limit their total influence on the metric
# calculate the distance matrix

tdis.13 <- gawdis(t.13, w.type = "optimized", opti.maxiter = 300,
               groups.weight = T, groups = c(1, 2, 2, 3, 3, 4, 5))
attr(tdis.13, "correls")
attr(tdis.13, "weights")

pco.13 <- dudi.pco(sqrt(tdis.13), scannf = FALSE, nf = 4) # select four axes
scatter(pco.13)

pco.13$li
sum(pco.13$eig[1:4]) / sum(pco.13$eig) # 0.61
sum(pco.13$eig[1:3]) / sum(pco.13$eig) # 0.54
sum(pco.13$eig[1:2]) / sum(pco.13$eig) # 0.42

rowSums(a13.fun[5:43])

# let's use first four axes of PCoA for functional diversity metrics
t.ax.13 <- as.matrix(pco.13$li[1:4])

# returns pairwise between-site values of each functional beta-diversity component
fun.b.13 <- functional.beta.pair(a13.fun[5:43], t.ax.13, index.family = "sorensen")
str(fun.b.13)

# run NMDS mdoel for total beta-diversity component (sor)
nmds.a13.f.sor <- metaMDS(fun.b.13$funct.beta.sor, trymax = 500, autotransform = TRUE)
nmds.a13.f.sor
stressplot(nmds.a13.f.sor)
goodness(nmds.a13.f.sor)
nmds.a13.f.sor$stress
plot(nmds.a13.f.sor)

a13.fun$Treatment <- factor(a13.fun$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))

# plot the NMDS model
ordiplot(nmds.a13.f.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a13.f.sor, a13.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("bottomleft", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!


## Test for differences in ground beetle functional composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(fun.b.13$funct.beta.sor ~ a13.fun$Treatment, permutations = 999)
adonis2(fun.b.13$funct.beta.sor ~ a13.fun$Canopy, permutations = 999)
adonis2(fun.b.13$funct.beta.sor ~ a13.fun$Veg, permutations = 999)
pairwise.adonis(fun.b.13$funct.beta.sor, a13.fun$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a13.fbeta.sor <- betadisper(fun.b.13$funct.beta.sor, a13.fun$Treatment, type = c("median"))
a13.fbeta.sor
anova(a13.fbeta.sor)
plot(a13.fbeta.sor)
boxplot(a13.fbeta.sor, ylab = "Distance to median")
TukeyHSD(a13.fbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(fun.b.13$funct.beta.sim ~ a13.fun$Treatment, permutations = 999)
adonis2(fun.b.13$funct.beta.sim ~ a13.fun$Canopy, permutations = 999)
adonis2(fun.b.13$funct.beta.sim ~ a13.fun$Veg, permutations = 999)

# nestedness
adonis2(fun.b.13$funct.beta.sne ~ a13.fun$Treatment, permutations = 999)
adonis2(fun.b.13$funct.beta.sne ~ a13.fun$Canopy, permutations = 999) # canopy gap nested within closed canopy
pairwise.adonis(fun.b.13$funct.beta.sne, a13.fun$Canopy)
adonis2(fun.b.13$funct.beta.sne ~ a13.fun$Veg, permutations = 999)

##################################################################################

# now 2014
t.14 <- t2
intersect(colnames(a14.6[5:44]), rownames(t.14))

setdiff(colnames(a14.6[5:44]), rownames(t.14))
a14.fun <- a14.6[,-25]
setdiff(colnames(a14.fun[5:43]), rownames(t.14))

setdiff(rownames(t.14), colnames(a14.fun[5:43]))
t.14 <- t.14[-3,]
t.14 <- t.14[-c(4:7),]
t.14 <- t.14[-c(5:6),]
t.14 <- t.14[-c(7:9),]
t.14 <- t.14[-11,]
t.14 <- t.14[-c(14:17),]
t.14 <- t.14[-20,]
t.14 <- t.14[-21,]
t.14 <- t.14[-23,]
t.14 <- t.14[-26,]
t.14 <- t.14[-29,]
t.14 <- t.14[-32,]
t.14 <- t.14[-c(33:34),]
t.14 <- t.14[-34,]
t.14 <- t.14[-c(35:37),]
t.14 <- t.14[-38,]
t.14 <- t.14[-40,]
setdiff(rownames(t.14), colnames(a14.fun[5:43]))

intersect(colnames(a14.fun[5:43]), rownames(t.14))
rownames(t.14) == colnames(a14.fun[5:43])

##################################################################################

# community-weighed means
cwm.14 <- functcomp(t.14, as.matrix(a14.fun[5:43]), CWM.type = "all")
cwm.14
cwm.14.full <- cbind(a14.fun[1:4], cwm.14)
write.csv(cwm.14.full, file = "CWMs_2014.csv")

# functional beta-diversity
# weight body length traits and traits on the head to limit their total influence on the metric
# calculate the distance matrix

tdis.14 <- gawdis(t.14, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 2, 3, 3, 4, 5))
attr(tdis.14, "correls")
attr(tdis.14, "weights")

pco.14 <- dudi.pco(sqrt(tdis.14), scannf = FALSE, nf = 4) # select four axes
scatter(pco.14)

pco.14$li
sum(pco.14$eig[1:4]) / sum(pco.14$eig) # 0.60
sum(pco.14$eig[1:3]) / sum(pco.14$eig) # 0.54
sum(pco.14$eig[1:2]) / sum(pco.14$eig) # 0.40

rowSums(a14.fun[5:43])

# let's use first four axes of PCoA for functional diversity metrics
t.ax.14 <- as.matrix(pco.14$li[1:4])

# returns pairwise between-site values of each functional beta-diversity component
fun.b.14 <- functional.beta.pair(a14.fun[5:43], t.ax.14, index.family = "sorensen")
str(fun.b.14)

# run NMDS mdoel for total beta-diversity component (sor)
nmds.a14.f.sor <- metaMDS(fun.b.14$funct.beta.sor, trymax = 500, autotransform = TRUE)
nmds.a14.f.sor
stressplot(nmds.a14.f.sor)
goodness(nmds.a14.f.sor)
nmds.a14.f.sor$stress
plot(nmds.a14.f.sor)

a14.fun$Treatment <- factor(a14.fun$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))

# plot the NMDS model
ordiplot(nmds.a14.f.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a14.f.sor, a14.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("bottomleft", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!


## Test for differences in ground beetle functional composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(fun.b.14$funct.beta.sor ~ a14.fun$Treatment, permutations = 999)
adonis2(fun.b.14$funct.beta.sor ~ a14.fun$Canopy, permutations = 999)
adonis2(fun.b.14$funct.beta.sor ~ a14.fun$Veg, permutations = 999)
pairwise.adonis(fun.b.14$funct.beta.sor, a14.fun$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a14.fbeta.sor <- betadisper(fun.b.14$funct.beta.sor, a14.fun$Treatment, type = c("median"))
a14.fbeta.sor
anova(a14.fbeta.sor)
plot(a14.fbeta.sor)
boxplot(a14.fbeta.sor, ylab = "Distance to median")
TukeyHSD(a14.fbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(fun.b.14$funct.beta.sim ~ a14.fun$Treatment, permutations = 999)
adonis2(fun.b.14$funct.beta.sim ~ a14.fun$Canopy, permutations = 999)
adonis2(fun.b.14$funct.beta.sim ~ a14.fun$Veg, permutations = 999)

# nestedness
adonis2(fun.b.14$funct.beta.sne ~ a14.fun$Treatment, permutations = 999)
adonis2(fun.b.14$funct.beta.sne ~ a14.fun$Canopy, permutations = 999) # canopy gap nested within closed canopy
pairwise.adonis(fun.b.14$funct.beta.sne, a14.fun$Canopy)
adonis2(fun.b.14$funct.beta.sne ~ a14.fun$Veg, permutations = 999)

##################################################################################

# now 2015
t.15 <- t2
intersect(colnames(a15.5[5:46]), rownames(t.15))

setdiff(colnames(a15.5[5:46]), rownames(t.15))
a15.fun <- a15.5[,-25]
setdiff(colnames(a15.fun[5:45]), rownames(t.15))

setdiff(rownames(t.15), colnames(a15.fun[5:45]))
t.15 <- t.15[-1,]
t.15 <- t.15[-c(5:7),]
t.15 <- t.15[-7,]
t.15 <- t.15[-c(9:11),]
t.15 <- t.15[-c(12:13),]
t.15 <- t.15[-c(15:17),]
t.15 <- t.15[-16,]
t.15 <- t.15[-19,]
t.15 <- t.15[-20,]
t.15 <- t.15[-24,]
t.15 <- t.15[-27,]
t.15 <- t.15[-30,]
t.15 <- t.15[-35,]
t.15 <- t.15[-c(37:39),]
t.15 <- t.15[-c(41:42),]
t.15 <- t.15[-c(42:43),]

intersect(colnames(a15.fun[5:45]), rownames(t.15))
rownames(t.15) == colnames(a15.fun[5:45])

##################################################################################

# community-weighed means
cwm.15 <- functcomp(t.15, as.matrix(a15.fun[5:45]), CWM.type = "all")
cwm.15
cwm.15.full <- cbind(a15.fun[1:4], cwm.15)
write.csv(cwm.15.full, file = "CWMs_2015.csv")

# functional beta-diversity
# weight body length traits and traits on the head to limit their total influence on the metric
# calculate the distance matrix

tdis.15 <- gawdis(t.15, w.type = "optimized", opti.maxiter = 300,
                  groups.weight = T, groups = c(1, 2, 2, 3, 3, 4, 5))
attr(tdis.15, "correls")
attr(tdis.15, "weights")

pco.15 <- dudi.pco(sqrt(tdis.15), scannf = FALSE, nf = 4) # select four axes
scatter(pco.15)

pco.15$li
sum(pco.15$eig[1:4]) / sum(pco.15$eig) # 0.62
sum(pco.15$eig[1:3]) / sum(pco.15$eig) # 0.56
sum(pco.15$eig[1:2]) / sum(pco.15$eig) # 0.43

rowSums(a15.fun[5:45])

# let's use first four axes of PCoA for functional diversity metrics
t.ax.15 <- as.matrix(pco.15$li[1:4])

# returns pairwise between-site values of each functional beta-diversity component
fun.b.15 <- functional.beta.pair(a15.fun[5:45], t.ax.15, index.family = "sorensen")
str(fun.b.15)

# run NMDS mdoel for total beta-diversity component (sor)
nmds.a15.f.sor <- metaMDS(fun.b.15$funct.beta.sor, trymax = 500, autotransform = TRUE)
nmds.a15.f.sor
stressplot(nmds.a15.f.sor)
goodness(nmds.a15.f.sor)
nmds.a15.f.sor$stress
plot(nmds.a15.f.sor)

a15.fun$Treatment <- factor(a15.fun$Treatment, levels = c("Light", "Veg", "Light+Veg", "Control"))

# plot the NMDS model
ordiplot(nmds.a15.f.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a15.f.sor, a15.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("bottomleft", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))

#orditorp(nmds.a15.t.sor, "sites") #used to double check the legend!


## Test for differences in ground beetle functional composition across disturbance treatments
## for total beta-diversity component

# PERMANOVA tests whether the group centroid of beetle communities among treatments
# differs in multivariate space (e.g. different community composition)
adonis2(fun.b.15$funct.beta.sor ~ a15.fun$Treatment, permutations = 999)
adonis2(fun.b.15$funct.beta.sor ~ a15.fun$Canopy, permutations = 999)
adonis2(fun.b.15$funct.beta.sor ~ a15.fun$Veg, permutations = 999)
pairwise.adonis(fun.b.15$funct.beta.sor, a15.fun$Treatment)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogenetiy of variances
a15.fbeta.sor <- betadisper(fun.b.15$funct.beta.sor, a15.fun$Treatment, type = c("median"))
a15.fbeta.sor
anova(a15.fbeta.sor)
plot(a15.fbeta.sor)
boxplot(a15.fbeta.sor, ylab = "Distance to median")
TukeyHSD(a15.fbeta.sor, which = "group", conf.level = 0.95)

# PERMANOVA tests for turnover and nestedness

# turnover
adonis2(fun.b.15$funct.beta.sim ~ a15.fun$Treatment, permutations = 999)
adonis2(fun.b.15$funct.beta.sim ~ a15.fun$Canopy, permutations = 999)
adonis2(fun.b.15$funct.beta.sim ~ a15.fun$Veg, permutations = 999)

# nestedness
adonis2(fun.b.15$funct.beta.sne ~ a15.fun$Treatment, permutations = 999)
adonis2(fun.b.15$funct.beta.sne ~ a15.fun$Canopy, permutations = 999) # canopy gap nested within closed canopy
pairwise.adonis(fun.b.15$funct.beta.sne, a15.fun$Canopy)
adonis2(fun.b.15$funct.beta.sne ~ a15.fun$Veg, permutations = 999)

##################################################################################

# Let's make a panel plot for taxonomic and functional beta-diversity

png("NMDS_Betadiversity.png", width = 1800, height = 2500, pointsize = 30)

par(mfrow=c(3,2))
par(mar=c(5,8,4,2))

# 2013 - taxonomic
ordiplot(nmds.a13.t.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.7), ylim = c(-0.6, 0.8),
         main = "Taxonomic", cex.main = 2)
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a13.t.sor, dis = "sites", select = which(a13.5$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a13.t.sor, a13.5$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)
text(-0.7, 0.65, "A", pos = 3, font = 2, cex = 2)


# 2013 - functional
ordiplot(nmds.a13.f.sor, disp = "sites", type = "n", xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.95),
         main = "Functional", cex.main = 2)
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a13.f.sor, dis = "sites", select = which(a13.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a13.f.sor, a13.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)

legend("topright", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), cex = 1.5, bty = "n", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))
text(-0.78, 0.79, "B", pos = 3, font = 2, cex = 2)


# 2014 - taxonomic
ordiplot(nmds.a14.t.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.7), ylim = c(-0.4, 0.4))
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a14.t.sor, dis = "sites", select = which(a14.6$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a14.t.sor, a14.6$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)
text(-0.65, 0.515, "C", pos = 3, font = 2, cex = 2)


# 2014 - functional
ordiplot(nmds.a14.f.sor, disp = "sites", type = "n", xlim = c(-0.7, 0.8), ylim = c(-0.4, 0.4))
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a14.f.sor, dis = "sites", select = which(a14.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a14.f.sor, a14.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)
text(-0.65, 0.55, "D", pos = 3, font = 2, cex = 2)


# 2015 - taxonomic
ordiplot(nmds.a15.t.sor, disp = "sites", type = "n", xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.65))
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a15.t.sor, dis = "sites", select = which(a15.5$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a15.t.sor, a15.5$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)
text(-0.75, 0.63, "E", pos = 3, font = 2, cex = 2)


# 2015 - functional
ordiplot(nmds.a15.f.sor, disp = "sites", type = "n", xlim = c(-0.6, 0.6), ylim = c(-0.4, 0.4))
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Light"), pch = 16, cex = 2, col = "#FDE725FF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Veg"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Light+Veg"), pch = 15, cex = 2, col = "#2D708EFF")
points(nmds.a15.f.sor, dis = "sites", select = which(a15.fun$Treatment=="Control"), pch = 18, cex = 2, col = "#481567FF")

ordiellipse(nmds.a15.f.sor, a15.fun$Treatment, draw = "lines", col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"), 
            lwd = 3, kind = "sd", conf = 0.95, label = FALSE)
text(-0.55, 0.45, "F", pos = 3, font = 2, cex = 2)

dev.off()



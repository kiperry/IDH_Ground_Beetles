###################################################################################
#
# PNR ground beetle data
#
# Changes in ground beetle species richness and diversity among disturbance treatments
# Years 2013, 2014, 2015
#
# Rarefaction and jackknife estimates
#
# KI Perry; 10 April 2023
#
###################################################################################

## Load necessary packages

library(reshape2)
library(ggplot2)
library(viridis)
library(vegan)
library(fossil)

#library(iNEXT)
#citation("iNEXT")

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

# need to pool species data across sampling intervals for community analyses
a13.2 <- melt(a13, id = c("Interval", "Quadrat", "Canopy", "Veg",
                          "Treatment", "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a13.2)

# double check only species in the variable column
levels(a13.2$variable)
names(a13.2)[10] <- "Species"

# now create data frame for rarefaction analyses
a13.3 <- dcast(a13.2, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same
colSums(a13.3[5:52])
colSums(a13.3[5:52]) == colSums(a13[10:57])

str(a13.3)

################################################################################
# need to remove extra control sites
a13.3.2 <- a13.3[1:24,]
colSums(a13.3.2[5:52])

# need to remove columns with zero values
# removes species not collected in these sites during this year
a13.4 <- a13.3.2[, colSums(a13.3.2 !=0) > 0]

colSums(a13.4[5:45])
# can use this format for community composition analyses

################################################################################

## Compare ground beetle species richness among treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

levels(a13.4$Treatment)
rare.a13.L <- a13.4[which(a13.4$Treatment == "Light"),]
rare.a13.V <- a13.4[which(a13.4$Treatment == "Veg"),]
rare.a13.LV <- a13.4[which(a13.4$Treatment == "Light+Veg"),]
rare.a13.C <- a13.4[which(a13.4$Treatment == "Control"),]

sp.a13.L <- specaccum(rare.a13.L[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.V <- specaccum(rare.a13.V[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.LV <- specaccum(rare.a13.LV[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a13.C <- specaccum(rare.a13.C[5:45], method = "rarefaction", permutations = 100, gamma = "jack2")

# make the plot

plot(sp.a13.C, pch = 19, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 320), ylim = c(0, 32))
plot(sp.a13.L, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 3, col = "#FDE725FF")
plot(sp.a13.V, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 3, col = "#73D055FF")
plot(sp.a13.LV, add = TRUE, pch = 9, xvar = c("individuals"), lty = 3, lwd = 3, col = "#2D708EFF")

legend("bottomright", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1.5, bty = "n", lwd = 3,
       col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))


# calculates species richness for each site
specnumber(a13.4[5:45])

# calculates species richness by treatment
specnumber(a13.4[5:45], groups = a13.4$Treatment)

# calculates the estimated species richness of a population using first- and second-order
# jackknife estimators

# start with the entire data set
jack1(a13.4[5:45], taxa.row = FALSE, abund = TRUE)# essentially gamma diversity
jack2(a13.4[5:45], taxa.row = FALSE, abund = TRUE)

# canopy
jack1(rare.a13.L[5:45], taxa.row = FALSE, abund = TRUE)
jack2(rare.a13.L[5:45], taxa.row = FALSE, abund = TRUE)

# understory
jack1(rare.a13.V[5:45], taxa.row = FALSE, abund = TRUE)
jack2(rare.a13.V[5:45], taxa.row = FALSE, abund = TRUE)

# canopy + understory
jack1(rare.a13.LV[5:45], taxa.row = FALSE, abund = TRUE)
jack2(rare.a13.LV[5:45], taxa.row = FALSE, abund = TRUE)

# undisturbed
jack1(rare.a13.C[5:45], taxa.row = FALSE, abund = TRUE)
jack2(rare.a13.C[5:45], taxa.row = FALSE, abund = TRUE)

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
rowSums(a14[10:62])
colSums(a14[10:62])
head(a14)
tail(a14)
dim(a14)

summary(a14)

# need to remove the first sampling interval with all the NAs
# couldn't sample in sites when treatments were being established
a14.1 <- a14[37:252,]

str(a14.1)
summary(a14.1)

################################################################################

# need to pool species data across sampling intervals for community analyses
a14.2 <- melt(a14.1, id = c("Interval", "Quadrat", "Canopy", "Veg",
                          "Treatment", "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a14.2)

# double check only species in the variable column
levels(a14.2$variable)
names(a14.2)[10] <- "Species"

# now create data frame for rarefaction analyses
a14.3 <- dcast(a14.2, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same
colSums(a14.3[5:57])
colSums(a14.3[5:57]) == colSums(a14.1[10:62])

str(a14.3)

################################################################################
# need to remove extra control sites
a14.3.2 <- a14.3[1:24,]
colSums(a14.3.2[5:57])

# need to remove columns with zero values
# removes species not collected in these sites during this year
a14.4 <- a14.3.2[, colSums(a14.3.2 !=0) > 0]

colSums(a14.4[5:44])
# can use this format for community composition analyses

################################################################################

## Compare ground beetle species richness among treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

levels(a14.4$Treatment)
rare.a14.L <- a14.4[which(a14.4$Treatment == "Light"),]
rare.a14.V <- a14.4[which(a14.4$Treatment == "Veg"),]
rare.a14.LV <- a14.4[which(a14.4$Treatment == "Light+Veg"),]
rare.a14.C <- a14.4[which(a14.4$Treatment == "Control"),]

sp.a14.L <- specaccum(rare.a14.L[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.V <- specaccum(rare.a14.V[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.LV <- specaccum(rare.a14.LV[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a14.C <- specaccum(rare.a14.C[5:44], method = "rarefaction", permutations = 100, gamma = "jack2")

# make the plot

plot(sp.a14.C, pch = 19, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 450), ylim = c(0, 32))
plot(sp.a14.L, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 3, col = "#FDE725FF")
plot(sp.a14.V, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 3, col = "#73D055FF")
plot(sp.a14.LV, add = TRUE, pch = 9, xvar = c("individuals"), lty = 3, lwd = 3, col = "#2D708EFF")

legend("bottomright", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1.5, bty = "n", lwd = 3,
       col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))


# calculates species richness for each site
specnumber(a14.4[5:44])

# calculates species richness by treatment
specnumber(a14.4[5:44], groups = a14.4$Treatment)

# calculates the estimated species richness of a population using first- and second-order
# jackknife estimators

# start with the entire data set
jack1(a14.4[5:44], taxa.row = FALSE, abund = TRUE) # essentially gamma diversity
jack2(a14.4[5:44], taxa.row = FALSE, abund = TRUE)

# canopy
jack1(rare.a14.L[5:44], taxa.row = FALSE, abund = TRUE)
jack2(rare.a14.L[5:44], taxa.row = FALSE, abund = TRUE)

# understory
jack1(rare.a14.V[5:44], taxa.row = FALSE, abund = TRUE)
jack2(rare.a14.V[5:44], taxa.row = FALSE, abund = TRUE)

# canopy + understory
jack1(rare.a14.LV[5:44], taxa.row = FALSE, abund = TRUE)
jack2(rare.a14.LV[5:44], taxa.row = FALSE, abund = TRUE)

# undisturbed
jack1(rare.a14.C[5:44], taxa.row = FALSE, abund = TRUE)
jack2(rare.a14.C[5:44], taxa.row = FALSE, abund = TRUE)


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

################################################################################

# need to pool species data across sampling intervals for community analyses
a15.2 <- melt(a15, id = c("Interval", "Quadrat", "Canopy", "Veg",
                            "Treatment", "Disturbance", "DateSet", "DateColl", "Trap_Days"))
str(a15.2)

# double check only species in the variable column
levels(a15.2$variable)
names(a15.2)[10] <- "Species"

# now create data frame for rarefaction analyses
a15.3 <- dcast(a15.2, Quadrat + Treatment + Canopy + Veg ~ Species, sum)

# double check the numbers are the same
colSums(a15.3[5:61])
colSums(a15.3[5:61]) == colSums(a15[10:66])

str(a15.3)

################################################################################
# need to remove extra control sites
a15.3.2 <- a15.3[1:24,]
colSums(a15.3.2[5:61])

# need to remove columns with zero values
# removes species not collected in these sites during this year
a15.4 <- a15.3.2[, colSums(a15.3.2 !=0) > 0]

colSums(a15.4[5:46])
# can use this format for community composition analyses

################################################################################

## Compare ground beetle species richness among treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

levels(a15.4$Treatment)
rare.a15.L <- a15.4[which(a15.4$Treatment == "Light"),]
rare.a15.V <- a15.4[which(a15.4$Treatment == "Veg"),]
rare.a15.LV <- a15.4[which(a15.4$Treatment == "Light+Veg"),]
rare.a15.C <- a15.4[which(a15.4$Treatment == "Control"),]

sp.a15.L <- specaccum(rare.a15.L[5:46], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a15.V <- specaccum(rare.a15.V[5:46], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a15.LV <- specaccum(rare.a15.LV[5:46], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a15.C <- specaccum(rare.a15.C[5:46], method = "rarefaction", permutations = 100, gamma = "jack2")

# make the plot

plot(sp.a15.C, pch = 19, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 300), ylim = c(0, 35))
plot(sp.a15.L, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 3, col = "#FDE725FF")
plot(sp.a15.V, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 3, col = "#73D055FF")
plot(sp.a15.LV, add = TRUE, pch = 9, xvar = c("individuals"), lty = 3, lwd = 3, col = "#2D708EFF")

legend("bottomright", legend = c("Canopy", "Understory", "Canopy+Understory", "Undisturbed"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1.5, bty = "n", lwd = 3,
       col = c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF"))


# calculates species richness for each site
specnumber(a15.4[5:46])

# calculates species richness by treatment
specnumber(a15.4[5:46], groups = a15.4$Treatment)

# calculates the estimated species richness of a population using first- and second-order
# jackknife estimators

# start with the entire data set
jack1(a15.4[5:46], taxa.row = FALSE, abund = TRUE) # essentially gamma diversity
jack2(a15.4[5:46], taxa.row = FALSE, abund = TRUE)

# canopy
jack1(rare.a15.L[5:46], taxa.row = FALSE, abund = TRUE)
jack2(rare.a15.L[5:46], taxa.row = FALSE, abund = TRUE)

# understory
jack1(rare.a15.V[5:46], taxa.row = FALSE, abund = TRUE)
jack2(rare.a15.V[5:46], taxa.row = FALSE, abund = TRUE)

# canopy + understory
jack1(rare.a15.LV[5:46], taxa.row = FALSE, abund = TRUE)
jack2(rare.a15.LV[5:46], taxa.row = FALSE, abund = TRUE)

# undisturbed
jack1(rare.a15.C[5:46], taxa.row = FALSE, abund = TRUE)
jack2(rare.a15.C[5:46], taxa.row = FALSE, abund = TRUE)

## Create panel figure


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

# set up color and point vectors for figures
colvec1 <- c("darkorange1", "darkgoldenrod1", "black", "chartreuse4")
colvec2 <- c("#FDE725FF", "#73D055FF", "#2D708EFF", "#481567FF")
pchvec1 <- c(16, 17, 15, 18)
pchvec2 <- c(19, 15, 4, 9)
ltyvec <- c(1, 2, 3, 4)

library(reshape2)

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

library(BiodiversityR)
citation("BiodiversityR")

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

bark.sp.j1 <- diversitycomp(a13.4[5:45], y = a13.4, factor1 = "Treatment", index = "jack1")
bark.sp.j1
bark.sp.j2 <- diversitycomp(a13.4[5:45], y = a13.4, factor1 = "Treatment", index = "jack2")
bark.sp.j2

bark.j1 <- diversityresult(a13.4[5:45], y=NULL, index = "jack1")
bark.j1
bark.j2 <- diversityresult(a13.4[5:45], y=NULL, index = "jack2")
bark.j2


#########################################################

## Partition beta-diversity of ground beetle communities among treatments

#change data set to presence/absence for this part of the analyses
a13.5 <- a13.4
a13.5[a13.5 > 0] <- 1

str(a13.5)
rowSums(a13.5[5:45])
colSums(a13.5[5:45])







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

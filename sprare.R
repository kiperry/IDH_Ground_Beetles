###################################################################################
#
# PNR ground beetle data
#
# Changes in ground beetle species richness and diversity among disturbance treatments
# Years 2013, 2014, 2015
#
# Rarefaction and extrapolation
# Hills numbers
#
# KI Perry; 10 April 2023
#
###################################################################################

## Load necessary packages

library(reshape2)

library(ggplot2)

library(viridis)

library(iNEXT)
citation("iNEXT")

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
colSums(a13.4[5:52])

a13.4 <- a13.4[, colSums(a13.4 != 0) > 0]
colSums(a13.4[5:45]) #can use this format for community composition analyses

# need to revise format for rarefaction analyses

a13.5 <- melt(a13.4, id = c("Treatment", "Quadrat", "Canopy", "Veg"))
str(a13.5)
# double check only species in the variable column
levels(a13.5$variable)
names(a13.5)[5] <- "Species"
a13.6 <- dcast(a13.5, Species ~ Treatment, sum)
str(a13.6)
a13.7 <- a13.6[,-1]
rownames(a13.7) <- a13.6[,1]
a13.7 <- a13.7[,c(1,2,4,3)] #change column order
str(a13.7)

## Compare ground beetle species richness among treatments
r13 <- iNEXT(a13.7, q = c(0, 1, 2), datatype = "abundance", nboot = 100)
r13 <- iNEXT(a13.7, q = 1, datatype = "abundance", nboot = 100)
r13 <- iNEXT(a13.7, q = 2, datatype = "abundance", nboot = 100)

r13$DataInfo
r13$iNextEst
r13$AsyEst

ggiNEXT(r13, type = 1, se = TRUE, facet.var = "Assemblage",
        color.var = "Assemblage", grey = TRUE) +
  theme_bw(base_size = 18) +
  theme(legend.position = "right")

ggiNEXT(r13, type = 1, se = TRUE, facet.var = "Order", color.var = "Assemblage", grey = FALSE) +
  theme_bw(base_size = 18) + scale_fill_viridis(alpha = 0.7, discrete = TRUE, option = "D") +
  scale_color_viridis(discrete = TRUE, option = "D") + ylab("Species diversity metrics")

ggiNEXT(r13, type = 2, se = TRUE, facet.var = "None", color.var = "Assemblage", grey = FALSE)
ggiNEXT(r13, type = 3, se = TRUE, facet.var = "None", color.var = "Assemblage", grey = FALSE)

#the confidence intervals of any standardized diversity are obtained by a bootstrap method.
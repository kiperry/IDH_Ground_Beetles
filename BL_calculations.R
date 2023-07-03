###################################################################################
#
# PNR ground beetle data
#
# Body length range calculations
# Years 2013, 2014, 2015 combined
#
#
# KI Perry; 3 July 2023
#
###################################################################################

# load the data
bl <- read.csv("./PNR_Carabidae_BodyLength.csv")
str(bl)

# change species to factor
bl$Species <- as.factor(bl$Species)

# change length to numeric
bl$Length <- as.numeric(bl$Length)

# calculate average body length by species
aggregate(bl$Length, list(bl$Species), FUN = mean, na.rm = TRUE)

# calculate body length range by species
aggregate(bl$Length, list(bl$Species), FUN = range, na.rm = TRUE)

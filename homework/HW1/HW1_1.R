##
##  HW1_1.r
##
##  Investigate the full PUMS dataset to find demographic features that could uniquely identify individuals for a reidentification attack.
##
##  BP 2.18.2019
##

#read in the PUMS
data <- read.csv("../../data/FultonPUMS5full.csv");

puma_regions <- unique(data$puma);


##
##  HW1_1.r
##
##  Investigate the full PUMS dataset to find demographic features that could uniquely identify individuals for a reidentification attack.
##
##  BP 2.18.2019
##

rm(list=ls()); #clear workspace of all variables
library(plyr);
library(dplyr);

#read in the PUMS
data <- read.csv("../../data/FultonPUMS5full.csv");

#groupby PUMA region, age, sex, and race in data
data_gb <- data %>% group_by(puma, sex, age, latino, black, asian) %>% summarise(n=n());
#get rows of groupby where we only have unique individuals
unique_indivs <- data_gb[data_gb$n == 1, ];
print(nrow(unique_indivs))
#calculate the percent of people that I could uniquely
perc_unq <- nrow(unique_indivs) / nrow(data);
print(perc_unq)

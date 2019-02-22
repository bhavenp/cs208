##
##  HW1_1.r
##
##  Investigate the full PUMS dataset to find demographic features that could uniquely identify individuals for a reidentification attack.
##
##  BP 2.18.2019
##

rm(list=ls()); #clear workspace of all variables

#read in the PUMS
data <- read.csv("../../data/FultonPUMS5full.csv");

puma_regions <- unique(data$puma);
puma_regions_table <- table(data$puma);
puma_regions_table

# puma_1101 <- data[data$puma == 1101, ];
dups <- data[duplicated(data), ]


ages = c();
education = c();
num_unique = c();
for(a in 18:100){
  for(e in 1:16){
    query <- data[data$puma == 1101 & data$age==a & data$educ==e, ];
    num_unique <- c( num_unique, nrow(query));
    ages <- c(ages, a);
    education <- c(education, e);
  }
}
query_df <- data.frame("Age" = ages, "Education" = education, "Unique" = num_unique);

unq_query_df <- query_df[query_df$Unique == 1, ]

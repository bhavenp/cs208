##
##  HW1_3.r
##
##  Implement a membership attack on the sample of 100 points from the PUMS dataset
##
##  BP 2.18.2019
##

rm(list=ls());
library(plyr); #import this library for doing rounding for noise
library(ggplot2); #import library for plotting
library(gridExtra);

#### Parameters ####
prime <- 113; # prime number for hashing creating random vectors to pass to query
n <- 100;        # Dataset size
num_attr <- n*n; #calculate the number of attributes needed
noise_type <- "Subsampling"; # What type of noise will be used as defense. Can be "Rounding", "Gaussian", or "Subsampling"
# noise_val <- 100; #noise to introduce for Rounding
# noise_val <- 11; #noise to introduce for Gaussian
noise_val <- 50; #noise to introduce for Subsampling


#### Import Data ####
pums_full <- read.csv(file = "../../data/FultonPUMS5full.csv"); #read in the full data
#trim the population to PUB categories
pums_full_trimmed <- pums_full[, c("sex", "age", "educ", "latino", "black", "asian",  "married", "divorced", "children", "employed", "militaryservice", "disability", "englishability")];
pums_100 <- read.csv(file="../../data/FultonPUMS5sample100.csv"); #read in sample data from data folder
#trim sample to only PUB cols
pums_100_trimmed <- pums_100[, c("sex", "age", "educ", "latino", "black", "asian",  "married", "divorced", "children", "employed", "militaryservice", "disability", "englishability")];
#--------------------------------------------------------------------#

#### Do some cleaning of the data ####
#dummify age by comparing each individual's age to the mean
pums_full_trimmed$age <- ifelse(pums_full_trimmed$age > mean(pums_full_trimmed$age), 1, 0);
pums_100_trimmed$age <- ifelse(pums_100_trimmed$age > mean(pums_full_trimmed$age), 1, 0);

#dummify education column in the population and sample
for(e in unique(pums_full_trimmed$educ)){
  pums_full_trimmed[paste("educ", e, sep = "_")] <- ifelse(pums_full_trimmed$educ == e, 1, 0);
  pums_100_trimmed[paste("educ", e, sep = "_")] <- ifelse(pums_100_trimmed$educ == e, 1, 0);
}
#remove the education column in the population and sample
pums_full_trimmed <- subset(pums_full_trimmed, select = -c(educ));
pums_100_trimmed <- subset(pums_100_trimmed, select = -c(educ)); #remove the education column

#generate random predicates
preds_to_add <- as.matrix( replicate(num_attr-ncol(pums_full_trimmed), sample(0:prime-1, size=ncol(pums_full_trimmed)), simplify=TRUE) );

#generate binary values for each row in sample using predicates
preds_to_add_pop <- (as.matrix(pums_full_trimmed) %*% preds_to_add) %% prime %% 2;
#generate binary values for each row in sample using predicates
preds_to_add_samp <- (as.matrix(pums_100_trimmed) %*% preds_to_add) %% prime %% 2;

#create final matrices, which represent our population and our sample of 100 with 10,000 attributes
pums_full_final <- cbind(as.matrix(pums_full_trimmed), preds_to_add_pop);
pums_100_final <- cbind(as.matrix(pums_100_trimmed), preds_to_add_samp);

## Generate underlying population attributes
pop_prob <- colMeans(pums_full_final);
#--------------------------------------------------------------------#

#### Bulid Null Distribution ####
## A utility function to create data from the population
rmvbernoulli <- function(n=1, prob){
  history <- matrix(NA, nrow=n, ncol=length(prob))
  for(i in 1:n){
    x<- rbinom(n=length(prob), size=1, prob=prob); #
    x[x==0] <- -1      # Placeholder for transformation
    history[i,] <- x
  }
  return(history)
}

#### function to generate noisy sample means from data
gen_sample_probs <- function(data, noise_type, noise_param){
  if(noise_type == "Rounding"){
    col_sums = colSums(data);
    noisy_means <- round_any(col_sums, noise_param) / nrow(data); #this function will round the sums according to a multiple of noise_param and divide by the number of data points
    return( 2*(noisy_means - 0.5) );
  }
  else if(noise_type == "Gaussian"){
    noisy_sums <- colSums(data) + rnorm(n=ncol(data), mean=0, sd=noise_param); #add noise sampled by from N(0, noise_param)
    noisy_means <- noisy_sums / nrow(data);
    return( 2*(noisy_means - 0.5) );
  }
  else{
    subsamp_ind <- sample(x=1:nrow(data), size=noise_param, replace=FALSE); #get a random sample of the indices w/o replacement
    noisy_means <- colMeans(data[subsamp_ind, ]);  #get column means from subsample
    return( 2*(noisy_means - 0.5) );
  }
}

## A null distribution and critical value generator. Taken from membershipAttack.r
nullDistribution <- function(null.sims=1000, alpha=0.05, test_stat, population.prob, sample_means){
  population.mean <- 2*(population.prob-0.5)
  hold <- rep(NA,null.sims);
  
  for(i in 1:null.sims){
    nullAlice <- rmvbernoulli(n=1, prob=population.prob); #get an Alice that is just from the population
    hold[i] <- eval(test_stat(alice=nullAlice, sample.mean=sample_means, population.mean=population.mean));
  }
  nullDistribution <- sort(hold, decreasing=TRUE);
  criticalValue <- nullDistribution[round(alpha*null.sims)];
  return(list(nullDist=nullDistribution, criticalVal=criticalValue));
}

#function that defines the Dwork test statistic. Taken from membershipAttack.r
test.Dwork <- function(alice, sample.mean, population.mean){
  test.statistic <- sum(alice * sample.mean) - sum(population.mean * sample.mean);
  return(test.statistic)
}
#--------------------------------------------------------------------#

#### Perform Attack ####
range_d = c(1, 5, seq(50, 10000, by=50)); #range of number of attributes to go through

history <- matrix(NA, nrow=length(range_d), ncol=3);
myalpha <- 1 / (10*n); #myalpha is 0.001
population.mean <- 2*(pop_prob - 0.5); #scale population means/probs to -1 and 1
#get noisy sample means for all 10,000 attributes
sample_means <- gen_sample_probs(pums_100_final, noise_type=noise_type, noise_param=noise_val);

print(Sys.time())
#loop through the number of attributes
row_counter = 1;
for(d in range_d ){ #loop through the different number of attributes
  print(d);
  sample_means_d <- sample_means[1:d]; #cut sample means to 1:d attributes
  #generate a new test statistic for a new null distribution based on 'd' attributes
  output <- nullDistribution(alpha = myalpha, test_stat = test.Dwork,
                             population.prob = pop_prob[1:d], sample_means = sample_means_d);
  # test_stats_for_samp <- c();
  true_pos = 0;
  for(a in 1:nrow(pums_100_final)){ #loop through the 100 people in the sample
    alice <- pums_100_final[a, 1:d]; #choose Alice from sample. Cut to 1:d columns
    alice[alice==0] <- -1; #change 0's for alice to -1's
    # Conduct test for test statistic
    test.alice.Dwork <- test.Dwork(alice=alice, sample.mean=sample_means_d,
                                   population.mean=population.mean[1:d]);
    # test_stats_for_samp <- c(test_stats_for_samp, test.alice.Dwork); #for debugging  
    if( test.alice.Dwork >= output$criticalVal){ #check if test statistic is greater than critical value
      true_pos = true_pos + 1;
    }
  }
  tp_rate = true_pos / nrow(pums_100_final); #calculate the true positive rate
  history[row_counter, ] <- c(d, tp_rate, output$criticalVal);
  row_counter = row_counter + 1;
}
print(Sys.time())



#### Plot results ####
f_size = 15;
history <- as.data.frame(history);
colnames(history) <- c("num_attr", "tpr", "crit_val");
# Plot average RMSE of reconstruction against noise input
p <- ggplot(data = history, aes(x=history$num_attr, y=history$tpr)) + geom_point();
p <- p + labs(title = paste("Membership attack with", noise_type, "noise =", noise_val), x="Number of attributes", y = "True positive rate") + theme(plot.title = element_text(hjust=0.5), text = element_text(size=f_size));

#### Export the graph
ggsave(filename = paste("./figs/memAttack", noise_type, "noise.jpg", sep = "_"), plot=p, width = 10, height = 6);
#### Save data so I don't have to run the simualtion again if I want to re-plot the data
write.csv(history, paste("./memAttack", noise_type, "noise.csv", sep = "_"));


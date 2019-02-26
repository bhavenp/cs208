##
##  HW1_2.r
##
##  Implement a reconstruction attack on the sample of 100 points from the PUMS dataset
##
##  BP 2.18.2019
##

rm(list=ls());
library(plyr); #import this library for doing rounding for noise
library(ggplot2); #import library for plotting
library(grid);

#### Parameters ####
prime <- 113; # prime number for hashing creating random vectors to pass to query
n <- 100;        # Dataset size
k.trials <- 2*n;  # Number of queries
num_exps <- 10; #number of experiments
noise_input <- "Rounding"; # What type of noise will be used as defense. Can be "Rounding", "Gaussian", or "Subsampling"

noise_vec <- c(1:100); #noise parameters for Rounding and Subsampling
# noise_vec <- c(seq(1, 1.9, 0.1), 2:100); #noise parameters for Gaussian


#### Import Data ####
pums_100 <- read.csv(file="../../data/FultonPUMS5sample100.csv"); #read in data from data folder
#trim dataset to only the PUB categories plus "uscitizen"
pums_100_trimmed <- pums_100[, c("sex", "age", "educ","income", "latino", "black", "asian",  "married", "divorced", "children", "employed", "militaryservice", "disability", "englishability", "uscitizen")];

#### Query function ####
query <- function(random_vec, data, prime, noise_type, noise_param){
  data_mat <- data.matrix(data[, 1:14]); # convert all rows and PUB categories into matrix
  #multiply matrix by random vector and 
  person_in_sum <- (data_mat %*% random_vec) %% prime %% 2; 
  us_citz_sum <- sum(person_in_sum * data[, 15]);
  
  if(noise_type == "Rounding"){
    noisy_sum <- round_any(us_citz_sum, noise_param); #this function will round the sum according to a multiple of noise_param 
  }
  else if(noise_type == "Gaussian"){
    noisy_sum <- us_citz_sum + rnorm(n=1, mean=0, sd=noise_param); #add noise sampled by from N(0, noise_param^2)
  }
  else{
    subsamp_ind <- sample(x=1:length(person_in_sum), size=noise_param, replace=FALSE); #get a random sample of the indices w/o replacement
    subsamp_sum <- sum(person_in_sum[subsamp_ind] * data[subsamp_ind, 15]); #get the US citizen count for the subsample
    noisy_sum <- subsamp_sum * (n / noise_param); #sum to report will be subsamp_sum x (scaling factor)
  }
  
  # return the actual sum, the noisy sum and indices. I am returning the indices here so that I don't have to recalculate which individuals were included in the sum.
  return(list(us_citz_sum=us_citz_sum, noisy_sum=noisy_sum, indices=person_in_sum));
}

#### Run experiment function ####
#### Give a prime number used for querying, data frame of PUB values, and how noise should be added to query
run_experiment <- function(prime, data_input, noise_type, noise_to_add){
  #### Here we run our query repeatedly and record results
  history <- matrix(NA, nrow=k.trials, ncol=n+2);  # a matrix to store results in
  
  for(i in 1:k.trials){
    rand_vec <- sample(0:prime-1, size = ncol(data_input)-1);
    res <- query(random_vec=rand_vec, data=data_input, prime=prime, noise_type=noise_type, noise_param=noise_to_add);
    history[i,] <- c(res$us_citz_sum, res$noisy_sum, res$indices);  # save into our history matrix
  }
  
  #### Convert matrix into data frame
  xnames <- paste("x", 1:n, sep="");
  varnames<- c("y", xnames);
  # convert noisy sum and indices in matrix into data frame
  releaseData <- as.data.frame(history[, 2:ncol(history) ]);
  names(releaseData) <- varnames; #add column names to data frame
  
  #### Run a linear regression
  formula <- paste(xnames, collapse=" + ");    # construct the formula, y ~ x1 ... xn -1
  formula <- paste("y ~ ", formula, "-1");
  formula <- as.formula(formula);
  # print(formula);
  
  output <- lm(formula, data=releaseData);                   # run the regression
  estimates <- output$coef;                                  # save the estimates                   
  estimate_conv <- (estimates>0.5); # convert estimates to binary values
  
  sensitiveData <- data_input[, "uscitizen"];
  exp_acc <- sum(estimate_conv == sensitiveData) / n; #calculate the fraction of USCITIZEN correctly reconstructed for this experiment
  #calculate RMSE between the exact value of our query and the noisy answer we passed back
  exp_rmse <- ( sum((history[, 2] - history[, 1]) ** 2) / nrow(history)) ** 0.5; 
  
  #return the 
  return(list(exp_rmse=exp_rmse, exp_acc=exp_acc))
}

print(Sys.time())
#### Run through all noise parameters, each with 10 experiments ####
final_results <- matrix(NA, nrow=length(noise_vec), ncol=3);  # a matrix to store results in

for(i in 1:length(noise_vec)){
  noise_to_add <- noise_vec[i];
  agg_rmse <- c(); #empty vector to hold the RMSE values of individual experiments
  agg_acc <- c(); #empty vector to hold the accuracy values of individual experiments
  #need to go through num_exps for each noise parameter
  for(e in 1:num_exps){
    exp_res = run_experiment(prime = prime, data_input=pums_100_trimmed, noise_type=noise_input, noise_to_add=noise_to_add);
    agg_rmse <- c(agg_rmse, exp_res$exp_rmse); #add RMSE to vector
    agg_acc <- c(agg_acc, exp_res$exp_acc); #add acc to vector
  }
  #put average of RMSEs and Accs from the experiments into the matrix
  final_results[i, ] <- c(noise_to_add, mean(agg_rmse), mean(agg_acc));
}
print(Sys.time())

final_results <- as.data.frame(final_results);
colnames(final_results) <- c("Param_vals", "RMSE", "Acc")
#### Plot results ####
f_size = 15;
# Plot average RMSE of reconstruction against noise input
p_rmse <- ggplot(data = final_results, aes(x=final_results$Param_vals, y=final_results$RMSE)) + geom_point();
p_rmse <- p_rmse + labs(x=paste(noise_input, "noise"), y = "Average RMSE") + theme(plot.title = element_text(hjust=0.5), text = element_text(size=f_size));
# Plot average accuracy of reconstruction against noise input
p_acc <- ggplot(data = final_results, aes(x=final_results$Param_vals, y=final_results$Acc)) + geom_point(); 
p_acc <- p_acc + geom_hline(yintercept = 0.96, linetype="dashed", color = "blue"); 
p_acc <- p_acc + labs(x=paste(noise_input, "noise"), y = "Average Accuracy") + theme(plot.title = element_text(hjust=0.5), text = element_text(size=f_size));
# Plot average RMSE vs average accuracy of reconstruction
p_rmse_acc <- ggplot(data = final_results, aes(x=final_results$RMSE, y=final_results$Acc)) + geom_point(); 
p_rmse_acc <- p_rmse_acc + geom_hline(yintercept = 0.96, linetype="dashed", color = "blue");
p_rmse_acc <- p_rmse_acc + labs(x="Average RMSE", y = "Average Accuracy") + theme(plot.title = element_text(hjust=0.5), text = element_text(size=f_size));

#create grid for plotting
gs <- grid.arrange(p_rmse, p_acc, p_rmse_acc, nrow=1, ncol=3, top=textGrob(paste("Average RMSE & Accuracy for",noise_input,"noise"), gp=gpar(fontsize=15)) );

#### Export the graph
ggsave(filename = paste("./figs/regAttack", noise_input, "noise.jpg", sep = "_"), plot=gs, width = 11, height = 6);

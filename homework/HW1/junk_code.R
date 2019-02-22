query1 <- function(random_vec, data, prime){
  us_citz_sum = 0;
  person_in_sum <- c(); #vector of 0s and 1s depending on if person was added to dataset
  #
  for(i in 1:nrow(data)){ # go through each row to hash it
    dot_prod <- random_vec * data[i, 1:14]; #get dot product of random vector and row from data
    in_subset <- sum(dot_prod) %% prime %% 2;
    if(in_subset){ #check if this person should be included in the subset
      us_citz_sum = us_citz_sum + data[i, 15];
      person_in_sum <- c(person_in_sum, 1);
    }
    else{
      person_in_sum <- c(person_in_sum, 0);
    }
  }
  
  #need to add noise
  
  return(us_citz_sum);
}
res1 <- query1(random_vec=rand_vec, data=pums_100_trimmed, prime=prime);
# Given v, getting k using the Algorithm 2 in the paper

source("Functions/GetCriticalValues/get_k_two_step.R")

## Define a function to get the Type-B v in the paper
get_v_B <- function(v_max){
  i <- 5
  v_vec <- ceiling((1:i)^2/2)
  while (tail(v_vec,1) < v_max) {
    i <- i+1
    v_vec <- ceiling((1:i)^2/2)
  }
  return(ceiling((1:(i-1))^2/2))
}

###
alpha <- 0.05
v_max = 65

v <- get_v_B(v_max)

##########
num_trial <- 200
B <- 100000
Bernoulli_mat <- matrix(sample(c(-1,1), B*num_trial, replace=TRUE),nrow=B, ncol=num_trial)

res <- get_k_two_step(v, Bernoulli_mat, alpha, step1_size_vec=c(1,0.5,0.1,0.05,0.01), 
               step2_size_vec=c(50,25,10,5,1))
k_twostep <- res$k_list[[3]]

save(k_twostep, file = "ExampleOfGettingVK/k_twostep.RData")
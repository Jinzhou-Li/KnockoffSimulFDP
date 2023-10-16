# Given v, getting critical k_ct for closed testing (cf. Section 4.2 and 5.2 in the paper)
library(parallel)
source("Functions/GetCriticalValues/get_criticalValue_ct.R")

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
p <- 200

v <- get_v_B(v_max)

##########
num_trial <- 200
# We used a small B below for illustration purpose. 
# One should use a large B, e.g., B=100000, in order to obtain an accurate k_mat_ct
B <- 1000    
Bernoulli_mat <- matrix(sample(c(-1,1), B*num_trial, replace=TRUE),nrow=B, ncol=num_trial)

# load k_twostep and use it as an initial vector for faster computation (optional).
load("ExampleOfGettingVK/k_twostep.RData")

#####
get_Critical_Weight1_simu <- function(I_length, v, alpha, Bernoulli_mat, k_twostep){
  if(I_length==1) return(rep(2,length(v)))
  Bernoulli_mat_temp <- Bernoulli_mat[,1:I_length]
  return( get_criticalValue_ct(v, Bernoulli_mat_temp, alpha, weight_type="1", 
                               step_size_vec=c(100,50,25,10,5,1), initial_vec=k_twostep)[[1]] )
}

k_mat_ct <- t(mcmapply(get_Critical_Weight1_simu, 1:p, 
                           MoreArgs=list(v, alpha, Bernoulli_mat, k_twostep), 
                           mc.cores=10) )
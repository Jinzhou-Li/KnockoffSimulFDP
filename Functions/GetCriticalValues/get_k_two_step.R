# Algorithm 2 in the paper: 
# get k such that the approximated/simulated P( NB^p(v_1)<=k_1-1, ..., NB^p(v_m)<=k_m-1 ) >= 1-alpha

# Bernoulli_mat' as input to make sure different types of tuning parameters are based on the same 
# simulated Bernoulli trials (B x p)

get_k_two_step <- function(v_vec, Bernoulli_mat, alpha, step1_size_vec=c(1,0.5,0.1,0.05,0.01), 
                           step2_size_vec=c(50,25,10,5,1)){
  if(length(v_vec)<2){stop("when length(v_vec)=1, can calculate k based on JS, no need to use this function!")}
  c_alpha <- -log(alpha)/log(2-alpha) 
  
  # based on the simulated Bernoulli trials, get 'stat_matrix' using 'v_vec' (B x m)
  Test_stat_mat <- mapply(function(i, Bernoulli_mat_input, v_vec_input){
    return(apply(Bernoulli_mat_input, 1, get_NBstat_simu, v_vec_input[i]))},
    1:length(v_vec), MoreArgs=list(Bernoulli_mat, v_vec))
  
  # find the smallest 'j' s.t. 'j-c_j+1=max(v_vec)': to obtain 'k' for 'max(v_vec)'
  j_max <- max(v_vec)
  c_temp <- floor(c_alpha*(1+j_max)/(1+c_alpha)) + 1
  v_max_temp <- j_max - c_temp + 1
  while (v_max_temp < max(v_vec)) {
    j_max <- j_max+1
    c_temp <- floor(c_alpha*(1+j_max)/(1+c_alpha)) + 1
    v_max_temp <- j_max - c_temp + 1
  }
  
  # Step0: get 'k^{raw}' and 'j^*' based on 'v_vec' (see definition of 'k^{raw}' in the paper)
  c_raw <- floor(c_alpha*(1+(1:j_max))/(1+c_alpha)) + 1
  b_raw <- (1:j_max) - c_raw + 1
  j_star <- mapply(function(i, b_input){
    return(which(b_input==i)[1])}, v_vec, MoreArgs=list(b_raw))
  k_raw <- pmin(c_raw[j_star], ncol(Bernoulli_mat)+1) # k-1 must be less than number of trials
  prob_raw <- 1 - mean(apply(t(t(Test_stat_mat) >= k_raw), 1, sum) > 0)  # original probability: didn't exhaust the alpha level.
  
  ### Two-step approach
  ### Step 1
  Step1_list <- get_k_updateC(Test_stat_mat, j_star, alpha, step1_size_vec)
  k_step1 <- pmin(Step1_list[[1]], ncol(Bernoulli_mat)+1)
  prob_step1 <- Step1_list[[2]]   # take a look: more or less exhaust the alpha level.
  c_step1 <- Step1_list[[3]]      # compare to the previous c(alpha)
  
  ### Step 2
  Step2_list <- get_k_greedy(Test_stat_mat, k_step1, v_vec, alpha, step2_size_vec)
  k_step2 <- pmin(Step2_list[[1]], ncol(Bernoulli_mat)+1)
  prob_step2 <- Step2_list[[2]]   # take a look: more or less exhaust the alpha level.
  
  return(list(k_list=list(k_raw, k_step1, k_step2), 
              #b_list=list(b_raw, b_step1, b_step2), 
              c_step1,
              prob_list=list(prob_raw, prob_step1, prob_step2)))
}

############################### Other functions
##### get NB test statistic based on one sequence of Bernoulli trials
get_NBstat_simu <- function(vec, para){
  if(sum(vec<0) >= para){return( sum(vec[1: which(vec<0)[para] ] > 0) )}
  if(sum(vec<0) < para){return( sum(vec>0) )}
}

##### step 1 update
get_k_updateC <- function(Test_stat_mat, j_star, alpha, step_size_vec){
  c_alpha <- -log(alpha)/log(2-alpha) 
  
  for (i in 1:length(step_size_vec)) {
    step_size <- step_size_vec[i]
  
    k_temp <- floor(c_alpha*(1+j_star)/(1+c_alpha)) + 1
    prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_temp), 1, sum) > 0)  
    
    c_alpha_temp <- c_alpha
    while (prob>=1-alpha) {
      c_alpha <- c_alpha_temp
      c_alpha_temp <- c_alpha_temp - step_size
      
      k_temp <- floor(c_alpha_temp*(1+j_star)/(1+c_alpha_temp)) + 1
      prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_temp), 1, sum) > 0)   
    }
  }
  
  k_step1 <- floor(c_alpha*(1+j_star)/(1+c_alpha)) + 1
  prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_step1), 1, sum) > 0)
  return(list(k_step1, prob, c_alpha))
}

######## step 2 update
get_k_greedy <- function(Test_stat_mat, k_vec, v_vec, alpha, step_size_vec){
  for (i in 1:length(k_vec)) {
    for (j in 1:length(step_size_vec)) {
      step_size <- step_size_vec[j]
      
      if(k_vec[i]-step_size < v_vec[i] ){next} # if step size is too large
      
      Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= k_vec), 1, sum) > 0)
      
      k_greedy_temp <- k_vec
      while (Prob_temp >= 1-alpha) {
        k_vec <- k_greedy_temp
        k_greedy_temp[i] <- k_greedy_temp[i] - step_size
        Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= k_greedy_temp), 1, sum) > 0)
      }
    }
  }
  
  prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_vec), 1, sum) > 0)
  return(list(k_vec, prob))
}
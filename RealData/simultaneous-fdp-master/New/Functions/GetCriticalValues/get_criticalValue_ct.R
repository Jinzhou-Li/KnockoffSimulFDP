# Get the critical values for the closed testing (weight=1 or local rank) via greedy approach

# For 'weight_type=="1"': 'para_vec' is 'v_vec'
# (want to compare to the interpolation version, so we take the same critical values 'k^{p}' as initial value)
# For 'weight_type=="LocalRank"': 'para_vec' is 'b_vec' 
# (want to compare to weight=1 version, so use the same 'Bernoulli_mat', and 'b' is from weight=1 version (b=v+k-1))

get_criticalValue_ct <- function(para_vec, Bernoulli_mat, alpha, weight_type, 
                                 step_size_vec=c(200,100,50,25,10,5,1), initial_vec=NULL){
  
  # based on the simulated Bernoulli trials, get 'stat_matrix' using 'b_vec' (B x m)
  if(weight_type=="1"){
    Test_stat_mat <- mapply(function(i, Bernoulli_mat_input, v_vec_input){
      return(apply(Bernoulli_mat_input, 1, get_NBstat_simu, v_vec_input[i]))},
      1:length(para_vec), MoreArgs=list(Bernoulli_mat, para_vec))
    if(is.null(initial_vec)){stop("Need to give 'initial_vec' for weight=1 version")}
  }
  if(weight_type=="LocalRank"){
    Test_stat_mat <- mapply(function(i, Bernoulli_mat_input, b_vec_input){
      return(apply(Bernoulli_mat_input, 1, get_RankStat_simu, b_vec_input[i]))},
      1:length(para_vec), MoreArgs=list(Bernoulli_mat, para_vec))
    initial_vec <- apply(Test_stat_mat, 2, max) + 1 
  }
  
  ### Greedily approximate critical values
  greedy_result_list <- get_Critical_greedy(Test_stat_mat, alpha, initial_vec, step_size_vec)
  critical_vec <- greedy_result_list[[1]]
  prob <- greedy_result_list[[2]]   # take a look: more or less exhaust the alpha level.
  
  return(list(critical_vec, prob))
}

############################### Other functions
##### get test statistic (with local rank as weight) based on one sequence of Bernoulli trials
get_RankStat_simu <- function(vec, b){
  I_size <- length(vec)
  sign_rank_vec <- (I_size:1)*vec[1:I_size]
  return(sum((sign_rank_vec[1:min(b,I_size)] > 0) * sign_rank_vec[1:min(b,I_size)]))
}

##### get test statistic (with weight 1) based on one sequence of Bernoulli trials
get_NBstat_simu <- function(vec, para){
  if(sum(vec<0) >= para){return( sum(vec[1: which(vec<0)[para] ] > 0) )}
  if(sum(vec<0) < para){return( sum(vec>0) )}
}

######## greedy approach to get critical values
get_Critical_greedy <- function(Test_stat_mat, alpha, initial_vec, step_size_vec){
  critical_vec <- initial_vec
  
  # greedily update critival values one by one
  for (i in 1:length(critical_vec)) {
    for (j in 1:length(step_size_vec)) {
      step_size <- step_size_vec[j]
      if(critical_vec[i]-step_size < 0 ){next} # if step size is too large
      
      Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= critical_vec), 1, sum) > 0)
      
      critical_vec_temp <- critical_vec
      while (Prob_temp >= 1-alpha) {
        critical_vec <- critical_vec_temp
        critical_vec_temp[i] <- critical_vec_temp[i] - step_size
        Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= critical_vec_temp), 1, sum) > 0)
      }
    }
  }
  prob <- 1 - mean(apply(t(t(Test_stat_mat) >= critical_vec), 1, sum) > 0)
  return(list(critical_vec, prob))
}
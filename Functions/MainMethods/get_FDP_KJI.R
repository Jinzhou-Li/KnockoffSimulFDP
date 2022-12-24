# Algorithm 1 in the paper: get simultaneous FDP based on joint-KFWER and interpolation

# inputs: W: satisfying |W_1| > ... > |W_m| > 0 
#         v, k satisfy the inequality (14) in the paper

# output: a vector of FDP bound for sets 'R: which(W >= abs(W[i]))'

get_FDP_KJI <- function(W, k_vec, v_vec){
  if(sum(W==0)>0 | length(unique(W))!=length(W) | !identical(W[order(abs(W), decreasing=TRUE)],W)){
    stop("The input W might have ties or zeros or not sorted!")}
  m <- length(k_vec)
  S_list <- lapply(1:m, function(i_S, W, v_vec){
    candidate_set <- which(cumsum(W<0)==v_vec[i_S])
    ifelse(length(candidate_set)==0, threshold <- min(abs(W)), threshold <- abs(W[candidate_set[1]]))
    S <- which(W>=threshold)
    return(S)
    }, W, v_vec)
  
  ###
  p <- length(W)
  FDP_bound_vec <- rep(1,p)
  for (i in 1:p) {
    R <- which(W >= abs(W[i]))
    
    if(length(R)==0){FDP_bound_vec[i] <- 0; next}  # this means all W are negative, so discovery set is empty, so FDP is 0.
    
    FDP_k_temp <- mapply(function(j, p, R, S_list, k_vec){
      S_temp <- S_list[[j]]
      return( min(length(R), k_vec[j]-1+length(setdiff(R,S_temp))) / max(1,length(R)) )
      }, 1:m, MoreArgs=list(p, R, S_list, k_vec))
    
    FDP_bound_vec[i] <- min(FDP_k_temp)
  }
  return(FDP_bound_vec)
}
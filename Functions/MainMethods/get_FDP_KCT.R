##### Get FDP bounds for sets 'which(W >= abs(W[i]))' using closed testing, i=1,...,p

get_FDP_KCT <- function(W, b_mat, critical_value_mat, weight_type){
  p <- length(W)
  
  FDP_bounds <- rep(1,p)
  # Set the FDP bound as 0 when no discovery: e.g., when all front $W_i's$ are negative
  if(which(W>0)[1] != 1) {FDP_bounds[1:(which(W>0)[1]-1)] <- 0}

  Sets_list <- lapply(1:p, function(thre_index,W){return(which(W >= abs(W[thre_index])))}, W)
  
  ######### Two initial steps to possibly reduce computation
  ### 1. locally test all H_{I} for I in 'Sets_list' 
  k_1 <- 1
  for(i_1 in p:1){
    I = Sets_list[[i_1]]
    local_test_I <- local_test_fun(W, I, b_mat[length(I),], critical_value_mat[length(I),], weight_type)
    if(local_test_I==0){k_1 <- i_1; break}
  }
  if(k_1==p) {return(FDP_bounds)} # If H_{[p]} is not locally rejected, all FDP bounds are 1.
  
  ### 2. Check whether H_{[k]}, k = p, ..., k_1, is rejected by closed testing
  k_2 <- k_1
  I <- Sets_list[[p]]
  closed_test_I <- closed_testing_fun(W, I, b_mat, critical_value_mat, weight_type)
  if(closed_test_I==0){
    return(FDP_bounds)  # if so, any 'Sets_list[[k]]' will not be rejected by closed testing
  }else{
    for (i_2 in (p-1):k_1) {
      if(identical(Sets_list[[i_2]], Sets_list[[i_2+1]])){next}
      I <- Sets_list[[i_2]]
      closed_test_I <- closed_testing_fun(W, I, b_mat, critical_value_mat, weight_type)
      # if some 'Sets_list[[k_2]]' is not rejected by closed testing, then all 'Sets_list[[i]]' with 'i < k_2' will not, so can stop
      if(closed_test_I==0){k_2 <- i_2; break} 
    }
  } 
  
  ### find FDP bounds for Sets_list[[k]], k=k_2+1,...,p 
  max_size_not_rej <- length(Sets_list[[k_2]]) # record the maximal size of set that is not rejected by closed testing
  for (i in (k_2+1):p) {
    # here it is fine to use 'i-1' as k_2>=4
    if(identical(Sets_list[[i]], Sets_list[[i-1]])){FDP_bounds[i] <- FDP_bounds[i-1]; next} # due to the form of R's
    I <- Sets_list[[i]]
    
    # (already tested above that I itself will not be rejected by closed testing)
    # Find the largest subset of I that is not rejected by closed testing
    # Due to the form of R's (nested), must exist a subset of size 'max_size_not_rej' that is not rejected, so can stop.
    for (t in (length(I)-1):max_size_not_rej) {  
      
      # only for checking, delete later
      if(length(I) <= max_size_not_rej) {print("length(I) <= max_size_not_rej")}
      
      B <- get_B_subset(t, W, I)
      closed_test_B <- closed_testing_fun(W, B, b_mat, critical_value_mat, weight_type)
      if(closed_test_B==0) {
        FDP_bounds[i] <- t/length(I)
        if(t>max_size_not_rej) {max_size_not_rej <- t} #; print(paste0("max_size_not_rej: ",max_size_not_rej))}
        break
      }
    }
  } # for-loop of i
  return(FDP_bounds)
}

#### Get set B: the set we look at for given size t \in [k-1,1] in (B)  
get_B_subset <- function(t, W, I){
  I_length <- length(I)
  W_I <- W[I]
  D_I <- sign(W_I)==1
  r_I <- rank(abs(W_I))
  B <- I[order(D_I*r_I)[1:t]]
  return(B)
}

##### Check whether H_I is locally rejected
local_test_fun <- function(W, I, b_vec, critical_value_vec, weight_type){
  test_stat_vec <- get_test_stat(W, I, b_vec, weight_type)
  return(max(test_stat_vec >= critical_value_vec)) # '1' is reject and '0' is not reject
}

##### Check whether H_I is rejected by closed testing
closed_testing_fun <- function(W, I, b_mat, critical_value_mat, weight_type){
  p <- length(W)
  k <- length(I)
  
  # Check whether all supersets of I are locally rejected from size k,...,p
  # First locally test I (superset of size k)
  local_test_I <- local_test_fun(W, I, b_mat[k,], critical_value_mat[k,], weight_type)
  if(local_test_I==0){return(0)}
  
  # If H_I is rejected: check for supersets with size k+1,...,p
  if(k==p){return(1)} # in this case, we already know from above that it is not rejected
  # otherwise if k<p
  for (t in 1:(p-k)) {
    S <- get_S_superset(t, W, I)
    I_super <- c(I,S)
    local_test_I_super <- local_test_fun(W, I_super, b_mat[k+t,], critical_value_mat[k+t,], weight_type)
    if(local_test_I_super==0) {return(0)}
  }
  return(1) # If rejected all supersets, we return 1
}

#### Get set S: the set we look at for given size t \in [k+1,p] in (1)
get_S_superset <- function(t, W, I){
  p <- length(W)
  I_c <- (1:p)[-I]
  W_I_c <- W[I_c]
  D_I_c <- sign(W_I_c)==1
  r_I_c <- rank(abs(W_I_c))
  I_c_length <- length(I_c)
  
  S <- I_c[order(D_I_c*r_I_c)[1:t]]
  return(S)
}

########### Function to calculate the test statistics
get_test_stat <- function(W, I, b_vec, weight_type){
  W_I <- W[I]
  I_length <- length(I)
  D_I <- sign(W_I)==1
  r_I <- rank(abs(W_I))
  
  test_stat_vec <- mapply(get_test_stat_sub, 1:length(b_vec), 
                          MoreArgs=list(b_vec, r_I, D_I, weight_type))
  return(test_stat_vec)
}

### 
get_test_stat_sub <- function(b_index, b_vec, r_I, D_I, weight_type){
  b <- b_vec[b_index]
  I_length <- length(r_I)
  if(weight_type=="1"){weight <- (r_I > I_length-b)}
  if(weight_type=="LocalRank"){weight <- r_I * (r_I > I_length-b)}
  return(sum(weight * D_I))
}
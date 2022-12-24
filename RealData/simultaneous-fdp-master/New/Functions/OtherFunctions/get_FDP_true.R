# Get true FDPs 

# Inputs: W: satisfying |W_1| > ... > |W_p| > 0
#         true_index: the vector indicating which elements in W are true
# Output: true FDPs: the i-th element is the FDP for set R_i: 'which(W >= abs(W[i]))'

get_FDP_true <- function(W, true_index){
  if(sum(W==0)>0 | length(unique(W))!=length(W) | !identical(W[order(abs(W), decreasing=TRUE)],W)){
    stop("The input W might have ties or zeros or not sorted!")}
  
  p <- length(W)
  true_num_disc <- rep(0,p)
  for (i in 1:p) {
    count_temp <- 0
    for (j in 1:i) {
      if(j %in% true_index & W[j]>0) {count_temp <- count_temp+1}
    }
    true_num_disc[i] <- count_temp
  }
  
  num_disc <- cumsum(W>0)
  FDP_true <- (num_disc-true_num_disc)/pmax(num_disc,1)
  return(FDP_true)
}
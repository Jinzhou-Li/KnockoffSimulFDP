# Method from Katsevich and Ramdas (2020)
# Compute FDP-bar based on pre-processed knockoff statistics (no tie, no zeros, sorted)
# Codes below are based on the codes in the file 'UKBB_utils.R' from Github "https://github.com/ekatsevi/simultaneous-fdp"

# Input: W: satisfying |W_1| > ... > |W_p| > 0 
#        alpha: nominal level

# Output: the FDP bounds. The i-th FDP bound is for set R_i: 'which(W >= abs(W[i]))'

get_FDP_bar_KR = function(W, alpha){
  
  if(sum(W==0)>0 | length(unique(W))!=length(W) | !identical(W[order(abs(W), decreasing=TRUE)],W)){
    stop("The input W might have ties or zeros or not sorted!")}
  
  C = -log(alpha)/log(2-alpha)   # compute knockoffs FDP multipler
  
  num_discoveries = cumsum(W>0) 
  # initial bound on V
  max_false_positives = pmin(num_discoveries, floor(C*(1+cumsum(W<0))))
  min_true_positives = num_discoveries - max_false_positives
  
  # interpolation
  max_false_positives = rev(cummin(rev(max_false_positives)))
  min_true_positives = cummax(min_true_positives)
  max_false_positives = pmin(max_false_positives, num_discoveries - min_true_positives)
  
  max_false_positives <- pmin(max_false_positives, num_discoveries)
  FDP_bound_kfwer_interp <- max_false_positives/pmax(1,num_discoveries)
  
  return(FDP_bound_kfwer_interp)
}
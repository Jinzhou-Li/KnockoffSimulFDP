# The pre-process function for W:
# given the original W, output the sorted W without ties and zeros (|W_sort[1]|>...>|W_sort[p]|>0)
# and the original position of W_sort[1], ..., W_sort[p] in the original W (i.e. which variable is W_sort[i])

# Input: original W
# Output: W_sort, W_order

preprocess_W_func <- function(W){    
  W_order <- order(abs(W), decreasing=TRUE)
  W_sort <- W[W_order]
  # delete zeros
  zero_index <- which(W_sort==0)
  if(length(zero_index)>0){
    W_sort <- W_sort[-zero_index]
    W_order <- W_order[-zero_index]
  }
  
  # break the tie if there is any without changing the order
  W_sort_abs <- abs(W_sort)
  for (i in 1:length(W_sort_abs)) {
    temp <- W_sort_abs[i]
    if(sum(W_sort_abs==temp) >= 2){     # if there is tie
      tie_index <- which(W_sort_abs==temp)
      first_index <- tie_index[1]
      last_index <- tail(tie_index,1)
      
      if(last_index != length(W_sort_abs)){
        max_value <- W_sort_abs[last_index] - W_sort_abs[last_index+1]
        W_sort_abs[tie_index] <- W_sort_abs[tie_index]- seq(0, max_value/2,length.out=length(tie_index))
      }else{ # if the tie contains the last element
        max_value <- W_sort_abs[first_index-1] - W_sort_abs[first_index]
        W_sort_abs[tie_index] <- W_sort_abs[tie_index] + seq(max_value/2,0,length.out=length(tie_index))
      }
    } # end if
  } # end for
  
  W_sort <- sign(W_sort) * W_sort_abs
  return(list(W_sort=W_sort,W_order=W_order))
}

### example
# W <- rnorm(15)
# W[c(3,9)] <- 0
# W[6] <- W[8] <- W[1]
# W[11] <- W[2]
# W
# preprocess_W_func(W)

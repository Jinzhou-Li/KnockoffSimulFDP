### Simulation example 1: generating W through a linear model, comparing KR and KJI
library(knockoff)

source("Functions/OtherFunctions/preprocess_W_func.R")
source("Functions/OtherFunctions/get_FDP_true.R")
source("Functions/OtherFunctions/get_FDP_bar_KR.R")
source("Functions/MainMethods/get_FDP_KJI.R")

######### Simulation setting:
set.seed(2022)
p = 200          # number of variables
n = 500         # number of observations
alpha = 0.05    

s <- 0.2
k = ceiling(s*p)    # number of variables with nonzero coefficients
amplitude = 6      # signal amplitude
rho = 0.6
Sigma = toeplitz(rho^(0:(p-1)))
nonzero = sample(p, k)  
beta = amplitude * (1:p %in% nonzero) / sqrt(n)

######### Generate data
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = X %*% beta + rnorm(n)

######### get original W using knockoff method (one can use different ways to get W. Below is just one option for illustration purpose)
W_ori <- knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax, fdr=0.1)$statistic
W_list <- preprocess_W_func(W_ori)  # Get the preprocessed W: |W_1|>...>|W_m|>0
W <- W_list[[1]]
W_order <- W_list[[2]]

#### get the true FDP
true_index <- mapply(function(i, W_order ){return(which(W_order==i))}, which(beta!=0), MoreArgs=list(W_order) )
FDP_true <- get_FDP_true(W, true_index)

########## Implement methods
FDP_KR <- get_FDP_bar_KR(W, alpha)

# Load the calculated tuning parameters v and k
load(paste0("Examples/TuningParameters/CriticalValuesVmax65.RData"))
FDP_JKI_A <- get_FDP_KJI(W, k_vec=k_list[[1]], v_vec=v_list[[1]])
FDP_JKI_B <- get_FDP_KJI(W, k_vec=k_list[[2]], v_vec=v_list[[2]])
FDP_JKI_C <- get_FDP_KJI(W, k_vec=k_list[[3]], v_vec=v_list[[3]])
FDP_JKI_D <- get_FDP_KJI(W, k_vec=k_list[[4]], v_vec=v_list[[4]])

### plot
plot(FDP_true, cex=0.1, ylim=c(0,1))
lines(FDP_KR, col="orange")
lines(FDP_JKI_A, col="skyblue1")
lines(FDP_JKI_B, col="blue")
lines(FDP_JKI_C, col="darkviolet")
lines(FDP_JKI_D, col="green")

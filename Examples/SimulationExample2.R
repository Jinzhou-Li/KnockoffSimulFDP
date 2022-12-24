### Simulation example 2: Directly generating W, comparing KJI and KCT
source("Functions/OtherFunctions/get_FDP_true.R")
source("Functions/MainMethods/get_FDP_KJI.R")
source("Functions/MainMethods/get_FDP_KCT.R")

######### Simulation setting:
set.seed(2022)
p = 50          # number of variables
alpha = 0.05        

##### Generate W_sort and the related FDP
W_sign <- rep(0,p)

### false null variables
position_false <- 10:20
position_true <- (1:p)[-position_false]
Berboulli_true <- rbinom(length(position_true), 1, 1)
Berboulli_true[which(Berboulli_true==0)] <- -1
W_sign[position_true] <- Berboulli_true

### true null variable
Berboulli_false <- rbinom(length(position_false), 1, 0.5)
Berboulli_false[which(Berboulli_false==0)] <- -1
W_sign[-position_true] <- Berboulli_false

W <- W_sign * p:1
FDP_true <- get_FDP_true(W, true_index=position_true)

########## Implement methods
# Load the calculated tuning parameters v and k
load(paste0("Examples/TuningParameters/CriticalValuesVmax50.RData"))

FDP_JKI_A <- get_FDP_KJI(W, k_vec=k_list[[1]], v_vec=v_list[[1]])
FDP_JKI_B <- get_FDP_KJI(W, k_vec=k_list[[2]], v_vec=v_list[[2]])
FDP_JKI_C <- get_FDP_KJI(W, k_vec=k_list[[3]], v_vec=v_list[[3]])
FDP_JKI_D <- get_FDP_KJI(W, k_vec=k_list[[4]], v_vec=v_list[[4]])

# get 'b_mat'(p x m) for CT: each row corresponds to 'b^{|I|}'
b_mat_A <- t(v_list[[1]] + t(CT_weight1_list[[1]]) - 1)  
FDP_KCT_A <- get_FDP_KCT(W, b_mat_A, CT_weight1_list[[1]], weight_type="1")

b_mat_B <- t(v_list[[2]] + t(CT_weight1_list[[2]]) - 1)  
FDP_KCT_B <- get_FDP_KCT(W, b_mat_B, CT_weight1_list[[2]], weight_type="1")

b_mat_C <- t(v_list[[3]] + t(CT_weight1_list[[3]]) - 1)  
FDP_KCT_C <- get_FDP_KCT(W, b_mat_C, CT_weight1_list[[3]], weight_type="1")

b_mat_D <- t(v_list[[4]] + t(CT_weight1_list[[4]]) - 1)  
FDP_KCT_D <- get_FDP_KCT(W, b_mat_D, CT_weight1_list[[4]], weight_type="1")

### plots
plot(FDP_true, cex=0.1, ylim=c(0,1))
lines(FDP_JKI_A, col="skyblue1")
lines(FDP_KCT_A, col="skyblue1", lty="dashed")
lines(FDP_JKI_B, col="blue")
lines(FDP_KCT_B, col="blue", lty="dashed")
lines(FDP_JKI_C, col="darkviolet")
lines(FDP_KCT_C, col="darkviolet", lty="dashed")
lines(FDP_JKI_D, col="green")
lines(FDP_KCT_D, col="green", lty="dashed")
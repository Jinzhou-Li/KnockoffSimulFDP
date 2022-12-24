######## The code is largely based on https://github.com/ekatsevi/simultaneous-fdp

library(knockoff)
library(ggpubr)
# set up workspace

setwd("RealData/simultaneous-fdp-master")
source("setup.R")

# analysis parameters
alpha = 0.05              # confidence level for simultaneous FDP bound
q = 0.1                   # FDR target level
phenotype = "platelet"    # phenotype of interest

# download KnockoffZoom data, if necessary
download_KZ_data(data_dir)

# read in processed data
df = read_KZ_data(data_dir, phenotype)

# compute knockoffs FDP multipler
C = -log(alpha)/log(2-alpha)

# compute FDP-hat and FDP-bar
FDP_hat_vec = get_FDP_hat(df$W)
FDP_bar_vec = get_FDP_bar(df$W, C)

######### 
source("New/Functions/MainMethods/get_FDP_KJI.R")
load(paste0("New/TuningParameters/CriticalValuesVmax1200.RData"))

### run KJI
FDP_KJI_A <- get_FDP_KJI(df$W, k_vec=k_list[[1]], v_vec=v_list[[1]])
FDP_KJI_B <- get_FDP_KJI(df$W, k_vec=k_list[[2]], v_vec=v_list[[2]])
FDP_KJI_C <- get_FDP_KJI(df$W, k_vec=k_list[[3]], v_vec=v_list[[3]])
FDP_KJI_D <- get_FDP_KJI(df$W, k_vec=k_list[[4]], v_vec=v_list[[4]])

############ use the number of rejections on the x-axis:
deleted <- which(df$W<=0)

#### save result
df_real <- data.frame(Num_rej=1:length(FDP_bar_vec[-deleted]), FDP_bar_vec=FDP_bar_vec[-deleted], 
                 FDP_hat_vec=FDP_hat_vec[-deleted], 
                 FDP_KJI_A=FDP_KJI_A[-deleted], FDP_KJI_B=FDP_KJI_B[-deleted],
                 FDP_KJI_C=FDP_KJI_C[-deleted], FDP_KJI_D=FDP_KJI_D[-deleted])

save(df_real, file = paste0("New/df_FDP_real.RData"))

############################## Plot:
Method_color <- c("FDP_hat"="black", "KR"="orange1", 
                  "KJI-A"="skyblue1", "KJI-B"="blue",
                  "KJI-C"="darkviolet", "KJI-D"="green3")

plot_all <- ggplot() + geom_point(data = df_real, aes(x = Num_rej, y = FDP_hat_vec, color = "FDP_hat")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_bar_vec, color = "KR")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_A, color = "KJI-A")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_B, color = "KJI-B")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_C, color = "KJI-C")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_D, color = "KJI-D")) +
  scale_color_manual("", values = Method_color) +
  labs( x = 'Number of rejections', y = 'FDP bound' ) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom")
  
plot_zoom_in <- ggplot() + geom_point(data = df_real, aes(x = Num_rej, y = FDP_hat_vec, color = "FDP_hat")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_bar_vec, color = "KR")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_A, color = "KJI-A")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_B, color = "KJI-B")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_C, color = "KJI-C")) +
  geom_point(data = df_real, aes(x = Num_rej, y = FDP_KJI_D, color = "KJI-D")) +
  scale_color_manual("", values = Method_color) +
  labs( x = 'Number of rejections', y = 'FDP bound' ) +
  xlim(c(1,1500)) +
  ylim(c(0,0.25)) +
  geom_hline(yintercept=0.10, linetype='dotted', col = 'black')+
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom")

plot_comb_real <- ggarrange(plot_zoom_in, plot_all,
                            nrow = 1, ncol = 2, common.legend=TRUE)
ggsave("New/real_plots1.png", plot=plot_comb_real, width = 16, height = 8)

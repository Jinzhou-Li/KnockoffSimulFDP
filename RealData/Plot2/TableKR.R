######## The code is largely based on https://github.com/ekatsevi/simultaneous-fdp

setwd("RealData/simultaneous-fdp-master")

# set up workspace
source("setup.R")

# analysis parameters
alpha = 0.05      # confidence level for simultaneous FDP bound
phenotypes = c("height", "bmi", "platelet",    # list of phenotypes
               "sbp", "cvd", "hypothyroidism", 
               "respiratory", "diabetes")

# download KnockoffZoom data, if necessary
download_KZ_data(data_dir)

# read in data
df = read_KZ_data(data_dir, phenotypes)

# compute knockoffs FDP multipler
C = -log(alpha)/log(2-alpha)

# compute FDP-hat and FDP-bar
FDP_hat_vec = unname(unlist(sapply(phenotypes, 
                                   function(phenotype)(get_FDP_hat(df %>% filter(Phenotype == phenotype) %>% pull(W))))))
FDP_bar_vec = unname(unlist(sapply(phenotypes, 
                                   function(phenotype)(get_FDP_bar(df %>% filter(Phenotype == phenotype) %>% pull(W), C)))))
df$FDP_hat = FDP_hat_vec
df$FDP_bar = FDP_bar_vec
df = df %>% filter(W > 0)
df$num_rejections = unname(unlist(sapply(phenotypes, function(phenotype)(1:nrow(df %>% filter(Phenotype == phenotype))))))

# compute numbers of rejections for each Type-I error target
FDR_thresh = sapply(phenotypes, function(phenotype)(get_FDR_thresh(df, phenotype, 0.1)))
FDP_thresh_05 = sapply(phenotypes, function(phenotype)(get_FDP_thresh(df, phenotype, 0.05)))
FDP_thresh_1 = sapply(phenotypes, function(phenotype)(get_FDP_thresh(df, phenotype, 0.1)))
FDP_thresh_2 = sapply(phenotypes, function(phenotype)(get_FDP_thresh(df, phenotype, 0.2)))
df_thresh = data.frame(phenotypes, FDR_thresh, FDP_thresh_2, FDP_thresh_1, FDP_thresh_05)
df_thresh$phenotypes = c("height", "body mass index", "platelet count", "systolic blood pressure", "cardiovascular disease",
                        "hypothyroidism", "respiratory disease", "diabetes")
names(df_thresh) = c("Trait", "FDR $\\leq$ 0.1", "FDP $\\leq$ 0.2", "FDP $\\leq$ 0.1", "FDP $\\leq$ 0.05")
rownames(df_thresh) = c()

# save table
table <- df_thresh %>% kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "") %>%
  add_header_above(c(" " = 2, "with probability 0.95" = 3)) 

save(table, file = paste0("New/Tables/table-KR.RData")) 

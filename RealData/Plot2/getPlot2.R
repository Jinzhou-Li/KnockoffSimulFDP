### extract numbers from tables and make a plot
library(latex2exp)
library(ggplot2)
library(ggpubr)

file_vec <- c("KR", "1", "2", "3", "4")
Trait_vec <- c("height", "body mass index", "platelet count", "systolic blood pressure", 
               "cardiovascular disease", "hypothyroidism", "respiratory disease", "diabetes")
PlotCase_vec <- c("FDP controlled at 0.2", "FDP controlled at 0.1", "FDP controlled at 0.05")
Method_vec <- c("KR", "KJI-A", "KJI-B", "KJI-C", "KJI-D","FDR-control")

Trait_df <- c()
PlotCase_df <- c()
Num_disc_df <- c()
Method_df <- c()

for (file_name_index in 1:length(file_vec)) {
  Table_name <- paste0("RealData/simultaneous-fdp-master/New/Tables/table-",
                       file_vec[file_name_index],".RData")
  load(Table_name)
  
  # FDR result is the same for all files, only need to record for the first one
  if(file_name_index==1){
    for (trait_index in 1:length(Trait_vec)) {
      taget_line <- attr(table, "kable_meta")$content[trait_index+1]
      gfg_numbers <- regmatches(taget_line, gregexpr("[[:digit:]]+", taget_line))
      as.numeric(unlist(gfg_numbers))[-1]
      
      # do not record FDR result
      Trait_df <- c(Trait_df, rep(Trait_vec[trait_index],3))
      PlotCase_df <- c(PlotCase_df, PlotCase_vec)
      Num_disc_df <- c(Num_disc_df, rep(as.numeric(unlist(gfg_numbers))[1],3))
      Method_df <- c(Method_df, rep(Method_vec[6],3))
    }
  }
  
  for (trait_index in 1:length(Trait_vec)) {
    taget_line <- attr(table, "kable_meta")$content[trait_index+1]
    gfg_numbers <- regmatches(taget_line, gregexpr("[[:digit:]]+", taget_line))
    
    # do not record FDR result
    Trait_df <- c(Trait_df, rep(Trait_vec[trait_index], 3) )
    PlotCase_df <- c(PlotCase_df, PlotCase_vec)
    Num_disc_df <- c(Num_disc_df, as.numeric(unlist(gfg_numbers))[-1])
    Method_df <- c(Method_df, rep(Method_vec[file_name_index],3))
  }
}

df <- data.frame(Trait=Trait_df, PlotCase=PlotCase_df, Num_disc=Num_disc_df, Method=Method_df)
df$Trait <- factor(df$Trait, levels=Trait_vec)
df$PlotCase <- factor(df$PlotCase, levels=PlotCase_vec)
df$Method <- factor(df$Method, levels=Method_vec)

Method_color <- c("FDR-control"="black", "KR"="orange1", 
                  "KJI-A"="skyblue1", "KJI-B"="blue", 
                  "KJI-C"="darkviolet",  "KJI-D"="green3")
Method_pointshape <- c("FDR-control"=2, "KR"=1, "KJI-A"=7, "KJI-B"=8, 
                       "KJI-C"=9,  "KJI-D"=10)

##### Plot
plot_gg <- ggplot(df, aes(x=Trait, y=Num_disc, color=Method, shape=Method)) +
  geom_point(size=2, position=position_dodge(width=0.4)) +
  scale_color_manual(values=Method_color) +
  scale_shape_manual(values = Method_pointshape) +
  labs(x = 'Trait', y = 'Number of discoveries') +
  theme_bw(base_size = 15) + 
  theme(legend.title= element_blank(),
        legend.justification = c(0, 1), legend.position = c(0.6, 0.995),
        legend.text = element_text(size=9), axis.text.x=element_text(angle=45, hjust=1)
        ) 

Facet_plot <- plot_gg + facet_grid(rows = vars(PlotCase)) + 
  guides(col = guide_legend(ncol = 2))

ggsave("RealData/Plot2/Plot2_comb.png", plot = Facet_plot, 
       width =18, height = 22, units = "cm", bg="white")

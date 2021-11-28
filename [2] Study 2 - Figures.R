###############################################################################################################
###############################################################################################################
###Uncertainty Models - Figures.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

#Loading libraries
library(tidyverse); library(ggpubr)

#Resolution and filter
Raw_DEM_reso = 10
Filter_type = 40

Dataset1_size = 300

#Set directory
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2/", Dataset1_size)

#Create directory for output
if(!dir.exists(paste0(wd, "/plots"))){dir.create(paste0(wd, "/plots"))}

#Soils variables analysed
Plain_Response_names <- c("M", "T")
Full_text_names <- c("Moisture Regime", "Textural Class")

#####
#####Quantifying model improvement between random and high entropy
if(!dir.exists(paste0(wd, "/plots/Accuracy_Kappa_comparison"))){dir.create(paste0(wd, "/plots/Accuracy_Kappa_comparison"))}

#Performance improvment over Round 1 - Individual plots
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)

for (i in 1:length(unique(df$Variable))){
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)

df <- df[df$Variable == unique(df$Variable)[i],]
df$PlotX <- paste(df$Uncertainty, df$Data_pts_added, sep = "_")



#Violin - Accuracy
#Get one observation of R1 accuracy per treatment per uncertainty metric - then you can compare them on the same graph
R1_df <- group_by(df, Time, Data_pts_added) %>% summarise(Round1_Accuracy = mean(Round1_Accuracy), Round1_Kappa = mean(Round1_Kappa)); R1_df$Round1 = "Round 1"

R1_A <- ggplot()+
  geom_violin(data = R1_df, aes(x = Round1, y = Round1_Accuracy*100), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Accuracy*100) - 5, max(df$Round2_Accuracy*100) + 5)+
  labs(x = " ", y = "Round 1 Accuracy (%)")+
  theme(plot.margin = unit(c(2.3,1,0,1), "mm"))

R2_A <- ggplot()+
  geom_violin(data = df, aes(x = factor(PlotX, levels = unique(paste(df$Uncertainty, rep(unique(df$Data_pts_added), each = 2), sep = "_"))),
                             y = Round2_Accuracy*100, fill = Uncertainty),
              draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Accuracy*100) - 5, max(df$Round2_Accuracy*100) + 5)+
  stat_summary(fun.y=mean, geom="point")+
  labs(y = "Accuracy (%)",
       x = paste(unique(df$Uncertainty_metric), "treatment of additional points"))

Aplot <- ggarrange(R1_A, R2_A, widths = c(1,7), nrow = 1)
Aplot <- annotate_figure(Aplot, top = text_grob(paste0(Full_text_names[i]), face = "bold", size = 16))

ggsave(paste0(wd, "/plots/Accuracy_Kappa_comparison/", unique(df$Uncertainty_metric), "_", Plain_Response_names[i], "_Accuracy_ComparisonViolin.png"), width = 10, height = 6)

#Violin - Kappa
R1_K <- ggplot()+
  geom_violin(data = R1_df, aes(x = Round1, y = Round1_Kappa), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Kappa) - 0.05, max(df$Round2_Kappa) + 0.05)+
  labs(x = " ", y = "Round 1 Kappa")+
  theme(plot.margin = unit(c(2.3,1,0,1), "mm"))

R2_K <- ggplot()+
  geom_violin(data = df, aes(x = factor(PlotX, levels = unique(paste(df$Uncertainty, rep(unique(df$Data_pts_added), each = 2), sep = "_"))),
                             y = Round2_Kappa, fill = Uncertainty),
              draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Kappa) - 0.05, max(df$Round2_Kappa) + 0.05)+
  stat_summary(fun.y=mean, geom="point")+
  labs(y = "Kappa",
       x = paste(unique(df$Uncertainty_metric), "treatment of additional points"))

Kplot <- ggarrange(R1_K, R2_K, widths = c(1,7), nrow = 1)
Kplot <- annotate_figure(Kplot, top = text_grob(paste0(Full_text_names[i]), face = "bold", size = 16))

ggsave(paste0(wd, "/plots/Accuracy_Kappa_comparison/", unique(df$Uncertainty_metric), "_", Plain_Response_names[i], "_Kappa_ComparisonViolin.png"), width = 10, height = 6)

}

#####
#####



#####
####
###
##
#END
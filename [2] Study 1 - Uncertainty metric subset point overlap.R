###############################################################################################################
###############################################################################################################
###[2] Similarity in Uncertainty metric points.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

library(tidyverse)

#####Specifying model and study parameters - user inputs
Raw_DEM_reso = 10
Filter_type = 40

#Set directory
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1")

if(!dir.exists(paste0(wd, "/plots"))){dir.create(paste0(wd, "/plots"))}
if(!dir.exists(paste0(wd, "/plots/D2_unc_similarity"))){dir.create(paste0(wd, "/plots/D2_unc_similarity"))}

df <- read.csv(paste0(wd, "/D2_unc_similarity/D2_uncertainty_similarity.csv"))

#M
M_df <- df[df$Variable == "M",]

M_df <- pivot_longer(M_df, cols = 3:6, names_to = "Metric_Combination", values_to = "Points_similar")
M_df$Metric_Combination <- gsub(M_df$Metric_Combination, pattern = "_Similarity", replacement = "")

ggplot()+
  geom_violin(data = M_df, aes(x = factor(Metric_Combination, levels = c("IU_EU_CI", "IU_CI", "EU_CI", "IU_EU")), y = Points_similar),
              draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(title = "Similarity in high uncertainty points selected between metrics - Moisture Regime",
       x = "Uncertainty metric combination", 
       y = "Points Similar")

ggsave(paste0(wd, "/plots/D2_unc_similarity/M_Similarity.png"), width = 10, height = 6)

write.csv(group_by(M_df, Variable, Metric_Combination) %>% summarize(mean(Points_similar)), paste0(wd, "/plots/D2_unc_similarity/M_Similarity.csv"), row.names = F) 

#T
T_df <- df[df$Variable == "T",]

T_df <- pivot_longer(T_df, cols = 3:6, names_to = "Metric_Combination", values_to = "Points_similar")
T_df$Metric_Combination <- gsub(T_df$Metric_Combination, pattern = "_Similarity", replacement = "")

ggplot()+
  geom_violin(data = T_df, aes(x = factor(Metric_Combination, levels = c("IU_EU_CI", "IU_CI", "EU_CI", "IU_EU")), y = Points_similar),
              draw_quantiles = c(0.25, 0.5, 0.75))+
  labs(title = "Similarity in high uncertainty points selected between metrics - Textural Class",
       x = "Uncertainty metric combination", 
       y = "Points Similar")

ggsave(paste0(wd, "/plots/D2_unc_similarity/T_Similarity.png"), width = 10, height = 6)

write.csv(group_by(T_df, Variable, Metric_Combination) %>% summarize(mean(Points_similar)), paste0(wd, "/plots/D2_unc_similarity/T_Similarity.csv"), row.names = F) 
#####
####
###
##
#END

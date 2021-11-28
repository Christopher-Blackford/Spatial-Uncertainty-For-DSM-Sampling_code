###############################################################################################################
###############################################################################################################
###Uncertainty Models - Figures.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

#Loading libraries
library(tidyverse); library(ggpubr)

#Resolution and filter
Raw_DEM_reso = 10
Filter_type = 40

#Set directory
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1")

#Create directory for output
if(!dir.exists(paste0(wd, "/plots"))){dir.create(paste0(wd, "/plots"))}

#Soils variables analysed
Plain_Response_names <- c("M", "T")
Full_text_names <- c("Moisture Regime", "Textural Class")

#The way the simulations were run, a random treatment was run for each uncertainty metric
#This allowed comparing each metric to a control regardless of the way the "high uncertainty" treatment was defined.
#However, since we settled on using the top 500 points for each uncertainty metric, random sampling 500 points should mean there is
#no difference between "random treatments" for each uncertainty metric. For example, if we had used an uncertainty threshold value (e.g 0.8),
#then the random datasets would be of different sizes and couldn't be pooled. Since we used the point threshold, we are randomly removing
#2/3 of our random treatments before analysis so each treatment can be compared to the same control.

#####
#####Quantifying model improvement between random and high entropy
if(!dir.exists(paste0(wd, "/plots/Accuracy_Kappa_comparison"))){dir.create(paste0(wd, "/plots/Accuracy_Kappa_comparison"))}

for (i in 1:length(Plain_Response_names)){
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)
df <- df[df$Variable == unique(df$Variable)[i],]
  
#Get one observation of R1 accuracy per treatment per uncertainty metric - then you can compare them on the same graph
R1_df <- group_by(df, Time) %>% summarise(Round1_Accuracy = mean(Round1_Accuracy), Round1_Kappa = mean(Round1_Kappa)); R1_df$Round1 = "Round 1"
  
#Finding out correct x-axis order of uncertainty metrics based on performance of high uncertainty treatment
level_order <- summarise(group_by(df[df$Uncertainty == "High",], Uncertainty_metric), level_order = mean(Round2_Accuracy))
level_order <- as.character(c("Random", level_order[order(level_order$level_order, decreasing = F),][[1]]))

R1_A <- ggplot()+
  geom_violin(data = R1_df, aes(x = Round1, y = Round1_Accuracy*100), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Accuracy*100) - 5, max(df$Round2_Accuracy*100) + 5)+
  labs(x = " ", y = "Round 1 Accuracy (%)")+
  theme(plot.margin = unit(c(2.3,1,0,1), "mm"))

R2_A <- ggplot()+
  geom_violin(data = df, aes(x = factor(df$Uncertainty_metric, levels = c("Random", "IU", "EU", "CI")), y = Round2_Accuracy*100, fill = Uncertainty), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Accuracy*100) - 5, max(df$Round2_Accuracy*100) + 5)+
  labs(x = "Uncertainty metric", y = "Round 2 Accuracy (%)")

R1_K <- ggplot()+
  geom_violin(data = R1_df, aes(x = Round1, y = Round1_Kappa), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Kappa) - 0.05, max(df$Round2_Kappa) + 0.05)+
  labs(x = " ", y = "Round 1 Kappa")+
  theme(plot.margin = unit(c(2.3,1,0,1), "mm"))

R2_K <- ggplot()+
  geom_violin(data = df, aes(x = factor(df$Uncertainty_metric, levels = c("Random", "IU", "EU", "CI")), y = Round2_Kappa, fill = Uncertainty), draw_quantiles = c(0.25, 0.5, 0.75), trim = F)+
  ylim(min(df$Round1_Kappa) - 0.05, max(df$Round2_Kappa) + 0.05)+
  labs(x = "Uncertainty metric", y = "Round 2 Kappa")
  
Aplot <- ggarrange(R1_A, R2_A, widths = c(1,5), nrow = 1)
Aplot <- annotate_figure(Aplot, top = text_grob(paste0(Full_text_names[i]), face = "bold", size = 16))
ggsave(plot = Aplot, paste0(wd, "/plots/Accuracy_Kappa_comparison/", Plain_Response_names[i], "_Accuracy_ComparisonViolin.png"), width = 10, height = 6)

Kplot <- ggarrange(R1_K, R2_K, widths = c(1,5), nrow = 1)
Kplot <- annotate_figure(Kplot, top = text_grob(paste0(Full_text_names[i]), face = "bold", size = 16))
ggsave(plot = Kplot, paste0(wd, "/plots/Accuracy_Kappa_comparison/", Plain_Response_names[i], "_Kappa_ComparisonViolin.png"), width = 10, height = 6)
}
#####
#####


#####
#####Uncertainty covaring with soil attribute class at a site
if(!dir.exists(paste0(wd, "/plots/Soil_Unc_graph"))){dir.create(paste0(wd, "/plots/Soil_Unc_graph"))}
df <- read.csv(paste0(wd, "/D2Subselection/M_Classes_in_Dsub.csv"))
#M
df$MR <-  gsub(pattern = "C", replacement = "", df$MR)
df$MR <-  gsub(pattern = "N", replacement = "-", df$MR)
df$MR <-  as.factor(df$MR)

for (j in 1:length(unique(df$Uncertainty_metric))){
df2 <- df[df$Uncertainty_metric == unique(df$Uncertainty_metric)[j],]

trendline <- df2 %>% 
  group_by(MR, Most_Frequent, Uncertainty) %>%
  summarise(D2_dist = mean(D2_dist)*100, DSub_freq_dif = mean(DSub_freq_dif)*100, Random_freq_dif = mean(Random_freq_dif)*100)

trendline <- trendline[order(trendline$Most_Frequent),]

write.csv(trendline, paste0(wd, "/plots/Soil_Unc_graph/", unique(df2$Uncertainty_metric), "_M_SoilClassSelectionByUnc_Avg.csv"), row.names = F)

colors <- c("Random" = "blue", "High" = "red")

ggplot()+
  geom_point(data = df2, aes(x = fct_reorder(df2$MR, df2$Most_Frequent), y = 100*DSub_freq_dif, col = "High"))+
  geom_line(data = trendline, aes(x = Most_Frequent, y = DSub_freq_dif, col = "High"))+
  
  geom_point(data = df2, aes(x = fct_reorder(df2$MR, df2$Most_Frequent), y = 100*Random_freq_dif, col = "Random"))+
  geom_line(data = trendline, aes(x = Most_Frequent, y = Random_freq_dif, col = "Random"))+
  labs(title = paste("Distribution of Moisture Regime values -", unique(df2$Uncertainty_metric), "metric"),
       x = "Moisture Regime Class (Most to Least Common)",
       y = "% in Dataset 2 - % in Treatment dataset",
       colour = "Treatment")+
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.25))

ggsave(paste0(wd, "/plots/Soil_Unc_graph/", unique(df2$Uncertainty_metric), "_M_SoilClassSelectionByUnc.png"), width = 10, height = 6)
}

#T
df <- read.csv(paste0(wd, "/D2Subselection/T_Classes_in_Dsub.csv"))
for (j in 1:length(unique(df$Uncertainty_metric))){
df2 <- df[df$Uncertainty_metric == unique(df$Uncertainty_metric)[j],]
  
trendline <- df2 %>% 
  group_by(TEXT, Most_Frequent, Uncertainty) %>%
  summarise(D2_dist = mean(D2_dist)*100, DSub_freq_dif = mean(DSub_freq_dif)*100, Random_freq_dif = mean(Random_freq_dif)*100)
  
trendline <- trendline[order(trendline$Most_Frequent),]

write.csv(trendline, paste0(wd, "/plots/Soil_Unc_graph/", unique(df2$Uncertainty_metric), "_T_SoilClassSelectionByUnc_Avg.csv"), row.names = F)

colors <- c("Random" = "blue", "High" = "red")

ggplot()+
  geom_point(data = df2, aes(x = fct_reorder(df2$TEXT, df2$Most_Frequent), y = 100*DSub_freq_dif, col = "High"))+
  geom_line(data = trendline, aes(x = Most_Frequent, y = DSub_freq_dif, col = "High"))+
  
  geom_point(data = df2, aes(x = fct_reorder(df2$TEXT, df2$Most_Frequent), y = 100*Random_freq_dif, col = "Random"))+
  geom_line(data = trendline, aes(x = Most_Frequent, y = Random_freq_dif, col = "Random"))+
  labs(title = paste("Distribution of Textural Class values -", unique(df2$Uncertainty_metric), "metric"),
       x = "Textural Class (Most to Least Common)",
       y = "% in Dataset 2 - % in Treatment dataset",
       colour = "Treatment")+
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.25))

ggsave(paste0(wd, "/plots/Soil_Unc_graph/", unique(df2$Uncertainty_metric), "_T_SoilClassSelectionByUnc.png"), width = 10, height = 6)
}
#####
#####


#############################################################################################
#############################################################################################
#####Supplementary

#####
#####Uncertainty histogram of Dataset 2
Uncertainty_metric_names <- c("IU", "EU", "CI")

if(!dir.exists(paste0(wd, "/plots/D2_unc_dist_graphs"))){dir.create(paste0(wd, "/plots/D2_unc_dist_graphs"))}
df <- read.csv(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), stringsAsFactors = FALSE)

for (i in 1:length(unique(df$Variable))){
  df <- read.csv(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), stringsAsFactors = FALSE)
  df <- df[df$Variable == unique(df$Variable)[i],]
  
  for(j in 1:length(Uncertainty_metric_names)){
  #Total distribution 
  ggplot(df, aes(x = get(Uncertainty_metric_names[j])))+ xlab(Uncertainty_metric_names[j])+
    geom_histogram(aes(), alpha = 0.5)
  #Standardized across multiple simulations
  ggplot(df, aes(x = get(Uncertainty_metric_names[j])))+
    geom_histogram(aes(y=..count../length(unique(df$time))), alpha = 0.5)+
    labs(title = Full_text_names[i], x = Uncertainty_metric_names[j], y = "Standardized frequency")
  
  ggsave(paste0(wd, "/plots/D2_unc_dist_graphs/", Uncertainty_metric_names[j], "_", Plain_Response_names[i], "_D2_entropy_dist.png"), width = 10, height = 6)
  }
}
#####
#####

#####
#####Uncertainty histogram of DSub comparing High and Random Uncertainty
if(!dir.exists(paste0(wd, "/plots/High_vs_Random_hist"))){dir.create(paste0(wd, "/plots/High_vs_Random_hist"))}
for(i in 1:length(Plain_Response_names)){
  for(j in 1:length(Uncertainty_metric_names)){
  df <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metric_names[j], ".csv"), stringsAsFactors = FALSE)
  df <- df[df$Variable == Plain_Response_names[i],]
  
  #Histogram
  ggplot(df)+
    geom_histogram(aes(x=DSub_HighUncValues), fill = "red", alpha = 0.5)+
    geom_histogram(aes(x=DSub_RandomUncValues), fill = "blue", alpha = 0.5)+
    labs(title = Full_text_names[i],
         x = paste0("DSub ", Uncertainty_metric_names[j]),
         y = "Count")
  
  ggsave(paste0(wd, "/plots/High_vs_Random_hist/", Uncertainty_metric_names[j], "_", Plain_Response_names[i], "_DSub_entropies.png"), width = 10, height = 6)
  }
}
#####
#####

#####
#####Histogram of entropy cutoff values using point threshold
if(!dir.exists(paste0(wd, "/plots/D2_points_added"))){dir.create(paste0(wd, "/plots/D2_points_added"))}

for(i in 1:length(Plain_Response_names)){
  for(j in 1:length(Uncertainty_metric_names)){
  df <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metric_names[j], ".csv"), stringsAsFactors = FALSE)
  df <- df[df$Variable == Plain_Response_names[i],]
  df <- distinct(df, time, Variable, DSub_Unc_threshold)
  
  #Histogram
  ggplot(df)+
    geom_histogram(aes(x=DSub_Unc_threshold), alpha = 0.5)+
    labs(title = paste(Plain_Response_names[i], "DSub ", Uncertainty_metric_names[j], " threshold values"),
         y = "Count")
  ggsave(paste0(wd, "/plots/D2_points_added/", Uncertainty_metric_names[j], "_", Plain_Response_names[i], "_DSub_unc_thresholds.png"), width = 10, height = 6)
}
}
#####
#####

#############################################################################################
#############################################################################################

#####
####
###
##
#END
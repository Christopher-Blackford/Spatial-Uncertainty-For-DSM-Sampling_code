###############################################################################################################
###############################################################################################################
###Uncertainty Models_Stats.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

#Loading libraries
library(tidyverse)

#Resolution and filter
Raw_DEM_reso = 10
Filter_type = 40

Dataset1_size = 300

#Set directory
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2/", Dataset1_size)

#Create directory for output
if(!dir.exists(paste0(wd, "/stats"))){dir.create(paste0(wd, "/stats"))}
if(!dir.exists(paste0(wd, "/stats/mlr_plots"))){dir.create(paste0(wd, "/stats/mlr_plots"))}

#####
#####Quantifying model accuracy improvement between random and high entropy
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)

#####
#####
#Moisture Regime
M_df <- df[df$Variable == "M",]

#We can manually convert the p-values to the one-sided case since we expect the Random treatment to be lower than the High treatment?
M_A_Reg <- lm(Round2_Accuracy*100 ~ Uncertainty + Data_pts_added + Uncertainty*Data_pts_added, data = M_df)
summary(M_A_Reg)
write.csv(cbind(summary(M_A_Reg)$coefficients, R2 = summary(M_A_Reg)$r.squared), paste0(wd,"/stats/MR_Accuracy_mlr.csv"))

x = c(300:1200)

#Random line
y_Rand = M_A_Reg$coefficients[[3]]*x + M_A_Reg$coefficients[[4]]*x + M_A_Reg$coefficients[[2]] + M_A_Reg$coefficients[[1]]

#High line
y_High = M_A_Reg$coefficients[[3]]*x + M_A_Reg$coefficients[[1]]

ggplot()+
  geom_point(data=M_df, aes(x = Data_pts_added, y = Round2_Accuracy*100, col = Uncertainty))+
  geom_line(aes(x = x, y = y_Rand), col = "blue")+
  geom_line(aes(x = x, y = y_High), col = "red")+
  labs(title = "Mositure Regime",
       y = "Accuracy",
       x = "Data pts added")
ggsave(paste0(wd, "/stats/mlr_plots/MR_Accuracy_mlr.png"), width = 10, height = 6)


M_K_Reg <- lm(Round2_Kappa ~ Uncertainty + Data_pts_added + Uncertainty*Data_pts_added, data = M_df)
summary(M_K_Reg)
write.csv(cbind(summary(M_K_Reg)$coefficients, R2 = summary(M_K_Reg)$r.squared), paste0(wd,"/stats/MR_Kappa_mlr.csv"))

#Random line
y_Rand = M_K_Reg$coefficients[[3]]*x + M_K_Reg$coefficients[[4]]*x + M_K_Reg$coefficients[[2]] + M_K_Reg$coefficients[[1]]

#High line
y_High = M_K_Reg$coefficients[[3]]*x + M_K_Reg$coefficients[[1]]

ggplot()+
  geom_point(data=M_df, aes(x = Data_pts_added, y = Round2_Kappa, col = Uncertainty))+
  geom_line(aes(x = x, y = y_Rand), col = "blue")+
  geom_line(aes(x = x, y = y_High), col = "red")+
  labs(title = "Mositure Regime",
       y = "Kappa",
       x = "Data pts added")
ggsave(paste0(wd, "/stats/mlr_plots/MR_Kappa_mlr.png"), width = 10, height = 6)



#####
#####
#Textural Class
T_df <- df[df$Variable == "T",]

T_A_Reg <- lm(Round2_Accuracy*100 ~ Uncertainty + Data_pts_added + Uncertainty*Data_pts_added, data = T_df)
summary(T_A_Reg)
write.csv(cbind(summary(T_A_Reg)$coefficients, R2 = summary(T_A_Reg)$r.squared), paste0(wd,"/stats/T_Accuracy_mlr.csv"))

#Random line
y_Rand = T_A_Reg$coefficients[[3]]*x + T_A_Reg$coefficients[[4]]*x + T_A_Reg$coefficients[[2]] + T_A_Reg$coefficients[[1]]

#High line
y_High = T_A_Reg$coefficients[[3]]*x + T_A_Reg$coefficients[[1]]

ggplot()+
  geom_point(data=T_df, aes(x = Data_pts_added, y = Round2_Accuracy*100, col = Uncertainty))+
  geom_line(aes(x = x, y = y_Rand), col = "blue")+
  geom_line(aes(x = x, y = y_High), col = "red")+
  labs(title = "Textural Class",
       y = "Accuracy",
       x = "Data pts added")
ggsave(paste0(wd, "/stats/mlr_plots/T_Accuracy_mlr.png"), width = 10, height = 6)

#
T_K_Reg <- lm(Round2_Kappa ~ Uncertainty + Data_pts_added + Uncertainty*Data_pts_added, data = T_df)
summary(T_K_Reg)
write.csv(cbind(summary(T_K_Reg)$coefficients, R2 = summary(T_K_Reg)$r.squared), paste0(wd,"/stats/T_Kappa_mlr.csv"))

#Random line
y_Rand = T_K_Reg$coefficients[[3]]*x + T_K_Reg$coefficients[[4]]*x + T_K_Reg$coefficients[[2]] + T_K_Reg$coefficients[[1]]

#High line
y_High = T_K_Reg$coefficients[[3]]*x + T_K_Reg$coefficients[[1]]

ggplot()+
  geom_point(data=T_df, aes(x = Data_pts_added, y = Round2_Kappa, col = Uncertainty))+
  geom_line(aes(x = x, y = y_Rand), col = "blue")+
  geom_line(aes(x = x, y = y_High), col = "red")+
  labs(title = "Textural Class",
       y = "Kappa",
       x = "Data pts added")
ggsave(paste0(wd, "/stats/mlr_plots/T_Kappa_mlr.png"), width = 10, height = 6)




######
#Generate summary MLR table

write.csv(rbind(cbind(summary(M_A_Reg)$coefficients, R2 = summary(M_A_Reg)$r.squared),
                cbind(summary(M_K_Reg)$coefficients, R2 = summary(M_K_Reg)$r.squared),
                cbind(summary(T_A_Reg)$coefficients, R2 = summary(T_A_Reg)$r.squared),
                cbind(summary(T_K_Reg)$coefficients, R2 = summary(T_K_Reg)$r.squared)
                ),
          paste0(wd, "/stats/MLR_results.csv")
          )

#####
####
###
##
#END
###############################################################################################################
###############################################################################################################
###Uncertainty Models_Stats.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

#Loading libraries
library(tidyverse)

#Resolution and filter
Raw_DEM_reso = 10
Filter_type = 40

Supplemental_stats = FALSE
Uncertainty_metric_names <- c("IU", "EU", "CI")
Plain_Response_names <- c("M", "T")

#Set directory
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1")

#Create directory for output
if(!dir.exists(paste0(wd, "/stats"))){dir.create(paste0(wd, "/stats"))}

############################
############################
#####Quantifying model improvement from Round 1 to Round 2 Random treatment - Paired t test
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)

#Moisture Regime
M_df <- df[df$Variable == "M",]
M_df <- M_df[M_df$Uncertainty == "Random",]

if(min(c(shapiro.test(M_df$Round2_Accuracy)$p.value, shapiro.test(M_df$Round1_Accuracy)$p.value)) > 0.05){print("Normality assumption met")}
if(bartlett.test(list(M_df$Round2_Accuracy, M_df$Round1_Accuracy))$p.value > 0.05){print("Homogeneity of variances met")} #safe to comebine standard deviations for cohen's d
t_results_Accuracy <- t.test(x = M_df$Round2_Accuracy, y = M_df$Round1_Accuracy, paired = TRUE, alternative = "greater")

if(min(c(shapiro.test(M_df$Round2_Kappa)$p.value, shapiro.test(M_df$Round1_Kappa)$p.value)) > 0.05){print("Normality assumption met")}
if(bartlett.test(list(M_df$Round2_Kappa, M_df$Round1_Kappa))$p.value > 0.05){print("Homogeneity of variances met")}
t_results_Kappa <- t.test(x = M_df$Round2_Kappa, y = M_df$Round1_Kappa, paired = TRUE, alternative = "greater")

#estimate of cohen's d based on d = t stat/sqrt(N)

M_A <- data.frame(Variable = "M", 
                  Metric = "Accuracy",
                  Round1 = mean(M_df$Round1_Accuracy),
                  Random_Round2 = mean(M_df$Round2_Accuracy),
                  Diff_in_means = mean(M_df$Round2_Accuracy) - mean(M_df$Round1_Accuracy),
                  t_stat = t_results_Accuracy$statistic[[1]],
                  p_value = t_results_Accuracy$p.value,
                  cohen_d = t_results_Accuracy$statistic[[1]]/sqrt(nrow(M_df)),
                  Notes = paste(capture.output(t_results_Accuracy), collapse = " "))

M_K <- data.frame(Variable = "M",
                  Metric = "Kappa",
                  Round1 = mean(M_df$Round1_Kappa),
                  Random_Round2 = mean(M_df$Round2_Kappa),
                  Diff_in_means = mean(M_df$Round2_Kappa) - mean(M_df$Round1_Kappa),
                  t_stat = t_results_Kappa$statistic[[1]],
                  p_value = t_results_Kappa$p.value,
                  cohen_d = t_results_Kappa$statistic[[1]]/sqrt(nrow(M_df)),
                  Notes = paste(capture.output(t_results_Kappa), collapse = " "))


#Textural Class
T_df <- df[df$Variable == "T",]
T_df <- T_df[T_df$Uncertainty == "Random",]

if(min(c(shapiro.test(T_df$Round2_Accuracy)$p.value, shapiro.test(T_df$Round1_Accuracy)$p.value)) > 0.05){print("Normality assumption met")}
if(bartlett.test(list(T_df$Round2_Accuracy, T_df$Round1_Accuracy))$p.value > 0.05){print("Homogeneity of variances met")}
t_results_Accuracy <- t.test(x = T_df$Round2_Accuracy, y = T_df$Round1_Accuracy, data = T_df, paired = TRUE, alternative = "greater")

if(min(c(shapiro.test(T_df$Round2_Kappa)$p.value, shapiro.test(T_df$Round1_Kappa)$p.value)) > 0.05){print("Normality assumption met")}
if(bartlett.test(list(T_df$Round2_Kappa, T_df$Round1_Kappa))$p.value > 0.05){print("Homogeneity of variances met")}
t_results_Kappa <- t.test(x = T_df$Round2_Kappa, y = T_df$Round1_Kappa, data = T_df, paired = TRUE, alternative = "greater")

T_A <- data.frame(Variable = "T", 
                  Metric = "Accuracy",
                  Round1 = mean(T_df$Round1_Accuracy),
                  Random_Round2 = mean(T_df$Round2_Accuracy),
                  Diff_in_means = mean(T_df$Round2_Accuracy) - mean(T_df$Round1_Accuracy),
                  t_stat = t_results_Accuracy$statistic[[1]],
                  p_value = t_results_Accuracy$p.value,
                  cohen_d = t_results_Accuracy$statistic[[1]]/sqrt(nrow(M_df)),
                  Notes = paste(capture.output(t_results_Accuracy), collapse = " "))

T_K <- data.frame(Variable = "T",
                  Metric = "Kappa",
                  Round1 = mean(T_df$Round1_Kappa),
                  Random_Round2 = mean(T_df$Round2_Kappa),
                  Diff_in_means = mean(T_df$Round2_Kappa) - mean(T_df$Round1_Kappa),
                  t_stat = t_results_Kappa$statistic[[1]],
                  p_value = t_results_Accuracy$p.value,
                  cohen_d = t_results_Kappa$statistic[[1]]/sqrt(nrow(M_df)),
                  Notes = paste(capture.output(t_results_Kappa), collapse = " "))

write.csv(rbind(M_A, M_K, T_A, T_K),
          paste0(wd, "/stats/Round1_to_Random_t_test.csv"),
          row.names = F)
############################
############################


############################
############################
#####Quantifying difference between control and treatment groups
df <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = F)

shapiro.test(df[df$Variable == "M" & df$Uncertainty_metric == "IU",]$Round2_Accuracy)$p.value #Almost violates normality
shapiro.test(df[df$Variable == "M" & df$Uncertainty_metric == "EU",]$Round2_Accuracy)$p.value
shapiro.test(df[df$Variable == "M" & df$Uncertainty_metric == "CI",]$Round2_Accuracy)$p.value #Almost violates normality

shapiro.test(df[df$Variable == "T" & df$Uncertainty_metric == "IU",]$Round2_Accuracy)$p.value
shapiro.test(df[df$Variable == "T" & df$Uncertainty_metric == "EU",]$Round2_Accuracy)$p.value
shapiro.test(df[df$Variable == "T" & df$Uncertainty_metric == "CI",]$Round2_Accuracy)$p.value

#Because normality is almost violated in 2 cases, I am using Wilcox test
#Tieing accuracy values in Wilcox test is fine, it means that the models are identifying the same proportion of soil classes as the correct value (e.g. 400/500 and 800/1000)

#Moisture Regime
M_df <- df[df$Variable == "M",]

wilcox_test_IU <- wilcox.test(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Accuracy, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")
wilcox_test_EU <- wilcox.test(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Accuracy, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")
wilcox_test_CI <- wilcox.test(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Accuracy, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")

#convert form p value to z score
#z<-qnorm(wilcox_test_IU$p.value/2)
#effect <- abs(z)/sqrt(30)
M_A <- data.frame(Variable = "M", 
                  Metric = "Accuracy",
                  Random_mean = mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  IU_mean = mean(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Accuracy),
                  IU_mean_dif = mean(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Accuracy) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  EU_mean = mean(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Accuracy),
                  EU_mean_dif = mean(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Accuracy) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  CI_mean = mean(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Accuracy),
                  CI_mean_dif = mean(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Accuracy) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  IU_p_value = wilcox_test_IU$p.value,
                  IU_effect_r = abs(qnorm(wilcox_test_IU$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "IU",])),
                  EU_p_value = wilcox_test_EU$p.value,
                  EU_effect_r = abs(qnorm(wilcox_test_EU$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "EU",])),
                  CI_p_value = wilcox_test_CI$p.value,
                  CI_effect_r = abs(qnorm(wilcox_test_CI$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "CI",]))
                  )

wilcox_test_IU <- wilcox.test(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Kappa, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")
wilcox_test_EU <- wilcox.test(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Kappa, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")
wilcox_test_CI <- wilcox.test(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Kappa, M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")


M_K <- data.frame(Variable = "M", 
                  Metric = "Kappa",
                  Random_mean = mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  IU_mean = mean(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Kappa),
                  IU_mean_dif = mean(M_df[M_df$Uncertainty_metric == "IU",]$Round2_Kappa) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  EU_mean = mean(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Kappa),
                  EU_mean_dif = mean(M_df[M_df$Uncertainty_metric == "EU",]$Round2_Kappa) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  CI_mean = mean(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Kappa),
                  CI_mean_dif = mean(M_df[M_df$Uncertainty_metric == "CI",]$Round2_Kappa) - mean(M_df[M_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  IU_p_value = wilcox_test_IU$p.value,
                  IU_effect_r = abs(qnorm(wilcox_test_IU$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "IU",])),
                  EU_p_value = wilcox_test_EU$p.value,
                  EU_effect_r = abs(qnorm(wilcox_test_EU$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "EU",])),
                  CI_p_value = wilcox_test_CI$p.value,
                  CI_effect_r = abs(qnorm(wilcox_test_CI$p.value/2))/sqrt(nrow(M_df[M_df$Uncertainty_metric == "CI",]))
                  )


#Textural Class
T_df <- df[df$Variable == "T",]

wilcox_test_IU <- wilcox.test(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Accuracy, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")
wilcox_test_EU <- wilcox.test(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Accuracy, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")
wilcox_test_CI <- wilcox.test(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Accuracy, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy, paired = TRUE, alternative = "greater")

#convert form p value to z score
#z<-qnorm(wilcox_test_IU$p.value/2)
#effect <- abs(z)/sqrt(30)
T_A <- data.frame(Variable = "T", 
                  Metric = "Accuracy",
                  Random_mean = mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  IU_mean = mean(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Accuracy),
                  IU_mean_dif = mean(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Accuracy) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  EU_mean = mean(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Accuracy),
                  EU_mean_dif = mean(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Accuracy) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  CI_mean = mean(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Accuracy),
                  CI_mean_dif = mean(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Accuracy) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Accuracy),
                  IU_p_value = wilcox_test_IU$p.value,
                  IU_effect_r = abs(qnorm(wilcox_test_IU$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "IU",])),
                  EU_p_value = wilcox_test_EU$p.value,
                  EU_effect_r = abs(qnorm(wilcox_test_EU$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "EU",])),
                  CI_p_value = wilcox_test_CI$p.value,
                  CI_effect_r = abs(qnorm(wilcox_test_CI$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "CI",]))
)

wilcox_test_IU <- wilcox.test(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Kappa, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")
wilcox_test_EU <- wilcox.test(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Kappa, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")
wilcox_test_CI <- wilcox.test(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Kappa, T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa, paired = TRUE, alternative = "greater")


T_K <- data.frame(Variable = "T", 
                  Metric = "Kappa",
                  Random_mean = mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  IU_mean = mean(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Kappa),
                  IU_mean_dif = mean(T_df[T_df$Uncertainty_metric == "IU",]$Round2_Kappa) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  EU_mean = mean(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Kappa),
                  EU_mean_dif = mean(T_df[T_df$Uncertainty_metric == "EU",]$Round2_Kappa) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  CI_mean = mean(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Kappa),
                  CI_mean_dif = mean(T_df[T_df$Uncertainty_metric == "CI",]$Round2_Kappa) - mean(T_df[T_df$Uncertainty_metric == "Random",]$Round2_Kappa),
                  IU_p_value = wilcox_test_IU$p.value,
                  IU_effect_r = abs(qnorm(wilcox_test_IU$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "IU",])),
                  EU_p_value = wilcox_test_EU$p.value,
                  EU_effect_r = abs(qnorm(wilcox_test_EU$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "EU",])),
                  CI_p_value = wilcox_test_CI$p.value,
                  CI_effect_r = abs(qnorm(wilcox_test_CI$p.value/2))/sqrt(nrow(T_df[T_df$Uncertainty_metric == "CI",]))
)

write.csv(rbind(M_A, M_K, T_A, T_K),
          paste0(wd, "/stats/Random_to_Treatments_wilcoxon_SR.csv"),
          row.names = F)
############################
############################






if(Supplemental_stats == TRUE){

#####
#####Quantifying differences in High and Random Entropy histograms
df <- read.csv(paste0(wd, "/D2Subselection/DSub_entropy.csv"), stringsAsFactors = FALSE)

#Entropy distributions aren't normal distribution (obviously) so can we tranform before testing of not?
#Assuming not, use Mann-Whiteny U test

#Moisture Regime
M_df <- df[df$Variable == "M",]
M_df$RAND <- sample(1:nrow(M_df), size = nrow(M_df), replace = F)

wilcox.test(M_df$DSub_HighEntropyValues, M_df$DSub_RandomEntropyValues, paired = TRUE)

#Textural Class
T_df <- df[df$Variable == "T",]
T_df$RAND <- sample(1:nrow(T_df), size = nrow(T_df), replace = F)

#I could output these results but honestly the differences are so striking I don't think it's worth it. It's not the main finding, these just justify our approach.
#####
#####

######How does entropy covary with soil attribute class at a site?
M_df <- read.csv(paste0(wd, "/D2Subselection/M_Classes_in_Dsub.csv"), stringsAsFactors = FALSE)

#I don't think the chisq test works here actually unless I convert y axis to means so I just have one count per class.
#or could do t-test on each class
#Or would need to smooth curve and do ks test to show how continuous distribution function is different?
chisq.test(tbl)


T_df <- read.csv(paste0(wd, "/D2Subselection/T_Classes_in_Dsub.csv"), stringsAsFactors = FALSE)

}

#####
####
###
##
#END


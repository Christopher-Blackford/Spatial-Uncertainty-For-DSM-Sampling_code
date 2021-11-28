##############################################################################################################################
####Quantify soil and environment value covariation with uncertainty.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

#####Look at if specific areas are chosen in these low/high/random entropy cells

###T-test and Kolmogorov-Smirnov test relationship between environmental covariates and entropy - note that some of these environmental covariates should not be looked at with a t-test because they are categorical
t_test_result <- NULL
ks_test_result <- NULL
for(t_test_loop in 1:length(names(Covariates))){
  col_index <- which(colnames(DSub_HighUnc) %in% names(Covariates)[t_test_loop])
  
  #If the variances of the 2 datasets are the same you can run a t-test
  if(var.test(Dataset2[,col_index], DSub_HighUnc[,col_index])$p.value > 0.05){
    t_test <- t.test(Dataset2[,col_index], DSub_HighUnc[,col_index])
    if(t_test$p.value < 0.01){ #Report if signifigance is met
      direction <- paste0(c(names(Covariates)[t_test_loop], as.character(round(t_test$estimate[[2]] - t_test$estimate[[1]], digits = 2))), collapse = "_")
      t_test_result <- append(t_test_result, direction)}
  }
  #If the variances of the 2 datasets are not the same, you can run the Kolmogorov-Smirnov
  if(ks.test(Dataset2[,col_index], DSub_HighUnc[,col_index])$p.value < 0.01){#If signifigance for KS test is met, record it
    ks_test_result <- append(ks_test_result, names(Covariates)[t_test_loop])
  } 
  #Many times the ks.test will find "ties" meaning that there were very similar values between the 2 datasets. This tricks the function into
  #thinking you are looking at a categorical rather than continuous variable. Annoying and could work around but since it's only a warning message it's not a big deal
}; rm(t_test); 
#Positive values indicate significantly higher covariate values in Dataset 2 compared to Dataset Subsequent
t_test_result <- paste(t_test_result, collapse = " "); ks_test_result <- paste(ks_test_result, collapse = " ")
###

#Soil Uncertainty analysis 
###Calculate differences in soil attribute class dominance in Dataset 2 vs Dataset Subsequent
Soil_class_FreqUnc <- data.frame(table(Dataset2[,loop_Name]))
Soil_class_FreqUnc$Freq <- Soil_class_FreqUnc$Freq/sum(Soil_class_FreqUnc$Freq)
Soil_class_FreqUnc <- plyr::rename(Soil_class_FreqUnc, c("Var1" = df_Response_names[i], "Freq" = "D2_dist"))
#Random treatment soil attribute distributions
Soil_class_FreqUnc$Random_freq_dist <- data.frame(table(DSub_RandomUnc[,loop_Name]))$Freq/sum(table(DSub_RandomUnc[,loop_Name]))
Soil_class_FreqUnc$Random_freq_dif <- Soil_class_FreqUnc$Random_freq_dist - Soil_class_FreqUnc$D2_dist
#High uncertainty treatment soil attribute distributions
Soil_class_FreqUnc$DSub_freq_dist <- data.frame(table(DSub_HighUnc[,loop_Name]))$Freq/sum(table(DSub_HighUnc[,loop_Name]))
Soil_class_FreqUnc$DSub_freq_dif <- Soil_class_FreqUnc$DSub_freq_dist - Soil_class_FreqUnc$D2_dist
#Ordering
Soil_class_FreqUnc <- Soil_class_FreqUnc[order(-Soil_class_FreqUnc$D2_dist),]
Soil_class_FreqUnc$Most_Frequent <- 1:nrow(Soil_class_FreqUnc)



#Add covariate difference in Dataset 2 vs Dataset Subsequent classes
Soil_class_FreqUnc$t_test_result <- t_test_result
Soil_class_FreqUnc$ks_test_result <- ks_test_result
Soil_class_FreqUnc$time <- time
Soil_class_FreqUnc$Uncertainty <- "High"
Soil_class_FreqUnc$Uncertainty_metric <- Uncertainty_metrics[Unc_metric]
Soil_class_FreqUnc <- Soil_class_FreqUnc %>% dplyr::select(Uncertainty, Uncertainty_metric, everything()) %>% dplyr::select(time, everything())

#Report the differences in soil attribute class dominance and covariate distribution in Dataset 2 vs Dataset Subsequent
if(file.exists(paste0(wd, "/D2Subselection/", Plain_Response_names[i], "_Classes_in_DSub.csv"))){
  DSub_Summary <- read.csv(paste0(wd, "/D2Subselection/", Plain_Response_names[i], "_Classes_in_DSub.csv"), stringsAsFactors = FALSE)
  DSub_Summary <- rbind(DSub_Summary, Soil_class_FreqUnc)
  write.csv(DSub_Summary, paste0(wd, "/D2Subselection/", Plain_Response_names[i], "_Classes_in_DSub.csv"), row.names = F)
}else{write.csv(Soil_class_FreqUnc, paste0(wd, "/D2Subselection/", Plain_Response_names[i], "_Classes_in_DSub.csv"), row.names = F)
}


#Old stuff
###Calculate differences in environmental covariates
#Holling's t-test can look at multivariate differences between sample and population but has issues running if there is multicollinearity between environmental covariates
#library(Hotelling)
#s1 <- data.matrix(Dataset2)
#s2 <- data.matrix(DSub_HighUnc)
#cols = c(1,2,3,4,6,7,8,12,13,14)
#s1_new <- s1[,cols]
#s2_new <- s2[,cols]
#new_temp=hotelling.test(s1_new,s2_new)
##Can also just do t-test loops

#####
####
###
##
#END

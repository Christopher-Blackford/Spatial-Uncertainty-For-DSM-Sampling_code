##############################################################################################################################
##############################################################################################################################
###Uncertainty Models.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)
start_time <- Sys.time()

#Loading libraries
library(raster); library(rgdal); library(rgeos)
library(caret); library(doParallel); library(tidyverse)

#####Specifying model and study parameters - user inputs
Raw_DEM_reso = 10
Filter_type = 40

repeats = 2
repeat_count = 0

Dataset1_size = 500
Dataset2_initsize = 6000
Cutoff = 500

#####
#Create directories for response dataset
if(!dir.exists(paste0("./Uncertainty_Models/"))){dir.create(paste0("./Uncertainty_Models/"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data"))}

#Create directories for model output depending on if you are using point of entropy threshold approach
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1"))}
#Short version of working directory which is useful for later
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1")

###Specific directory structure for uncertainty model output
#Model comparison
if(!dir.exists(paste0(wd, "/Model_Comparison"))){dir.create(paste0(wd, "/Model_Comparison"))}
if(!dir.exists(paste0(wd, "/Model_Comparison/models"))){dir.create(paste0(wd, "/Model_Comparison/models"))}
if(!dir.exists(paste0(wd, "/Model_Comparison/models/Round1"))){dir.create(paste0(wd, "/Model_Comparison/models/Round1"))}
if(!dir.exists(paste0(wd, "/Model_Comparison/models/Round2"))){dir.create(paste0(wd, "/Model_Comparison/models/Round2"))}
#Understanding dataset 2 properties and dataset 2 subselection
if(!dir.exists(paste0(wd, "/D2_unc_distribution"))){dir.create(paste0(wd, "/D2_unc_distribution"))}
if(!dir.exists(paste0(wd, "/D2_unc_similarity"))){dir.create(paste0(wd, "/D2_unc_similarity"))}
if(!dir.exists(paste0(wd, "/D2Subselection"))){dir.create(paste0(wd, "/D2Subselection"))}

#Running in parallel
detectCores() #Should be 12 on SOIL-BOT
registerDoParallel(makeCluster(detectCores()-1))
getDoParWorkers()
#How to stop cluster?
#stopImplicitCluster() or #endCluster(cl) or #stopCluster(cl)

###Acquiring continuous environmental covariates
source("./[subcode] - Load soil pedon dataset and summary histogram.R")

###Setting up models
#Number of decision trees for the random forest
num_trees <- 751

#For each tree, how many variables are randomly sampled at each node
mtry <- c(3,5,7,9,11,13,15) #Rule of thumb is that mtry is sqrt(# of predictors) = round(sqrt(38))
tune_grid <- data.frame(mtry = mtry) #Creating a tuning grid for input into train(tuneGrid = tune)

#How many folds of training data to analyse
num_folds <- 10
#How many times per fold you analyse the data
num_repeats <- 10

#Train function settings
fitControl <- trainControl(
  method = "repeatedcv",    #repeated cross-validation      
  number = num_folds,       #number of folds  
  repeats = num_repeats,    #number of times you divide datasets into folds
  returnResamp = "all",     #tells you which sample in train() function were used to train (bootstrapped) and which were left out (see: RF_Train$control$index)       
  allowParallel = TRUE,
  savePredictions = "final",
  classProbs = TRUE
  )

#Making legitimate Moisture regime class
#df$MR <- as.character(df$MOISTURE_R)
#for (i in 1:nrow(df)){
 # if(df$MR[i] != ""){df$MR[i] <- paste0("C", df$MR[i])}
 #}

#df$MR <- gsub(pattern = "C-1", replacement = "CN1", df$MR) 
#df$MR <- as.factor(df$MR)

####Looping over predictor dataset
Plain_Response_names <- c("M", "T")
df_Response_names <- c("MR", "TEXT")

##################################################################################
##################################################################################
#####Repeating for multiple simulations of drawing and redrawing soil points
repeat{
  repeat_count = repeat_count+1

  #Storing time to record each unique simulation
  time = gsub(Sys.time(), pattern = ":", replacement = "_"); time = gsub(time, pattern = " ", replacement = "_")
  
##################################################################################
##################################################################################
#####Running Uncertainty Models
print("Running Random Forest for Uncertainty Models")
for(i in 1:(length(Plain_Response_names))){ 
  Train_df <- df
  loop_Name <- which(colnames(Train_df) %in% df_Response_names[i]) #Column position of MR, TEXT
  
  #Universal cleaning (removing spaces and converting to factor)
  Train_df[,loop_Name] <- as.character(Train_df[,loop_Name])
  Train_df <- Train_df[(Train_df[,loop_Name] != "") & (Train_df[,loop_Name] != " "),]
  Train_df[,loop_Name] <- as.factor(Train_df[,loop_Name])
  
  #Splitting data for training/testing
  partition <- createDataPartition(Train_df[,loop_Name], p=Dataset1_size/nrow(Train_df), list=F) 
  #May need to add a line here to continue reducing dataset size so you aren't including excluded points in later parts of the analysis
  Dataset1 <- Train_df[partition,]; NextRound <- Train_df[-partition,]
  partition <- createDataPartition(NextRound[,loop_Name], p=Dataset2_initsize/nrow(NextRound), list=F) #Sample.init error if you include MC
  Dataset2 <- NextRound[partition,]; Validation_Data <- NextRound[-partition,]
  rm(NextRound)
  
  ##################################################################################
  ######Round 1
  ###Creating Round 1 Model
  Equation <- as.formula(paste(df_Response_names[i], " ~", paste(Predictor_variable_names)))
  start_train <- Sys.time()
  Round1_Model <- train( 
    data = Dataset1,
    Equation,                 #Equation needs to be blank term I think?
    method = "parRF",         #parRF - parallel random forest
    ntree = num_trees,
    importance=T,             #Ranks order of importance for predictor variables
    tuneGrid = tune_grid,     #Parameters to tune/optimize
    trControl = fitControl)
  end_train <- Sys.time(); capture.output(end_train-start_train)
  Round1_Model; confusionMatrix(Round1_Model)
  #Save location of Round 1 model
  saveRDS(Round1_Model, file = paste0(wd, "/Model_Comparison/models/Round1/", time, "_", Plain_Response_names[i], ".rds"))
  
  ###Calculate predicted values of round 2 points based on round 1 model
  D2_uncertainty <- Dataset2[c("ID", names(Covariates))]
  D2_uncertainty <- stats::predict(Round1_Model, D2_uncertainty, type = "prob")
  
  ###Calculate uncertainty metrics
  source("./[subcode] - Calculate uncertainty metrics.R")
  
  ###Calculate the entropy distribution for Dataset 2
  #If file exists, add to it; if not, create it.
  if(file.exists(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"))){D2_unc_dist <- read.csv(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), stringsAsFactors = F)
    D2_unc_dist <- rbind(D2_unc_dist, data.frame(time = time, ID = Dataset2$ID, Variable = Plain_Response_names[i], IU = Dataset2$IU, EU = Dataset2$EU, CI = Dataset2$CI))
    write.csv(D2_unc_dist, paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), row.names = F)
  }else{
    write.csv(data.frame(time = time, ID = Dataset2$ID, Variable = Plain_Response_names[i], IU = Dataset2$IU, EU = Dataset2$EU, CI = Dataset2$CI), paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), row.names = F)
    }
  rm(D2_uncertainty)
  
  
  ##################################################################################
  ######Round 2
  ###Similarity in high uncertainty points
  #Point threshold
  if(file.exists(paste0(wd, "/D2_unc_similarity/D2_uncertainty_similarity.csv"))){D2_unc_sim <- read.csv(paste0(wd, "/D2_unc_similarity/D2_uncertainty_similarity.csv"), stringsAsFactors = F)
  D2_unc_sim <- rbind(D2_unc_sim, data.frame(time = time,
                                             Variable = Plain_Response_names[i],
                                             "IU_EU_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,],Dataset2[order(-Dataset2$EU),][1:Cutoff,])), #IU-EU
                                             "IU_CI_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,],Dataset2[order(-Dataset2$CI),][1:Cutoff,])), #IU-CI
                                             "EU_CI_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$EU),][1:Cutoff,],Dataset2[order(-Dataset2$CI),][1:Cutoff,])), #EU-CI
                                             "IU_EU_CI_Similarity" = nrow(dplyr::intersect(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,],Dataset2[order(-Dataset2$EU),][1:Cutoff,]), Dataset2[order(-Dataset2$CI),][1:Cutoff,])) #IU_EU_CI
                                             )
                      )
  write.csv(D2_unc_sim, paste0(wd, "/D2_unc_similarity/D2_uncertainty_similarity.csv"), row.names = F)
  }else{write.csv(data.frame(time = time,
                             Variable = Plain_Response_names[i],
                             "IU_EU_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,],Dataset2[order(-Dataset2$EU),][1:Cutoff,])), #IU-EU
                             "IU_CI_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,],Dataset2[order(-Dataset2$CI),][1:Cutoff,])), #IU-CI
                             "EU_CI_Similarity" = nrow(dplyr::intersect(Dataset2[order(-Dataset2$EU),][1:Cutoff,],Dataset2[order(-Dataset2$CI),][1:Cutoff,])), #EU-CI
                             "IU_EU_CI_Similarity" = nrow(dplyr::intersect(dplyr::intersect(Dataset2[order(-Dataset2$IU),][1:Cutoff,], Dataset2[order(-Dataset2$EU),][1:Cutoff,]), Dataset2[order(-Dataset2$CI),][1:Cutoff,]))),
                  paste0(wd, "/D2_unc_similarity/D2_uncertainty_similarity.csv"), row.names = F)
    }
  
  ###############
  ###############
  #ROUND 2 Models
  
  ##########
  #RANDOM TREATMENT
  #Reducing to smallest sample size and creating Random uncertainty dataset - can do this before treatment loop for point threshold but not for uncertainty threshold
  DSub_RandomUnc <- Dataset2[sample(1:nrow(Dataset2), Cutoff, replace = F),]
  
  Dataset3_Random <- base::merge(Dataset1, DSub_RandomUnc, all=T)
  Dataset3_Random <- Dataset3_Random[,!(names(Dataset3_Random) %in% c("IU", "EU", "CI"))]
  
  start_train <- Sys.time()
  Round2_Random_Model <- train( 
    data = Dataset3_Random,
    Equation,                 #Equation needs to be blank term I think?
    method = "parRF",         #parRF - parallel random forest
    ntree = num_trees,
    importance=T,             #Ranks order of importance for predictor variables
    tuneGrid = tune_grid,     #Parameters to tune/optimize
    trControl = fitControl)
  end_train <- Sys.time(); capture.output(end_train-start_train)
  Round2_Random_Model; confusionMatrix(Round2_Random_Model)
  #Save location of Round 2 model
  saveRDS(Round2_Random_Model, file = paste0(wd, "/Model_Comparison/models/Round2/", time, "_", "Random", "_", Plain_Response_names[i], ".rds"))
  print(paste0("Done ", Plain_Response_names[i], " for Random treatment, repeat ", repeat_count))
  
  ##########
  #HIGH UNCERTAINTY TREATMENT
  ###Looping for different uncertainty metric selection
  Uncertainty_metrics <- c("IU", "EU", "CI")
  
  for(Unc_metric in 1:length(Uncertainty_metrics)){
  
  DSub_HighUnc <- Dataset2[order(-Dataset2[Uncertainty_metrics[Unc_metric]]),]
  DSub_HighUnc$Unc_count <- 1:nrow(DSub_HighUnc)
  DSub_HighUnc <- DSub_HighUnc[(DSub_HighUnc$Unc_count <= Cutoff),]
  DSub_HighUnc$Unc_count <- NULL
  rownames(DSub_HighUnc) <- 1:nrow(DSub_HighUnc)
  
  #Record uncertainty distribution of Dataset 2 Subselection for each uncertainty metric
  if(file.exists(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metrics[Unc_metric], ".csv"))){DSub_Unc <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metrics[Unc_metric], ".csv"), stringsAsFactors = F)
  #Need extra square brackets around sort to turn from dataframe to vector
  DSub_Unc <- rbind(DSub_Unc,
                        data.frame(time = time, 
                                   Variable = Plain_Response_names[i],
                                   Uncertainty_metric = Uncertainty_metrics[Unc_metric],
                                   DSub_Unc_threshold = min(DSub_HighUnc[Uncertainty_metrics[Unc_metric]]), 
                                   DSub_HighUncValues = sort(DSub_HighUnc[[Uncertainty_metrics[Unc_metric]]]), 
                                   DSub_RandomUncValues = sort(DSub_RandomUnc[[Uncertainty_metrics[Unc_metric]]])
                                   )
                    )
  write.csv(DSub_Unc, paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metrics[Unc_metric], ".csv"), row.names = F)
  
  }else{write.csv(data.frame(time = time, 
                             Variable = Plain_Response_names[i],
                             Uncertainty_metric = Uncertainty_metrics[Unc_metric],
                             DSub_Unc_threshold = min(DSub_HighUnc[Uncertainty_metrics[Unc_metric]]), 
                             DSub_HighUncValues = sort(DSub_HighUnc[[Uncertainty_metrics[Unc_metric]]]), 
                             DSub_RandomUncValues = sort(DSub_RandomUnc[[Uncertainty_metrics[Unc_metric]]])
                             ),
                  paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Uncertainty_metrics[Unc_metric], ".csv"), row.names = F)
    }
  
  ##################################################################################
  ##################################################################################
  
  ##########
  ##########
  #Covariate Uncertainty analysis 
  ####Quantify soil and environment value covariation with uncertainty
  source("./[subcode] - Quantify soil and environment value covariation with uncertainty - Same random.R")
  ##########
  ##########
  
  ###Creating Round 2 Model
  Dataset3_High <- base::merge(Dataset1, DSub_HighUnc, all=T)
  Dataset3_High <- Dataset3_High[,!(names(Dataset3_High) %in% c("IU", "EU", "CI"))]
  
  start_train <- Sys.time()
  Round2_High_Model <- train( 
    data = Dataset3_High,
    Equation,                 #Equation needs to be blank term I think?
    method = "parRF",         #parRF - parallel random forest
    ntree = num_trees,
    importance=T,             #Ranks order of importance for predictor variables
    tuneGrid = tune_grid,     #Parameters to tune/optimize
    trControl = fitControl)
  end_train <- Sys.time(); capture.output(end_train-start_train)
  Round2_High_Model; confusionMatrix(Round2_High_Model)
  #Save location of Round 2 model
  saveRDS(Round2_High_Model, file = paste0(wd, "/Model_Comparison/models/Round2/", time, "_", Uncertainty_metrics[Unc_metric], "_", "High", "_", Plain_Response_names[i], ".rds"))
  
  
  
  ##################################################################################
  ######Validation data
  ###Predict on withheld validation points

  #It's unnecessary to calculate R1 and R2 Random treatment predictions multiple times per uncertainty metric and due to stochasticity in predict function results will vary slightly
  if(Unc_metric==1){
  #Round 1 - Only calculate Round 1 results on High Uncertainty treatment. 
  Validation_Data$Model_Predictions_R1 <- predict(Round1_Model, newdata = Validation_Data, type = "raw", na.action=na.exclude)
  Round1_conMat <- confusionMatrix(Validation_Data$Model_Predictions_R1, Validation_Data[,loop_Name])
  
  #Round 2 - Random
  Validation_Data$Model_Predictions_R2_Random <- predict(Round2_Random_Model, newdata = Validation_Data, type = "raw", na.action=na.exclude)
  Round2_Random_conMat <- confusionMatrix(Validation_Data$Model_Predictions_R2_Random, Validation_Data[,loop_Name])

  if(file.exists(paste0(wd, "/Model_Comparison/Model_Comparison.csv"))){
    Model_Comparison <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = FALSE)
    Model_Comparison <- rbind(Model_Comparison,
                              data.frame(Time = time,
                                         Uncertainty_metric = "Random",
                                         Uncertainty = "Random",
                                         Variable = Plain_Response_names[i],
                                         Datapoints_added_to_dataset_1 = nrow(DSub_RandomUnc),
                                         Round1_Accuracy = Round1_conMat$overall[["Accuracy"]],
                                         Round2_Accuracy = Round2_Random_conMat$overall[["Accuracy"]],
                                         Accuracy_dif = Round2_Random_conMat$overall[["Accuracy"]] - Round1_conMat$overall[["Accuracy"]],
                                         Round1_Kappa = Round1_conMat$overall[["Kappa"]],
                                         Round2_Kappa = Round2_Random_conMat$overall[["Kappa"]],
                                         Kappa_dif = Round2_Random_conMat$overall[["Kappa"]] - Round1_conMat$overall[["Kappa"]],
                                         stringsAsFactors=F)
                              )
    write.csv(Model_Comparison, paste0(wd, "/Model_Comparison/Model_Comparison.csv"), row.names = F)
    }else{Model_Comparison <- data.frame(Time = time,
                                         Uncertainty_metric = "Random",
                                         Uncertainty = "Random",
                                         Variable = Plain_Response_names[i],
                                         Datapoints_added_to_dataset_1 = nrow(DSub_RandomUnc),
                                         Round1_Accuracy = Round1_conMat$overall[["Accuracy"]],
                                         Round2_Accuracy = Round2_Random_conMat$overall[["Accuracy"]],
                                         Accuracy_dif = Round2_Random_conMat$overall[["Accuracy"]] - Round1_conMat$overall[["Accuracy"]],
                                         Round1_Kappa = Round1_conMat$overall[["Kappa"]],
                                         Round2_Kappa = Round2_Random_conMat$overall[["Kappa"]],
                                         Kappa_dif = Round2_Random_conMat$overall[["Kappa"]] - Round1_conMat$overall[["Kappa"]],
                                         stringsAsFactors=F)
    write.csv(Model_Comparison, paste0(wd, "/Model_Comparison/Model_Comparison.csv"), row.names = F)
    }
  }
  
  #Round 2 - High
  Validation_Data$Model_Predictions_R2_High <- predict(Round2_High_Model, newdata = Validation_Data, type = "raw", na.action=na.exclude)
  Round2_High_conMat <- confusionMatrix(Validation_Data$Model_Predictions_R2_High, Validation_Data[,loop_Name])
  
  ####write out results to csv
  if(!dir.exists(paste0(wd, "/Model_Comparison/models/csv"))){dir.create(paste0(wd, "/Model_Comparison/models/csv"))}
  write.csv(Validation_Data, paste0(wd, "/Model_Comparison/models/csv/", time, "_", Uncertainty_metrics[Unc_metric], "_", "High", "_", Plain_Response_names[i], ".csv"), row.names = F)
  
  ######Create summary of multiple run accuracies
  Model_Comparison <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = FALSE)
  Model_Comparison <- rbind(Model_Comparison,
                            data.frame(Time = time,
                                       Uncertainty_metric = Uncertainty_metrics[Unc_metric],
                                       Uncertainty = "High",
                                       Variable = Plain_Response_names[i],
                                       Datapoints_added_to_dataset_1 = nrow(DSub_RandomUnc),
                                       Round1_Accuracy = Round1_conMat$overall[["Accuracy"]],
                                       Round2_Accuracy = Round2_High_conMat$overall[["Accuracy"]],
                                       Accuracy_dif = Round2_High_conMat$overall[["Accuracy"]] - Round1_conMat$overall[["Accuracy"]],
                                       Round1_Kappa = Round1_conMat$overall[["Kappa"]],
                                       Round2_Kappa = Round2_High_conMat$overall[["Kappa"]],
                                       Kappa_dif = Round2_High_conMat$overall[["Kappa"]] - Round1_conMat$overall[["Kappa"]],
                                       stringsAsFactors=F)
                            )
  write.csv(Model_Comparison, paste0(wd, "/Model_Comparison/Model_Comparison.csv"), row.names = F)
    
  print(paste0("Done ", Plain_Response_names[i], " for ", Uncertainty_metrics[Unc_metric], " High treatment, repeat ", repeat_count))
  }#Close Uncertainty metric loop
  
  } #Close soil response variable loop
  ##################################################################################
  ##################################################################################
  ##################################################################################
  ##################################################################################
  
  if(repeats == repeat_count){break}
} #Close repeat simulation loop


#Write out file summarizing uncertainty model parameters

#Loading in uncertainty distribution files for each uncertainty metric
DSub_Unc_Dis_IU <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_IU.csv"))
DSub_Unc_Dis_EU <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_EU.csv"))
DSub_Unc_Dis_CI <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_CI.csv"))

df = data.frame(Number_of_simulations = length(unique(Model_Comparison$Time)), #This is number of simulations for each uncertainty metric-treatment-soil attribute (1 simulation = 2 soil attributes * 3 uncertainty metrics * 2 treatments = 12 rows)
                Dataset1_Size = Dataset1_size,
                Dataset2_initSize = Dataset2_initsize,
                Cutoff = Cutoff,
                
                Average_DSub_Unc_IU_M = mean(unique(DSub_Unc_Dis_IU[DSub_Unc_Dis_IU$Variable == "M",]$DSub_Unc_threshold)),
                Average_DSub_Unc_IU_T = mean(unique(DSub_Unc_Dis_IU[DSub_Unc_Dis_IU$Variable == "T",]$DSub_Unc_threshold)),
                Var_DSub_Unc_IU_M = var(unique(DSub_Unc_Dis_IU[DSub_Unc_Dis_IU$Variable == "M",]$DSub_Unc_threshold)),
                Var_DSub_Unc_IU_T = var(unique(DSub_Unc_Dis_IU[DSub_Unc_Dis_IU$Variable == "T",]$DSub_Unc_threshold)),
                
                Average_DSub_Unc_EU_M = mean(unique(DSub_Unc_Dis_EU[DSub_Unc_Dis_EU$Variable == "M",]$DSub_Unc_threshold)),
                Average_DSub_Unc_EU_T = mean(unique(DSub_Unc_Dis_EU[DSub_Unc_Dis_EU$Variable == "T",]$DSub_Unc_threshold)),
                Var_DSub_Unc_EU_M = var(unique(DSub_Unc_Dis_EU[DSub_Unc_Dis_EU$Variable == "M",]$DSub_Unc_threshold)),
                Var_DSub_Unc_EU_T = var(unique(DSub_Unc_Dis_EU[DSub_Unc_Dis_EU$Variable == "T",]$DSub_Unc_threshold)),
                
                Average_DSub_Unc_CI_M = mean(unique(DSub_Unc_Dis_CI[DSub_Unc_Dis_CI$Variable == "M",]$DSub_Unc_threshold)),
                Average_DSub_Unc_CI_T = mean(unique(DSub_Unc_Dis_CI[DSub_Unc_Dis_CI$Variable == "T",]$DSub_Unc_threshold)),
                Var_DSub_Unc_CI_M = var(unique(DSub_Unc_Dis_CI[DSub_Unc_Dis_CI$Variable == "M",]$DSub_Unc_threshold)),
                Var_DSub_Unc_CI_T = var(unique(DSub_Unc_Dis_CI[DSub_Unc_Dis_CI$Variable == "T",]$DSub_Unc_threshold))
                )

write.csv(df, paste0(wd, "/Model_Summary.csv"), row.names = F)

end_time <- Sys.time()
print(paste("End of Uncertainty Models.R ", capture.output(end_time - start_time)))
#####
####
###
##
#END
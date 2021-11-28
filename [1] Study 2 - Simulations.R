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

repeats = 7
repeat_count = 0

Dataset1_size = 300
Dataset2_initsize = 6000

Dataset2_subset = c(300,600,900,1200)

Uncertainty_metric_abbrev = "CI" #Option are "IU" - Ignorance Uncertainty; "EU" - Exaggeration Uncertainty; "CI" - Confidence Index

#####
#Create directories for response dataset
if(!dir.exists(paste0("./Uncertainty_Models/"))){dir.create(paste0("./Uncertainty_Models/"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data"))}

#Create directories for model output depending on if you are using point of entropy threshold approach
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2/", Dataset1_size))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2/", Dataset1_size))}
#Short version of working directory which is useful for later
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study2/", Dataset1_size)

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
if(length(names(Covariates)) %% 2 == 1){mtry = seq(from = 3, to = length(names(Covariates)), by = 2)
}else{mtry = seq(from = 2, to = length(names(Covariates)), by = 2)} #if even
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
  
  ##########
  ###Calculate uncertainty metric - just choose the best performing one
  if(Uncertainty_metric_abbrev == "IU"){
  ###Ignorance uncertainty - Ranges from 0(high certainty) to 1 (low certainty)
  f1 <- function(x){x*log(x+0.000001)}
  Ignorance_Uncertainty <- as.data.frame(apply(D2_uncertainty, MARGIN = c(1,2), FUN = f1))
  Ignorance_Uncertainty$IU <- apply(Ignorance_Uncertainty, MARGIN = 1, FUN = sum)
  Ignorance_Uncertainty$IU <- -1*Ignorance_Uncertainty$IU/log(nlevels(Round1_Model))
  Ignorance_Uncertainty$IU[Ignorance_Uncertainty$IU < 0] <- 0
  Ignorance_Uncertainty <- Ignorance_Uncertainty[("IU")]
  Uncertainty_calculation <- Ignorance_Uncertainty
  }
  
  if(Uncertainty_metric_abbrev == "EU"){
  ###Exaggeration uncertainty
  Exaggeration_Uncertainty <- data.frame(EU = apply(D2_uncertainty, MARGIN = 1, FUN = max))
  Exaggeration_Uncertainty$EU <- 1 - Exaggeration_Uncertainty$EU
  Uncertainty_calculation <- Exaggeration_Uncertainty
  }
  
  if(Uncertainty_metric_abbrev == "CI"){
  ###Confusion index
  f2 <- function(x){1 - (max(x) - max(x[x != max(x)]))}
  Confusion_Index <- data.frame(CI = apply(D2_uncertainty, MARGIN = 1, FUN = f2))
  Uncertainty_calculation <- Confusion_Index
  }
  ##########
  
  #Combine metrics
  Dataset2 <- merge(Dataset2, Uncertainty_calculation, by = "row.names"); Dataset2$Row.names = NULL
  
  ###Calculate the entropy distribution for Dataset 2
  #If file exists, add to it; if not, create it.
  if(file.exists(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"))){D2_unc_dist <- read.csv(paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), stringsAsFactors = F)
    D2_unc_dist <- rbind(D2_unc_dist, data.frame(time = time, ID = Dataset2$ID, Variable = Plain_Response_names[i], Uncertainty_metric = Uncertainty_calculation))
    write.csv(D2_unc_dist, paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), row.names = F)
  }else{
    write.csv(data.frame(time = time, ID = Dataset2$ID, Variable = Plain_Response_names[i], Uncertainty_metric = Uncertainty_calculation), paste0(wd, "/D2_unc_distribution/D2_uncertainty_dist.csv"), row.names = F)
    }
  rm(D2_uncertainty)
  
  
  ##################################################################################
  ######Round 2

  ###############
  ###############
  ###Looping for different uncertainty metric selection
  
  #Reducing to smallest sample size and creating Random uncertainty dataset - can do this before treatment loop for point threshold but not for uncertainty threshold
  for(treatment in 1:length(Dataset2_subset)){
  DSub_RandomUnc <- Dataset2[sample(1:nrow(Dataset2), Dataset2_subset[treatment], replace = F),]
  
  DSub_HighUnc <- Dataset2[order(-Dataset2[Uncertainty_metric_abbrev]),]
  DSub_HighUnc$Unc_count <- 1:nrow(DSub_HighUnc)
  DSub_HighUnc <- DSub_HighUnc[(DSub_HighUnc$Unc_count <= Dataset2_subset[treatment]),]
  DSub_HighUnc$Unc_count <- NULL
  
  #for(Unc_metric in 1:length(Uncertainty_metrics)){}
  #Record size of Dataset 2 Subselection
  if(file.exists(paste0(wd, "/D2Subselection/DSub_size.csv"))){DSub_size <- read.csv(paste0(wd, "/D2Subselection/DSub_size.csv"), stringsAsFactors = F)
    DSub_size <- rbind(DSub_size, data.frame(time = time, Variable = Plain_Response_names[i], Uncertainty_metric = Uncertainty_metric_abbrev, DSub_size = nrow(DSub_HighUnc)))
    write.csv(DSub_size, paste0(wd, "/D2Subselection/DSub_size.csv"), row.names = F)
  }else{write.csv(data.frame(time = time, Variable = Plain_Response_names[i], Uncertainty_metric = Uncertainty_metric_abbrev, DSub_size = nrow(DSub_HighUnc)), 
                  paste0(wd, "/D2Subselection/DSub_size.csv"), row.names = F)
    }
  
  #Record uncertainty distribution of Dataset 2 Subselection for each uncertainty metric
  if(file.exists(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Dataset2_subset[treatment], ".csv"))){DSub_Unc <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Dataset2_subset[treatment], ".csv"), stringsAsFactors = F)
  #Need extra square brackets around sort to turn from dataframe to vector
  DSub_Unc <- rbind(DSub_Unc,
                        data.frame(time = time, 
                                   Variable = Plain_Response_names[i],
                                   Uncertainty_metric = Uncertainty_metric_abbrev,
                                   DSub_Unc_threshold = min(DSub_HighUnc[Uncertainty_metric_abbrev]), 
                                   DSub_HighUncValues = sort(DSub_HighUnc[[Uncertainty_metric_abbrev]]), 
                                   DSub_RandomUncValues = sort(DSub_RandomUnc[[Uncertainty_metric_abbrev]])
                                   )
                    )
  write.csv(DSub_Unc, paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Dataset2_subset[treatment], ".csv"), row.names = F)
  
  }else{write.csv(data.frame(time = time, 
                             Variable = Plain_Response_names[i],
                             Uncertainty_metric = Uncertainty_metric_abbrev,
                             DSub_Unc_threshold = min(DSub_HighUnc[Uncertainty_metric_abbrev]), 
                             DSub_HighUncValues = sort(DSub_HighUnc[[Uncertainty_metric_abbrev]]), 
                             DSub_RandomUncValues = sort(DSub_RandomUnc[[Uncertainty_metric_abbrev]])
                             ),
                  paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Dataset2_subset[treatment], ".csv"), row.names = F)
    }
  
  ##################################################################################
  ##################################################################################
  DSub_uncertainties <- list(DSub_RandomUnc, DSub_HighUnc)
  DSub_uncertainties_names <- c("RANDOM", "HIGH")

  for(j in 1:length(DSub_uncertainties)){
  
  DSub_selection <- DSub_uncertainties[[j]]
  rownames(DSub_selection) <- 1:nrow(DSub_selection)

  ##########
  ##########
  #Covariate Uncertainty analysis 
  ####Quantify soil and environment value covariation with uncertainty
  
  #source("./[subcode] - Quantify soil and environment value covariation with uncertainty.R")
  
  ##########
  ##########
  
  ###Creating Round 2 Model
  Dataset3 <- base::merge(Dataset1, DSub_selection, all=T)
  Dataset3 <- Dataset3[,!(names(Dataset3) %in% Uncertainty_metric_abbrev)]
  
  start_train <- Sys.time()
  Round2_Model <- train( 
    data = Dataset3,
    Equation,                 #Equation needs to be blank term I think?
    method = "parRF",         #parRF - parallel random forest
    ntree = num_trees,
    importance=T,             #Ranks order of importance for predictor variables
    tuneGrid = tune_grid,     #Parameters to tune/optimize
    trControl = fitControl)
  end_train <- Sys.time(); capture.output(end_train-start_train)
  Round2_Model; confusionMatrix(Round2_Model)
  #Save location of Round 2 model
  saveRDS(Round2_Model, file = paste0(wd, "/Model_Comparison/models/Round2/", time, "_", Uncertainty_metric_abbrev, "_",  Dataset2_subset[treatment], "_", DSub_uncertainties_names[j], "_", Plain_Response_names[i], ".rds"))
  
  ##################################################################################
  ######Validation data
  ###Predict on withheld validation points
  #Round 1 - Only calculate Round 1 results on High Uncertainty treatment. 
  #It's unnecessary to do so twice (i.e. for Random Uncertainty as well) and do to stochasticity in predict function results will vary slightly
  if(j==1){
  Validation_Data$Model_Predictions_R1 <- predict(Round1_Model, newdata = Validation_Data, type = "raw", na.action=na.exclude)
  Round1_conMat <- confusionMatrix(Validation_Data$Model_Predictions_R1, Validation_Data[,loop_Name])
  }
  #Round 2
  Validation_Data$Model_Predictions_R2 <- predict(Round2_Model, newdata = Validation_Data, type = "raw", na.action=na.exclude)
  Round2_conMat <- confusionMatrix(Validation_Data$Model_Predictions_R2, Validation_Data[,loop_Name])
  
  ####write out results to csv
  if(!dir.exists(paste0(wd, "/Model_Comparison/models/csv"))){dir.create(paste0(wd, "/Model_Comparison/models/csv"))}
  write.csv(Validation_Data, paste0(wd, "/Model_Comparison/models/csv/", time, "_", Uncertainty_metric_abbrev, "_", Dataset2_subset[treatment], "_", DSub_uncertainties_names[j], "_", Plain_Response_names[i], ".csv"), row.names = F)
  
  ######Create summary of multiple run accuracies
  if(file.exists(paste0(wd, "/Model_Comparison/Model_Comparison.csv"))){
    Model_Comparison <- read.csv(paste0(wd, "/Model_Comparison/Model_Comparison.csv"), stringsAsFactors = FALSE)
    Model_Comparison <- rbind(Model_Comparison,
                              data.frame(Time = time,
                                         Uncertainty_metric = Uncertainty_metric_abbrev,
                                         Uncertainty = DSub_uncertainties_names[j],
                                         Data_pts_added = Dataset2_subset[treatment],
                                         Variable = Plain_Response_names[i],
                                         Datapoints_added_to_dataset_1 = nrow(DSub_HighUnc),
                                         Round1_Accuracy = Round1_conMat$overall[["Accuracy"]],
                                         Round2_Accuracy = Round2_conMat$overall[["Accuracy"]],
                                         Accuracy_dif = Round2_conMat$overall[["Accuracy"]] - Round1_conMat$overall[["Accuracy"]],
                                         Round1_Kappa = Round1_conMat$overall[["Kappa"]],
                                         Round2_Kappa = Round2_conMat$overall[["Kappa"]],
                                         Kappa_dif = Round2_conMat$overall[["Kappa"]] - Round1_conMat$overall[["Kappa"]], 
                                         stringsAsFactors=F)
                              )
    write.csv(Model_Comparison, paste0(wd, "/Model_Comparison/Model_Comparison.csv"), row.names = F)
    
    }else{Model_Comparison <- data.frame(Time = time,
                                         Uncertainty_metric = Uncertainty_metric_abbrev,
                                         Uncertainty = DSub_uncertainties_names[j],
                                         Data_pts_added = Dataset2_subset[treatment],
                                         Variable = Plain_Response_names[i],
                                         Datapoints_added_to_dataset_1 = nrow(DSub_HighUnc),
                                         Round1_Accuracy = Round1_conMat$overall[["Accuracy"]],
                                         Round2_Accuracy = Round2_conMat$overall[["Accuracy"]],
                                         Accuracy_dif = Round2_conMat$overall[["Accuracy"]] - Round1_conMat$overall[["Accuracy"]],
                                         Round1_Kappa = Round1_conMat$overall[["Kappa"]],
                                         Round2_Kappa = Round2_conMat$overall[["Kappa"]],
                                         Kappa_dif = Round2_conMat$overall[["Kappa"]] - Round1_conMat$overall[["Kappa"]],
                                         stringsAsFactors=F)
    #Adding rows to a dataframe is a pain
    write.csv(Model_Comparison, paste0(wd, "/Model_Comparison/Model_Comparison.csv"), row.names = F)
    }
  
  print(paste0("Done ", Plain_Response_names[i], " for ", DSub_uncertainties_names[j], " ", Dataset2_subset[treatment], "pts treatment, repeat ", repeat_count))
  }#Close Uncertainty treatment loop
  
  }#Close Data pts loop
  
  } #Close soil response variable loop
  ##################################################################################
  ##################################################################################
  ##################################################################################
  ##################################################################################
  
  if(repeats == repeat_count){break}
} #Close repeat simulation loop


#Write out file summarizing uncertainty model parameters
#Loading in uncertainty distribution files for each uncertainty metric
df_combine <- vector()

for(index in 1:length(Dataset2_subset)){
  temp_df <- read.csv(paste0(wd, "/D2Subselection/DSub_Unc_Dis_", Dataset2_subset[index], ".csv"))
  
  temp_df <- data.frame(mean(unique(temp_df[temp_df$Variable == "M",]$DSub_Unc_threshold)),
                          mean(unique(temp_df[temp_df$Variable == "T",]$DSub_Unc_threshold)),
                          var(unique(temp_df[temp_df$Variable == "M",]$DSub_Unc_threshold)),
                          var(unique(temp_df[temp_df$Variable == "T",]$DSub_Unc_threshold))
                          )
  
  colnames(temp_df) <- c(paste0("Average_DSub_Unc_M_", Dataset2_subset[index]),
                           paste0("Average_DSub_Unc_T_", Dataset2_subset[index]),
                           paste0("Var_DSub_Unc_M_", Dataset2_subset[index]),
                           paste0("Var_DSub_Unc_T_", Dataset2_subset[index])
                           )
  
  df_combine = append(df_combine, temp_df)
  }

df_combine <- data.frame(df_combine)

df <- data.frame(Number_of_simulations = length(unique(Model_Comparison$Time)), #This is number of simulations for each uncertainty metric-treatment-soil attribute (1 simulation = 2 soil attributes * 3 uncertainty metrics * 2 treatments = 12 rows)
                Dataset1_Size = Dataset1_size,
                Dataset2_initSize = Dataset2_initsize,
                Treatments = paste0(Dataset2_subset, collapse = " ")
                )

df <- cbind(df, df_combine)

write.csv(df, paste0(wd, "/Model_Summary.csv"), row.names = F)

end_time <- Sys.time()
print(paste("End of Uncertainty Models.R ", capture.output(end_time - start_time)))
#####
####
###
##
#END
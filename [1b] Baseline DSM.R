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
Create_prob_raster = TRUE

#####
#Create directories for response dataset
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers/Raw"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers/Raw"))}
if(!dir.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers/Probability"))){dir.create(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/layers/Probability"))}

#Short version of working directory which is useful for later
wd = paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Baseline/")

#Running in parallel
detectCores() #Should be 12 on SOIL-BOT
beginCluster(detectCores()-1)
#How to stop cluster?
#stopImplicitCluster() or #endCluster(cl) or #stopCluster(cl)

###Acquiring continuous environmental covariates
Topo_Covariates <- list.files(path= paste0("./Predictor_Variables/Topography_metrics/", Raw_DEM_reso, "DEM", Filter_type), pattern="\\.tif$", full.names = TRUE)
Hydro_Covariates <- list.files(path= paste0("./Predictor_Variables/Hydrology_metrics/", Raw_DEM_reso, "DEM"), pattern="\\.tif$", full.names = TRUE)
Bio_continuousCovariates <- list.files(path= paste0("./Predictor_Variables/Biota_metrics/", Raw_DEM_reso, "DEM/continuous"), pattern="\\.tif$", full.names = TRUE)
Dist_Covariates <- list.files(path= paste0("./Predictor_Variables/DistanceFields/", Raw_DEM_reso, "DEM"), pattern="\\.tif$", full.names = TRUE)
###Acquiring categorical environmental covariates
Geo_Covariates <- list.files(path= paste0("./Predictor_Variables/Geology_metrics/", Raw_DEM_reso, "DEM"), pattern="\\.tif$", full.names = TRUE)
Bio_categoricalCovariates <- list.files(path= paste0("./Predictor_Variables/Biota_metrics/", Raw_DEM_reso, "DEM/categorical"), pattern="\\.tif$", full.names = TRUE)

#Stacking covariates
Covariates <- c(Topo_Covariates, Hydro_Covariates, Bio_continuousCovariates, Dist_Covariates, Geo_Covariates, Bio_categoricalCovariates)
Covariates <- raster::stack(Covariates)
rm(Bio_categoricalCovariates, Bio_continuousCovariates, Geo_Covariates, Hydro_Covariates, Dist_Covariates) #Clean up R global environment

df <- read.csv(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data/", Raw_DEM_reso, "DEM", Filter_type,"Response_df.csv"))
df$X <- as.numeric(df$X); df$Y <- as.numeric(df$Y)
Predictor_variable_names <- paste(names(Covariates), collapse="+")

#Making legitimate Moisture regime class
df$MR <- as.character(df$MOISTURE_R)

df$MR <- gsub(x = df$MR, pattern = "x", replacement = "Xeric")
df$MR <- gsub(x = df$MR, pattern = "s", replacement = "Saturated")
df$MR <- gsub(x = df$MR, pattern = "-1", replacement = "Dry") 
df$MR <- gsub(x = df$MR, pattern = "0", replacement = "Moderately_Dry")
df$MR <- gsub(x = df$MR, pattern = "1", replacement = "Moderately_Fresh")
df$MR <- gsub(x = df$MR, pattern = "2", replacement = "Fresh")
df$MR <- gsub(x = df$MR, pattern = "3", replacement = "Very_Fresh")
df$MR <- gsub(x = df$MR, pattern = "4", replacement = "Moderately_Moist")
df$MR <- gsub(x = df$MR, pattern = "5", replacement = "Moist")
df$MR <- gsub(x = df$MR, pattern = "6", replacement = "Very_Moist")
df$MR <- gsub(x = df$MR, pattern = "7", replacement = "Moderately_Wet")
df$MR <- gsub(x = df$MR, pattern = "8", replacement = "Wet")
df$MR <- gsub(x = df$MR, pattern = "9", replacement = "Very_Wet")

df$MR <- as.factor(df$MR)

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

####Looping over predictor dataset
Plain_Response_names <- c("M", "T")
df_Response_names <- c("MR", "TEXT")
Random_Forest_results <- data.frame(NULL)

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

  ###Creating RF Model
  Equation <- as.formula(paste(df_Response_names[i], " ~", paste(Predictor_variable_names)))
  start_train <- Sys.time()
  RF_Train <- train( 
    data = Train_df,
    Equation,                 #Equation needs to be blank term I think?
    method = "parRF",         #parRF - parallel random forest
    ntree = num_trees,
    importance=T,             #Ranks order of importance for predictor variables
    tuneGrid = tune_grid,     #Parameters to tune/optimize
    trControl = fitControl)
  end_train <- Sys.time(); capture.output(end_train-start_train)
  #Save location of RF model
  saveRDS(RF_Train, file = paste0(wd, Plain_Response_names[i], "_RF.rds"))
  
  #Saving accuracy output
  if(nrow(RF_Train$results[RF_Train$results$Accuracy == max(RF_Train$results$Accuracy),]) > 1){print(paste0("Multiple solutions for RF ", Plain_Response_names[i]))}
  
  #params_score <- data.frame(RF_Train$bestTune)
  Accuracy_score <- RF_Train$results[rownames(RF_Train$results) == rownames(RF_Train$bestTune),]
  Accuracy_score$Variable <- Plain_Response_names[i]
  Random_Forest_results <- rbind(Random_Forest_results, Accuracy_score)
  
  #^This takes a long time even with multiple cores
  print(paste0("Done ", Plain_Response_names[i], " for Random Forest"))
  
}

#Write out csv
Random_Forest_results$Model <- "RF"
Random_Forest_results$ntree <- num_trees
write.csv(Random_Forest_results, paste0(wd, "/RF_Model_performances.csv"), row.names = F)

#Write out predictive raster
Soil_Variable <- c("M", "T")

M_legend_order <- c("Xeric", "Dry", "Moderately_Dry", "Moderately_Fresh", "Fresh", 
                    "Very_Fresh", "Moderately_Moist", "Moist", "Very_Moist", 
                    "Moderately_Wet", "Wet", "Very_Wet", "Saturated")

T_legend_order <- c("Organic", "Silty_Clay", "Silt", "Silt_Loam",
                    "Silty_Clay_Loam", "Loam", "Clay_Loam", 
                    "Sandy_Clay_Loam", "Silty_Sand", 
                    "Sandy_Loam","Loamy_Sand", "Sand")

for (i in 1:length(Plain_Response_names)){
  RF_Train <- readRDS(paste0(wd, Plain_Response_names[i], "_RF.rds"))
  Map <- clusterR(Covariates, fun = raster::predict, args = list(RF_Train, type = "raw"))
  
  #Temporary reclass to organize for ordinal reclass
  reclass_matrix <- matrix(c(0.5:(length(levels(RF_Train$pred$obs))-0.5), 1.5:(length(levels(RF_Train$pred$obs))+0.5), 101:(length(levels(RF_Train$pred$obs))+100)),
                           nrow = length(levels(RF_Train$pred$obs)), ncol = 3)
  Map_reclass <- reclassify(Map, reclass_matrix)
  
  Current_classification <- data.frame(Current_ID = reclass_matrix[,3], Attribute = levels(RF_Train$pred$obs))
  if(Soil_Variable[i] == "M"){New_classification <- data.frame(New_ID = c(1:length(levels(RF_Train$pred$obs))), Attribute = M_legend_order)}
  if(Soil_Variable[i] == "T"){New_classification <- data.frame(New_ID = c(1:length(levels(RF_Train$pred$obs))), Attribute = T_legend_order)}                                   
  
  reclass_df <- merge(Current_classification, New_classification)
  reclass_df$from <- reclass_df$Current_ID - 0.5; reclass_df$to <- reclass_df$Current_ID + 0.5
  reclass_df <- reclass_df[,c("from", "to", "New_ID"),]
  
  Map_reclass <- reclassify(Map_reclass, reclass_df)
  
  Map_reclass_attributes <- list(data.frame(ID = 1:length(levels(RF_Train$pred$obs)), category = New_classification$Attribute))
  Map_reclass@data@attributes <- Map_reclass_attributes; Map_reclass@data@isfactor <- TRUE

  writeRaster(Map_reclass, filename = paste0(wd, "/layers/Raw/", Plain_Response_names[i], "_RF.tif"), format="GTiff", overwrite=TRUE)
  
  if(Create_prob_raster == TRUE){
  prob_start <- Sys.time()
  writeRaster(clusterR(Covariates, fun = raster::predict, args = list(RF_Train, type = "prob", index = 1:nlevels(RF_Train))),
              filename= paste0(wd, "/layers/Probability/", Plain_Response_names[i], "_prob.tif"),
              format="GTiff", overwrite=TRUE, bylayer=TRUE, suffix = levels(RF_Train))
  #writeRaster(x=raster::predict(Covariates, RF_Train, type = "prob", index = 1:nlevels(RF_Train)), filename= paste0(wd, "/layers/Probability/", Plain_Response_names[i], "_prob.tif"),
   #           format="GTiff", overwrite=TRUE, bylayer=TRUE, suffix = levels(RF_Train))
  prob_end <- Sys.time(); print(capture.output(prob_end - prob_start))
  }

}

endCluster()
end_time <- Sys.time()
print(paste("End of Uncertainty Models.R ", capture.output(end_time - start_time)))
#####
####
###
##
#END
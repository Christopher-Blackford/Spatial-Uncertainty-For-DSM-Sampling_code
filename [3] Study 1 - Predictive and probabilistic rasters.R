###############################################################################################################
###############################################################################################################
#Uncertainty Models predictive raster
#Code by: Christopher Blackford (christopher.blackford@canada.ca)
start_time <- Sys.time()

#Loading libraries
library(raster); library(rgdal); library(rgeos)
library(caret)
library(doParallel)

#####Specifying resolution and filter
Raw_DEM_reso = 10
Filter_type = 40

M_legend_order <- c("Xeric", "Dry", "Moderately_Dry", "Moderately_Fresh", "Fresh", 
                    "Very_Fresh", "Moderately_Moist", "Moist", "Very_Moist", 
                    "Moderately_Wet", "Wet", "Very_Wet", "Saturated")

T_legend_order <- c("Organic", "Silty_Clay", "Silt", "Silt_Loam",
                    "Silty_Clay_Loam", "Loam", "Clay_Loam", 
                    "Sandy_Clay_Loam", "Silty_Sand", 
                    "Sandy_Loam","Loamy_Sand", "Sand")

#Get accuracy from a list of model object
get_Accuracy <- function(x){max(x$results$Accuracy)}
#wd
wd <- paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/Study1/")

#Create directories
if(!dir.exists(paste0(wd, "predicted_layers"))){dir.create(paste0(wd, "predicted_layers"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round1"))){dir.create(paste0(wd, "predicted_layers/Round1"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round2"))){dir.create(paste0(wd, "predicted_layers/Round2"))}

if(!dir.exists(paste0(wd, "predicted_layers/Round1/Raw"))){dir.create(paste0(wd, "predicted_layers/Round1/Raw"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round2/Raw"))){dir.create(paste0(wd, "predicted_layers/Round2/Raw"))}

if(!dir.exists(paste0(wd, "predicted_layers/Round1/Probabilities"))){dir.create(paste0(wd, "predicted_layers/Round1/Probabilities"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round2/Probabilities"))){dir.create(paste0(wd, "predicted_layers/Round2/Probabilities"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round1/Probabilities/layers"))){dir.create(paste0(wd, "predicted_layers/Round1/Probabilities/layers"))}
if(!dir.exists(paste0(wd, "predicted_layers/Round2/Probabilities/layers"))){dir.create(paste0(wd, "predicted_layers/Round2/Probabilities/layers"))}

if(!dir.exists(paste0("./temp"))){dir.create(paste0("./temp"))} #These lines are useful as they store potentially large raster temporary files in a recognizable location instead of default C drive
rasterOptions(tmpdir="./temp")

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

#Running in parallel to predict raster on covariates
detectCores() #Should be 12 on SOIL-BOT
beginCluster(detectCores()-1)

##########
##########
###Round 1 maps
Soil_Variable = c("M", "T")

for(i in 1:length(Soil_Variable)){
R1_Models <- list.files(paste0(wd, "Model_Comparison/models/Round1"), pattern = paste0(Soil_Variable[i], ".rds"))
R1_Models_readme <- list.files(paste0(wd, "Model_Comparison/models/Round1"), pattern = paste0(Soil_Variable[i], ".rds"), full.names = T)

Models_list = lapply(FUN = readRDS, R1_Models_readme)
Model_accuracy = unlist(lapply(FUN = get_Accuracy, Models_list))

#if divisible by 2, choose higher model to get median value
if(length(Model_accuracy) %% 2 == 0){Model_select = length(Model_accuracy)/2 + 1}else{Model_select = length(Model_accuracy)/2}
Model_map <- Models_list[[order(Model_accuracy)[Model_select]]]

C_start <- Sys.time()
Map <- clusterR(Covariates, fun = raster::predict, args = list(Model_map, type = "raw"))
C_end <- Sys.time(); capture.output(C_end - C_start)


#Temporary reclass to organize for ordinal reclass
reclass_matrix <- matrix(c(0.5:(length(levels(Model_map$pred$obs))-0.5), 1.5:(length(levels(Model_map$pred$obs))+0.5), 101:(length(levels(Model_map$pred$obs))+100)),
                    nrow = length(levels(Model_map$pred$obs)), ncol = 3)
Map_reclass <- reclassify(Map, reclass_matrix)

Current_classification <- data.frame(Current_ID = reclass_matrix[,3], Attribute = levels(Model_map$pred$obs))
if(Soil_Variable[i] == "M"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = M_legend_order)}
if(Soil_Variable[i] == "T"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = T_legend_order)}                                   
                                   
reclass_df <- merge(Current_classification, New_classification)
reclass_df$from <- reclass_df$Current_ID - 0.5; reclass_df$to <- reclass_df$Current_ID + 0.5
reclass_df <- reclass_df[,c("from", "to", "New_ID"),]

Map_reclass <- reclassify(Map_reclass, reclass_df)

Map_reclass_attributes <- list(data.frame(ID = 1:length(levels(Model_map$pred$obs)), category = New_classification$Attribute))
Map_reclass@data@attributes <- Map_reclass_attributes; Map_reclass@data@isfactor <- TRUE

writeRaster(x=Map_reclass, 
            filename = paste0(wd, "predicted_layers/Round1/Raw/", gsub(R1_Models[order(Model_accuracy)[Model_select]], pattern = ".rds", replacement = ""), ".tif"),
            format="GTiff", overwrite=TRUE)
}


##########
##########
###Round 2 maps
Uncertainty_metric_names <- c("IU", "EU", "CI")

for(i in 1:length(Soil_Variable)){
#Random
R2_Random <- list.files(paste0(wd, "Model_Comparison/models/Round2"), pattern = paste0("Random_", Soil_Variable[i], ".rds"))
R2_Random_readme <- list.files(paste0(wd, "Model_Comparison/models/Round2"), pattern = paste0("Random_", Soil_Variable[i], ".rds"), full.names = T)

Models_list = lapply(FUN = readRDS, R2_Random_readme)
Model_accuracy = unlist(lapply(FUN = get_Accuracy, Models_list))

#if divisible by 2, choose higher model to get median value
if(length(Model_accuracy) %% 2 == 0){Model_select = length(Model_accuracy)/2 + 1}else{Model_select = length(Model_accuracy)/2}
Model_map <- Models_list[[order(Model_accuracy)[Model_select]]]

Map <- clusterR(Covariates, fun = raster::predict, args = list(Model_map, type = "raw"))

#Temporary reclass to organize for ordinal reclass
reclass_matrix <- matrix(c(0.5:(length(levels(Model_map$pred$obs))-0.5), 1.5:(length(levels(Model_map$pred$obs))+0.5), 101:(length(levels(Model_map$pred$obs))+100)),
                         nrow = length(levels(Model_map$pred$obs)), ncol = 3)
Map_reclass <- reclassify(Map, reclass_matrix)

Current_classification <- data.frame(Current_ID = reclass_matrix[,3], Attribute = levels(Model_map$pred$obs))
if(Soil_Variable[i] == "M"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = M_legend_order)}
if(Soil_Variable[i] == "T"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = T_legend_order)}                                   

reclass_df <- merge(Current_classification, New_classification)
reclass_df$from <- reclass_df$Current_ID - 0.5; reclass_df$to <- reclass_df$Current_ID + 0.5
reclass_df <- reclass_df[,c("from", "to", "New_ID"),]

Map_reclass <- reclassify(Map_reclass, reclass_df)

Map_reclass_attributes <- list(data.frame(ID = 1:length(levels(Model_map$pred$obs)), category = New_classification$Attribute))
Map_reclass@data@attributes <- Map_reclass_attributes; Map_reclass@data@isfactor <- TRUE

writeRaster(x=Map_reclass, 
            filename = paste0(wd, "predicted_layers/Round2/Raw/", gsub(R2_Random[order(Model_accuracy)[Model_select]], pattern = ".rds", replacement = ""), ".tif"),
            format="GTiff", overwrite=TRUE)

#Treatment
for(j in 1:length(Uncertainty_metric_names)){
R2_Treatment <- list.files(paste0(wd, "Model_Comparison/models/Round2"), pattern = paste0(Uncertainty_metric_names[j], "_High_", Soil_Variable[i], ".rds"))
R2_Treatment_readme <- list.files(paste0(wd, "Model_Comparison/models/Round2"), pattern = paste0(Uncertainty_metric_names[j], "_High_", Soil_Variable[i], ".rds"), full.names = T)

Models_list = lapply(FUN = readRDS, R2_Treatment_readme)
Model_accuracy = unlist(lapply(FUN = get_Accuracy, Models_list))

#if divisible by 2, choose higher model to get median value
if(length(Model_accuracy) %% 2 == 0){Model_select = length(Model_accuracy)/2 + 1}else{Model_select = length(Model_accuracy)/2}
Model_map <- Models_list[[order(Model_accuracy)[Model_select]]]

Map <- clusterR(Covariates, fun = raster::predict, args = list(Model_map, type = "raw"))

#Temporary reclass to organize for ordinal reclass
reclass_matrix <- matrix(c(0.5:(length(levels(Model_map$pred$obs))-0.5), 1.5:(length(levels(Model_map$pred$obs))+0.5), 101:(length(levels(Model_map$pred$obs))+100)),
                         nrow = length(levels(Model_map$pred$obs)), ncol = 3)
Map_reclass <- reclassify(Map, reclass_matrix)

Current_classification <- data.frame(Current_ID = reclass_matrix[,3], Attribute = levels(Model_map$pred$obs))
if(Soil_Variable[i] == "M"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = M_legend_order)}
if(Soil_Variable[i] == "T"){New_classification <- data.frame(New_ID = c(1:length(levels(Model_map$pred$obs))), Attribute = T_legend_order)}                                   

reclass_df <- merge(Current_classification, New_classification)
reclass_df$from <- reclass_df$Current_ID - 0.5; reclass_df$to <- reclass_df$Current_ID + 0.5
reclass_df <- reclass_df[,c("from", "to", "New_ID"),]

Map_reclass <- reclassify(Map_reclass, reclass_df)

Map_reclass_attributes <- list(data.frame(ID = 1:length(levels(Model_map$pred$obs)), category = New_classification$Attribute))
Map_reclass@data@attributes <- Map_reclass_attributes; Map_reclass@data@isfactor <- TRUE

writeRaster(x=Map_reclass, 
            filename = paste0(wd, "predicted_layers/Round2/Raw/", gsub(R2_Treatment[order(Model_accuracy)[Model_select]], pattern = ".rds", replacement = ""), ".tif"),
            format="GTiff", overwrite=TRUE)
}
}

#####
#####
#####
#####

endCluster()
end_time <- Sys.time()
print(paste("End of [3] Study 1 - Predictive and probabilistic rasters.R ", capture.output(end_time - start_time)))
#####
####
###
##
#END
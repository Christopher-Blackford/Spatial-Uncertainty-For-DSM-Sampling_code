##############################################################################################################################
###Load soil pedon dataset and summary histogram .R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

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

#####
#####Load in soil dataset if it exists, create it if it doesn't
if(file.exists(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data/", Raw_DEM_reso, "DEM", Filter_type,"Response_df.csv"))){print("Loading previous filter dataset"); df <- read.csv(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data/", Raw_DEM_reso, "DEM", Filter_type,"Response_df.csv"))
}else{print("Creating filter dataset")
  #Loading response variables (moisture, texture, depth class, depth mottle, mottle class)
  Response_df <- read.csv("./Models/Soil_data/Soil_response_data.csv")
  Response_df$ID <- 1:nrow(Response_df)
  
  xy <- subset(Response_df, select = c(X, Y))
  Response_points <- SpatialPointsDataFrame(coords = xy, data = Response_df, proj4string = crs(Covariates)) #Covariates@crs may need to be replaced with CRS(Covariates)
  
  #Extract covariate points to response data
  rasValue <- raster::extract(Covariates, Response_points, method = "simple")
  df <- as.data.frame(rasValue); df$ID <- 1:nrow(df)
  df <- merge(Response_df, df) #Quite a couple of non-overlapping points actually, double check this?
  df <- df[complete.cases(df),]
  
  #Remove non Moisture Regime and Textural class attributes
  df <- within(df, rm("DepthClass", "DepthMottl", "MottlClass"))
  
  row.names(df) <- 1:nrow(df)
  write.csv(df, paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data/", Raw_DEM_reso, "DEM", Filter_type,"Response_df.csv"), row.names=FALSE)
  df <- read.csv(paste0("./Uncertainty_Models/", Raw_DEM_reso, "DEM", Filter_type, "/data/", Raw_DEM_reso, "DEM", Filter_type,"Response_df.csv")) #Need to read csv so it treats geo and bio colums as integers
  
  rm(xy,Response_df,rasValue,Response_points)
}

df$X <- as.numeric(df$X); df$Y <- as.numeric(df$Y)
Predictor_variable_names <- paste(names(Covariates), collapse="+")
#####
#####

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

###Plot histograms of soil dataset values
#Create directory for output
if(!dir.exists(paste0(wd, "/plots"))){dir.create(paste0(wd, "/plots"))}

#####Histogram of soil attribute classes
#MR
MR_df <- as.data.frame(table(df[df$MR != "",][,"MR"])*100/sum(table(df[df$MR != "",][,"MR"])))
MR_df <- MR_df[MR_df$Var1 != "",]

ggplot(data=MR_df, aes(x=reorder(Var1, -Freq), y=Freq))+
  geom_bar(stat="identity")+
  labs(title = "Moisture Regime Class Histogram", x = "Class", y = "Percent of Soil Dataset")+
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.25))

ggsave(paste0(wd, "/plots/M_Hist.png"), width = 10, height = 6)

#T
T_df <- as.data.frame(table(df[df$TEXT != "",][,"TEXT"])*100/sum(table(df[df$TEXT != "",][,"TEXT"])))
T_df <- T_df[T_df$Var1 != "",]

ggplot(data=T_df, aes(x=reorder(Var1, -Freq), y=Freq))+
  geom_bar(stat="identity")+
  labs(title = "Textural Class Histogram", x = "Class", y = "Percent of Soil Dataset")+
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.25))

ggsave(paste0(wd, "/plots/T_Hist.png"), width = 10, height = 6)

#####Clean up R global environment
rm(Bio_categoricalCovariates, Bio_continuousCovariates, Geo_Covariates, Hydro_Covariates, Dist_Covariates)
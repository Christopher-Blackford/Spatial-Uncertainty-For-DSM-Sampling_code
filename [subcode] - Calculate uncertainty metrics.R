##############################################################################################################################
####Calculate uncertainty metrics.R
#Code by: Christopher Blackford (christopher.blackford@canada.ca)

###Ignorance uncertainty - Ranges from 0(high certainty) to 1 (low certainty)
f1 <- function(x){x*log(x+0.000001)}
Ignorance_Uncertainty <- as.data.frame(apply(D2_uncertainty, MARGIN = c(1,2), FUN = f1))
Ignorance_Uncertainty$IU <- apply(Ignorance_Uncertainty, MARGIN = 1, FUN = sum)
Ignorance_Uncertainty$IU <- -1*Ignorance_Uncertainty$IU/log(nlevels(Round1_Model))
Ignorance_Uncertainty$IU[Ignorance_Uncertainty$IU < 0] <- 0
Ignorance_Uncertainty <- Ignorance_Uncertainty[("IU")]

###Exaggeration uncertainty
Exaggeration_Uncertainty <- data.frame(EU = apply(D2_uncertainty, MARGIN = 1, FUN = max))
Exaggeration_Uncertainty$EU <- 1 - Exaggeration_Uncertainty$EU

###Confusion index
f2 <- function(x){1 - (max(x) - max(x[x != max(x)]))}
Confusion_Index <- data.frame(CI = apply(D2_uncertainty, MARGIN = 1, FUN = f2))

#Combine metrics
Dataset2 <- merge(Dataset2, cbind(Ignorance_Uncertainty, Exaggeration_Uncertainty, Confusion_Index), by = "row.names"); Dataset2$Row.names = NULL
#####
rm(Ignorance_Uncertainty, Exaggeration_Uncertainty, Confusion_Index)
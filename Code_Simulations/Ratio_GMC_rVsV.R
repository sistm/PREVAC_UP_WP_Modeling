# ----------------------- #
# OBJECTIVE: Calculation of the ratio of geometric mean concentration (GMC) of anti-EBOV GP antibody, either at the peak or at 1 year
#            Focus on the vaccine regimen rVsV
# Date (last update): 03/09/2023
# Author: Marie Alexandre
# ----------------------- #



rm(list=ls())

# ----- LIBRARIES ----- #
library(lme4) ; library(plyr); library(psych)
library(parallel) ; library(snow) ; library(doSNOW) ; library(foreach)
# --------------------- #


# ----- FUNCTIONS ----- #
Bootstrap_dataset_creation <- function(original_data,nb_ind){
  list_ids <- unique(original_data$IDs)
  # sampling with replacement at individual level
  sampled_ind <- sample(list_ids,replace=TRUE)
  # creation of the resulting dataset
  boot_dataset <- do.call("rbind",lapply(seq(1,length(sampled_ind)),function(i,data){
    # data <- original_data ; i <- 1
    tmp_data <- subset(data,IDs == sampled_ind[i])
    # rename of SUBJID_ID
    tmp_data$IDs_BOOT <- i
    tmp_data <- tmp_data[,c("IDs_BOOT",colnames(data))]
    return(tmp_data)
  },data=original_data))
  return(boot_dataset)
}
Ratios_GMT_Estimation <- function(data,targeted_variable,reference){
  library(plyr)
  
  # Estimation of geometric mean for each category of the targeted variable
  Summary_GMT <- ddply(.data=data,.variables = c(targeted_variable), summarise, GMT=geometric.mean(PRED))
  colnames(Summary_GMT)[1] <- "Category"
  
  # Estimation of ratios of GMT category of interest Vs reference
  reference_value <- Summary_GMT$GMT[which(Summary_GMT$Category == reference)]
  Summary_GMT$Ratio <- Summary_GMT$GMT/reference_value
  Summary_GMT <- cbind(Variable=targeted_variable,Summary_GMT,stringsAsFactors = FALSE)
  return(Summary_GMT)
}
# --------------------- #


'%notin%' <- Negate('%in%')

# -- Load of data ####
load("Data/dfprevacrVSV.rda")
dfprevac <- dfprevac[,c("SUBJID_N", "VISIT","COUNTRYID", "ARM", "SEX", "CATAGE", "ABBn","TEST_MED_FINAL_LOG")]
colnames(dfprevac) <- c("IDs","TIME","COUNTRYID", "ARM", "SEX", "CATAGE", "ABBn","LOG10AB")

dfprevac$RESCALED_TIME <- dfprevac$TIME-28  # rescale of time

rVSV_participants <- unique(dfprevac$IDs)
rVSV_Nbpart <- length(rVSV_participants)
# ----- #


# -- Calculation of the ratio by Bootstrap
Nb_bootstrap <- 2000
nb_cores <- 5

# Variables of parallel calculation
NCores <- min(detectCores()-1,nb_cores)
Cluster <- parallel::makeCluster(NCores,type="SOCK")
registerDoSNOW(Cluster)
pb <- txtProgressBar(max=Nb_bootstrap,style=3)         # We add a progress bar
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)
clusterExport(Cluster,c("Bootstrap_dataset_creation","dfprevac","rVSV_Nbpart","Ratios_GMT_Estimation"))

Bootstraps_ratios <- foreach(n=seq(1,Nb_bootstrap),.combine="rbind",.errorhandling = "stop",.packages = c("parallel","foreach","doSNOW","plyr","lme4","psych"),.options.snow = opts)%dopar%{
  
  # -- Creation of the bootstrap dataset
  boot_dataset <- Bootstrap_dataset_creation(original_data=dfprevac,nb_ind=rVSV_Nbpart)
  
  # -- Estimation of the linear mixed-effect model
  boot_MEM_rVSV <- lmer(data = boot_dataset,
                        formula = LOG10AB ~ RESCALED_TIME*ARM*SEX + COUNTRYID + ABBn + CATAGE + COUNTRYID*ARM + CATAGE*ARM + ABBn*ARM + (1 + RESCALED_TIME | IDs_BOOT), REML = T)
  
  # -- Addition of model predictions in the bootstrap dataset
  boot_dataset$LOG10PRED <- predict(boot_MEM_rVSV)
  boot_dataset$PRED <- 10^(boot_dataset$LOG10PRED)
  boot_dataset$ABBn <- ifelse(boot_dataset$ABBn >= log10(200), 1, 0)
  # We only consider vaccinated participants
  boot_dataset <- subset(boot_dataset,ARM == 1)
  
  # -- Estimation of ratios of GMT
  # Ratios Age category
  Ratios_CATAGE <- rbind(cbind(Boot=n,Time=0,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 0),targeted_variable="CATAGE",reference="Adults")),
                         cbind(Boot=n,Time=337,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 365-28),targeted_variable="CATAGE",reference="Adults")))
  
  # Ratios Sex 
  Ratios_SEX <- rbind(cbind(Boot=n,Time=0,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 0),targeted_variable="SEX",reference="Man")),
                      cbind(Boot=n,Time=337,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 365-28),targeted_variable="SEX",reference="Man")))
  
  # Ratios Country 
  Ratios_COUNTRY <- rbind(cbind(Boot=n,Time=0,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 0),targeted_variable="COUNTRYID",reference="Mali")),
                          cbind(Boot=n,Time=337,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 365-28),targeted_variable="COUNTRYID",reference="Mali")))
  
  # Ratios Pre-vaccine antibody concentration
  Ratios_ABBn <- rbind(cbind(Boot=n,Time=0,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 0),targeted_variable="ABBn",reference=0)),
                       cbind(Boot=n,Time=337,Ratios_GMT_Estimation(data=subset(boot_dataset,RESCALED_TIME == 365-28),targeted_variable="ABBn",reference=0)))
  Ratios_ABBn$Category <- as.character(Ratios_ABBn$Category) 
  
  Results_Ratio <- rbind(Ratios_CATAGE,Ratios_SEX,Ratios_COUNTRY,Ratios_ABBn)
  Results_Ratio$Category <- as.character(Results_Ratio$Category)
  
  return(Results_Ratio)
}
close(pb)
stopCluster(Cluster)
# ----- #


# -- Estimation of the distribution of the ratio for each model variable
GMC_Ratio_Distribution <- ddply(.data=Bootstraps_ratios,.variables = .(Time,Variable,Category),summarise,
                                Mean=mean(Ratio,na.rm=TRUE),ICMIN=quantile(Ratio,probs=c(0.025),na.rm=TRUE),
                                ICMAX=quantile(Ratio,probs=c(0.975),na.rm=TRUE))
# ----- #

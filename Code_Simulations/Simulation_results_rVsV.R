# ----------------------- #
# OBJECTIVE: Generate results to estimate both the time requested to reach a given Ab threshold and the probability to reach it at day 365
#            Focus on the vaccine regimen rVsV
# Date (last update): 03/09/2023
# Author: Marie Alexandre

# Description: This file is built in 2 parts:
#                1) Estimation of the mixed-effects model 
#                2) Generation of results
#                     a) Simulation of data
#                     b) Calculation of probabilities to reach a given Ab threshold
#                     c) Calculation of the time requested to reach a given Ab threshold 

# Additional information: This code has been written to be launch on a calculation server (use of batch file)
# ----------------------- #





rm(list=ls())




# ----- LIBRARIES ----- #
library(slurmR)   # for calculation of remote server

library(lme4) # for the function lmer
library(plyr) # for the function ddply
library(mvtnorm)  # for the function rmvnorm
# libraries for parallel calculation
library(parallel) ; library(snow) ; library(doSNOW) ; library(foreach)
# --------------------- #




# ----- FUNCTIONS ----- #
Simulation_database_Function <- function(obs_data,profiles,nb_ind){
  
  # -- Simulation of the categorical covariates
  profiles_selection <- sample(seq(1,nrow(profiles)),size=nb_ind,replace=TRUE,prob = profiles$Proportion)
  simulated_individuals <- profiles[profiles_selection,]
  
  # -- Addition of continuous covariates (Abbn)
  # simulated_individuals$ABBn <- sapply(seq(1,nrow(simulated_individuals)),function(i,data) rnorm(1,mean=data$Mean_Abbn,sd=data$sd_Abbn),data=simulated_individuals)
  simulated_individuals$ABBn <- sapply(seq(1,nrow(simulated_individuals)),function(i,data){
    tmp_data <- data[which(data$COUNTRY == simulated_individuals$COUNTRY[i] & data$ARM == simulated_individuals$ARM[i] & data$SEX == simulated_individuals$SEX[i] &
                             data$CATAGE == simulated_individuals$CATAGE[i] & data$ABBnInf200 == simulated_individuals$ABBnInf200[i]),]
    tmp_data <- unique(tmp_data[,c("ID","COUNTRY","ARM","SEX","CATAGE","ABBn","ABBnInf200")])
    
    if(nrow(tmp_data) == 1){
      tmp <-tmp_data$ABBn[1]
    }else{
      tmp <- sample(x=tmp_data$ABBn,size =1,replace = T)
    }
    return(tmp)
  },data=obs_data)
  
  simulated_individuals <- cbind(ID=seq(1,nb_ind),simulated_individuals[,c("COUNTRY","ARM","SEX","CATAGE","ABBnInf200","ABBn")])
  rownames(simulated_individuals) <- seq(1,nrow(simulated_individuals))
  
  return(simulated_individuals)
}
Simulation_Population_Parameters_Function <- function(parameters,varcov_mat,nb=1){
  simulated_parameters <- rmvnorm(nb,mean=parameters,sigma = varcov_mat)
  return(simulated_parameters)
}
Simulation_Random_Effects_Function <- function(varcov_Sd_RE,sd_intercept,sd_slope,correlation,nb_ind){
  
  # Simulation of standard deviation
  simulated_Sd_RE <- as.numeric(rmvnorm(1,mean=c(sd_intercept,sd_slope),sigma=varcov_Sd_RE))
  # Simulation of random effects for each individual
  varcov_RE <- rbind(c(simulated_Sd_RE[1]^2,simulated_Sd_RE[1]*simulated_Sd_RE[2]*correlation),
                     c(simulated_Sd_RE[1]*simulated_Sd_RE[2]*correlation,simulated_Sd_RE[2]^2))
  simulated_RE <- rmvnorm(nb_ind,mean=c(0,0),sigma=varcov_RE)
  colnames(simulated_RE) <- c("Gam0","Gam1")
  
  return(simulated_RE)
}
Model_prediction_Function <- function(fixed_params,random_effects,time,sd_Err,covariate){
  
  with(as.list(c(fixed_params)),{
    
    if(nrow(covariate)>1){
      dfsim <- cbind(covariate,random_effects)
    }else{
      dfsim <- covariate
      dfsim$Gam0 <- random_effects["Gam0"]
      dfsim$Gam1 <- random_effects["Gam1"]
    }
    
    dfsim$ARM <- as.numeric(dfsim$ARM)
    
    dfsim$SEX <- ifelse(dfsim$SEX == "Man", 0, 1)
    
    dfsim$CATAGE1 <- ifelse(dfsim$CATAGE == "1-4", 1, 0)
    dfsim$CATAGE5 <- ifelse(dfsim$CATAGE == "5-11", 1, 0)
    dfsim$CATAGE12 <- ifelse(dfsim$CATAGE == "12-17", 1, 0)
    
    dfsim$SIERRA <- ifelse(dfsim$COUNTRY == "Sierra Leone", 1, 0)
    dfsim$LIBERIA <- ifelse(dfsim$COUNTRY == "Liberia", 1, 0)
    dfsim$GUINEA <- ifelse(dfsim$COUNTRY == "Guinea", 1, 0)
    
    Model_prediction <- NULL
    nb_ind <- nrow(dfsim)
    for(i in 1:nb_ind){
      # i <- 1
      model_pred_ind <- data.frame(dfsim[i,],Time=Time_simulation-28)
      
      model_pred_ind$log10Ab <- B0 + model_pred_ind$Gam0 + (B1+model_pred_ind$Gam1)*model_pred_ind$Time +
        B2*model_pred_ind$ARM + B3*model_pred_ind$ARM*model_pred_ind$Time + B4*model_pred_ind$CATAGE1 + B5*model_pred_ind$CATAGE5 + B6*model_pred_ind$CATAGE12 +
        B7*model_pred_ind$CATAGE1*model_pred_ind$ARM + B8*model_pred_ind$CATAGE5*model_pred_ind$ARM + B9*model_pred_ind$CATAGE12*model_pred_ind$ARM +
        B10*model_pred_ind$SEX + B11*model_pred_ind$SEX*model_pred_ind$ARM + B12*model_pred_ind$SEX*model_pred_ind$Time +
        B13*model_pred_ind$GUINEA + B14*model_pred_ind$SIERRA + B15*model_pred_ind$LIBERIA +
        B16*model_pred_ind$GUINEA*model_pred_ind$ARM + B17*model_pred_ind$SIERRA*model_pred_ind$ARM + B18*model_pred_ind$LIBERIA*model_pred_ind$ARM +
        B19*model_pred_ind$ABBn + B20*model_pred_ind$ABBn*model_pred_ind$ARM +
        B21*model_pred_ind$ARM*model_pred_ind$SEX*model_pred_ind$Time
      # Addition of the noise
      model_pred_ind$log10Ab <- model_pred_ind$log10Ab + rnorm(nrow(model_pred_ind),mean=0,sd=sd_Err)
      
      Model_prediction <- rbind(Model_prediction,model_pred_ind)
    }
    return(Model_prediction)
  })
}
Proba_Above_AbThres_Function <- function(data,time,thres){
  
  return(length(unique(data$ID[which(data$Time == time & data$Ab >= thres)]))/length(unique(data$ID)))
}
Time_AbTrhreshold_Estimation_Function <- function(fixed_params,random_effects,sd_Err,covariate,Ab_thresh){
  
  with(as.list(c(fixed_params)),{
    dfsim <- cbind(covariate,random_effects)
    dfsim$ARM <- as.numeric(dfsim$ARM)
    dfsim$SEX <- ifelse(dfsim$SEX == "Man", 0, 1)
    dfsim$CATAGE1 <- ifelse(dfsim$CATAGE == "1-4", 1, 0)
    dfsim$CATAGE5 <- ifelse(dfsim$CATAGE == "5-11", 1, 0)
    dfsim$CATAGE12 <- ifelse(dfsim$CATAGE == "12-17", 1, 0)
    dfsim$SIERRA <- ifelse(dfsim$COUNTRY == "Sierra Leone", 1, 0)
    dfsim$LIBERIA <- ifelse(dfsim$COUNTRY == "Liberia", 1, 0)
    dfsim$GUINEA <- ifelse(dfsim$COUNTRY == "Guinea", 1, 0)
    
    Time_prediction <- NULL
    nb_ind <- nrow(dfsim)
    for(j in 1:nb_ind){
      time_pred_ind <- data.frame(dfsim[j,],Ab_thresh=Ab_thresh)
      ind_slope <- as.numeric((B1 + B3*time_pred_ind$ARM + B12*time_pred_ind$SEX + B21*time_pred_ind$ARM*time_pred_ind$SEX + time_pred_ind$Gam1))
      sign_ind_slope <- sign(ind_slope)
      
      if(sign_ind_slope == 1){
        time_pred_ind$time <- NA
        time_pred_ind$objective_time <- NA
      }else{
        estimated_time <- optimize(function(time) abs((B0 + time_pred_ind$Gam0 + (B1+time_pred_ind$Gam1)*time +
                                                         B2*time_pred_ind$ARM + B3*time_pred_ind$ARM*time + B4*time_pred_ind$CATAGE1 + B5*time_pred_ind$CATAGE5 + B6*time_pred_ind$CATAGE12 +
                                                         B7*time_pred_ind$CATAGE1*time_pred_ind$ARM + B8*time_pred_ind$CATAGE5*time_pred_ind$ARM + B9*time_pred_ind$CATAGE12*time_pred_ind$ARM +
                                                         B10*time_pred_ind$SEX + B11*time_pred_ind$SEX*time_pred_ind$ARM + B12*time_pred_ind$SEX*time +
                                                         B13*time_pred_ind$GUINEA + B14*time_pred_ind$SIERRA + B15*time_pred_ind$LIBERIA +
                                                         B16*time_pred_ind$GUINEA*time_pred_ind$ARM + B17*time_pred_ind$SIERRA*time_pred_ind$ARM + B18*time_pred_ind$LIBERIA*time_pred_ind$ARM +
                                                         B19*time_pred_ind$ABBn + B20*time_pred_ind$ABBn*time_pred_ind$ARM +
                                                         B21*time_pred_ind$ARM*time_pred_ind$SEX*time + rnorm(1,mean=0,sd=sd_Err)) - Ab_thresh),interval = c(0,1e10),tol=1e-6)
        
        time_pred_ind$time <- estimated_time$minimum
        time_pred_ind$objective_time <- estimated_time$objective
      }
      Time_prediction <- rbind(Time_prediction,time_pred_ind)
    }
    return(Time_prediction)
  })
}
# --------------------- #



# ----- Arguments for calculation server ----- #
# See the batch file "Simulation_results_rVsV.sh" for an example of run
args <- commandArgs(trailingOnly = TRUE)             # List of arguments called in the batch file 
slurm_arrayid <- Slurm_env(x='SLURM_ARRAY_TASK_ID')  # calculation was splitted into multiple jobs
arr <- as.numeric(slurm_arrayid)
# --------------------- #







# ****************************************************** ####
  #       PART 1: Estimation of the MEM on raw data        ####
# ***************************************************** ####
# -- Load of data ####
load("Data/dfprevacrVSV.rda")

dfprevac$New_VISIT <- dfprevac$VISIT-28  # rescale of time
# ----- #

# -- Estimation of the Mixed effects model ####
MEM_rVsV <- lmer(TEST_MED_FINAL_LOG ~ New_VISIT*ARM*SEX + COUNTRYID + ABBn + CATAGE + COUNTRYID*ARM + CATAGE*ARM + ABBn*ARM + (1 + New_VISIT | SUBJID_N) , dfprevac, REML = T)

# Fixed parameters
Fixed_parameters_order <- c("B0", "B1", "B2", "B10", "B13", "B15", "B14", "B19", "B4", "B6", "B5", "B3",
                            "B12", "B11", "B16", "B18", "B17", "B7", "B9", "B8", "B20", "B21")
MEM_Fixed_parameters <-  summary(MEM_rVsV)$coefficients[,1]
names(MEM_Fixed_parameters) <- Fixed_parameters_order
MEM_VCoV_Mat <- as.matrix(vcov(MEM_rVsV))
colnames(MEM_VCoV_Mat) <- Fixed_parameters_order ; rownames(MEM_VCoV_Mat) <- Fixed_parameters_order
# Random effects
Sd_Intercept <- attr(VarCorr(MEM_rVsV)$SUBJID_N, "stddev")[1]
Sd_Visit <- attr(VarCorr(MEM_rVsV)$SUBJID_N, "stddev")[2]
Corr_RE <- attr(VarCorr(MEM_rVsV)$SUBJID_N, "correlation")[1,2]
RE_VCoC_Mat <- rbind(c(1.646655e-05,3.691077e-11),c(3.691077e-11,1.758331e-15))  # Values extracted from the estimated covariance matrix
# # Error model
Sd_Residual <- attr(VarCorr(MEM_rVsV),"sc")
Res_Var <- 8.308073e-07 # Values extracted from the estimated covariance matrix
# ----- #

# -- Extraction of the database required for the simulations
Data <- dfprevac[,c("SUBJID_N","VISIT","New_VISIT","COUNTRYID","ARM","SEX","TEST_MED_FINAL_LOG","CATAGE","ABBn")]
colnames(Data) <- c("ID","TIME","New_TIME","COUNTRY","ARM","SEX","LOG10Ab","CATAGE","ABBn")
Data$ABBnInf200 <- 1*(Data$ABBn < log10(200))
# ----- #

# -- Estimation of the proportion of each type of profile in the population
Profiles <- ddply(.data=Data,.variables = .(COUNTRY,ARM,SEX,CATAGE,ABBnInf200),summarise,
                  N=length(unique(ID)),Proportion = length(unique(ID))/length(unique(Data$ID)),
                  Mean_Abbn = mean(ABBn,na.rm=T),sd_Abbn = sd(ABBn,na.rm=T))
# ----- #








# ****************************************** ####
#       PART 2: Generation of results        ####
# **************************************** ####
Folder_results <- "Rcode/Results/rVSV_2"
dir.create(Folder_results,recursive = TRUE)

Nb_Individual <- as.numeric(args[1]) 
Ab_thresholds <- seq(100,5000,by=100)


simus <- ((arr-1)*10+1):(arr*10)

for(n in simus){

  # -- Data simulation ####
  
  print("-- Simulation of individuals")
  # - Simulation of individuals
  Simulated_Individual <- Simulation_database_Function(obs_data = Data,profiles=Profiles,nb_ind=Nb_Individual)
  
  # - Simulation of model parameters
  # Time of simulation
  Time_simulation <- sort(unique(Data$TIME))
  
  print("-- Simulation of population parameters")
  # Fixed parameters
  Simulated_parameters <- Simulation_Population_Parameters_Function(parameters=MEM_Fixed_parameters,varcov_mat=MEM_VCoV_Mat)[1,]
  
  print("-- Simulation of random effects")
  # Random Effects
  Simulated_RE <- Simulation_Random_Effects_Function(varcov_Sd_RE=RE_VCoC_Mat,sd_intercept=Sd_Intercept,sd_slope=Sd_Visit,correlation=Corr_RE,nb_ind=Nb_Individual)
  
  print("-- Simulation of residual error")
  # Sd Error model
  Simulated_sd_ErrMod <- rnorm(1,Sd_Residual,Res_Var/(4*(sqrt(Sd_Residual))))
  
  # - Model prediction
  print("-- Simulation of data")
  NCores <- detectCores()-1
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  clusterExport(Cluster,c("Model_prediction_Function","Simulated_parameters","Simulated_RE","Time_simulation","Simulated_sd_ErrMod","Simulated_Individual"))
  
  Simulated_data <- foreach(i=seq(1,Nb_Individual),.combine="rbind",.errorhandling = "remove",.packages = c("parallel","foreach","doSNOW","plyr","mvtnorm"))%dopar%{
    ind_simulated_data <- Model_prediction_Function(fixed_params=Simulated_parameters,random_effects=Simulated_RE[i,],time=Time_simulation,sd_Err=Simulated_sd_ErrMod,covariate=Simulated_Individual[i,])
    ind_simulated_data$SEX <- ifelse(ind_simulated_data$SEX == 1,"Woman","Man")
    ind_simulated_data$New_Time <- ind_simulated_data$Time
    ind_simulated_data$Time <- ind_simulated_data$Time+28
    return(ind_simulated_data)
  }
  stopCluster(Cluster)
  # ----- #
  
  
  # -- Calculation of the probability to reach an Ab threshold at day 365 ####
  
  print("Probability estimation")
  Profile_Vac <- unique(Profiles[which(Profiles$ARM == 1),c("ARM","CATAGE","SEX")])
  
  Simulated_data$Ab <- 10^Simulated_data$log10Ab
  Simulated_data_Vac <- subset(Simulated_data,ARM == 1)
  Simulated_data_Vac <- Simulated_data_Vac[,c("ID","SEX","CATAGE","Time","log10Ab","Ab")]
  
  NCores <- detectCores()-1
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  
  clusterExport(Cluster,c("Profile_Vac","Ab_thresholds","Simulated_data_Vac","Proba_Above_AbThres_Function"))
  
  Probability_results <- foreach(i=seq(1,length(Ab_thresholds)),.errorhandling = "remove",.packages = c("parallel","foreach","doSNOW","plyr","mvtnorm"))%dopar%{
    
    Proba_Profile_Vs_Thresh <- Profile_Vac
    Proba_Profile_Vs_Thresh <- cbind(Proba_Profile_Vs_Thresh,Ab_thresh=Ab_thresholds[i],stringsAsFactors=F)
    Proba_Profile_Vs_Thresh$Proba_D365 <- sapply(seq(1,nrow(Proba_Profile_Vs_Thresh)),function(j,profile,data,time){
      tmp_data <- subset(data,CATAGE == as.character(profile$CATAGE[j]) & SEX == as.character(profile$SEX[j]))
      proba <- Proba_Above_AbThres_Function(data=tmp_data,time=time,thres=profile$Ab_thresh[j])
      return(proba)
    },data=Simulated_data_Vac,profile=Proba_Profile_Vs_Thresh,time=365)
    return(Proba_Profile_Vs_Thresh)
  }
  stopCluster(Cluster)
  
  Probability_results <- cbind(Simu=n,Probability_results)
  write.table(Probability_results,file = paste(Folder_results,paste("rVSV_proba_Simu",n,".csv",sep=""),sep="/"),row.names=FALSE,sep="\t",dec=".")
  # ----- #
  

  # -- Calculation of the time requested to reach an Ab threshold ####
  
  print("Time estimation")
  # -- Estimation of the time to reach an Ab threshold
  NCores <- detectCores()-1
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  clusterExport(Cluster,c("Profile_Vac","Ab_thresholds","Time_AbTrhreshold_Estimation_Function","Simulated_parameters","Simulated_RE","Simulated_sd_ErrMod","Simulated_Individual"))
  
  Estimated_time <- foreach(i=seq(1,length(Ab_thresholds)),.combine="rbind",.errorhandling = "remove",.packages = c("parallel","foreach","doSNOW","plyr","mvtnorm"))%dopar%{
    
    Ab_threshold <- Ab_thresholds[i]
    Log10_Abthresold <- log10(Ab_threshold)
    estimated_time <- Time_AbTrhreshold_Estimation_Function(Ab_thresh=Log10_Abthresold,fixed_params=Simulated_parameters,random_effects=Simulated_RE,sd_Err=0,covariate=Simulated_Individual)
    return(estimated_time)
  }
  stopCluster(Cluster)
  
  Median_Time <- ddply(.data=Estimated_time,.variables = .(ARM,CATAGE,SEX,Ab_thresh),summarise,Time=median(time,na.rm=TRUE))
  Median_Time <- cbind(Simu=n,Median_Time)
  write.table(Median_Time,file = paste(Folder_results,paste("rVSV_Time_Simu",n,".csv",sep=""),sep="/"),row.names=FALSE,sep="\t",dec=".")
  # ----- #
}


# ----------------------- #
# OBJECTIVE: Generate results to estimate both the time requested to reach a given Ab threshold and the probability to reach it at day 365
#            Focus on the vaccine regimen Ad26
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
Simulation_Random_Effects_Function <- function(Var_Var_Intercept, Sd_Intercept,nb_ind){
  
  # Simulation of intercept random effect
  simulated_var_RE <- as.numeric(rnorm(1,mean=c(Sd_Intercept^2), sqrt(Var_Var_Intercept)))
  simulated_RE   <- as.numeric(rnorm(nb_ind, mean=0, sqrt(simulated_var_RE)))
  return(simulated_RE)
}
Model_prediction_Function <- function(fixed_params,random_effects,sd_Err,covariate){
  
  with(as.list(c(fixed_params)),{
    
    dfsim <- covariate
    dfsim$Gam <- random_effects
    dfsim$ARM <- as.numeric(dfsim$ARM)
    dfsim$SEX <- ifelse(dfsim$SEX == "Man", 0, 1)
    dfsim$CATAGE1 <- ifelse(dfsim$CATAGE == "1-4", 1, 0)
    dfsim$CATAGE5 <- ifelse(dfsim$CATAGE == "5-11", 1, 0)
    dfsim$CATAGE12 <- ifelse(dfsim$CATAGE == "12-17", 1, 0)
    dfsim$SIERRA <- ifelse(dfsim$COUNTRY == "Sierra Leone", 1, 0)
    dfsim$LIBERIA <- ifelse(dfsim$COUNTRY == "Liberia", 1, 0)
    dfsim$GUINEA <- ifelse(dfsim$COUNTRY == "Guinea", 1, 0)
    
    Time_simulationa <- c(0,90,90)
    Time_simulationb <- c(0,0,185)
    
    Model_prediction <- NULL
    nb_ind <- nrow(dfsim)
    for(i in 1:nb_ind){
      model_pred_ind <- data.frame(dfsim[i,], Time=Time_simulation, VISITa=Time_simulationa, VISITb=Time_simulationb)
      
      model_pred_ind$log10Ab = B0 + model_pred_ind$Gam + B1*model_pred_ind$VISITa + B2*model_pred_ind$VISITb +
        B3*model_pred_ind$ARM + B4*model_pred_ind$CATAGE1 + B5*model_pred_ind$CATAGE5 + B6*model_pred_ind$CATAGE12 +
        B7*model_pred_ind$CATAGE1*model_pred_ind$ARM + B8*model_pred_ind$CATAGE5*model_pred_ind$ARM + B9*model_pred_ind$CATAGE12*model_pred_ind$ARM +
        B10*model_pred_ind$GUINEA + B11*model_pred_ind$SIERRA + B12*model_pred_ind$LIBERIA +
        B13*model_pred_ind$GUINEA*model_pred_ind$ARM + B14*model_pred_ind$SIERRA*model_pred_ind$ARM + B15*model_pred_ind$LIBERIA*model_pred_ind$ARM +
        B16*model_pred_ind$ABBn + B17*model_pred_ind$ABBn*model_pred_ind$ARM +
        B18*model_pred_ind$ARM*model_pred_ind$VISITa +
        B19*model_pred_ind$CATAGE1*model_pred_ind$VISITa + B20*model_pred_ind$CATAGE5*model_pred_ind$VISITa + B21*model_pred_ind$CATAGE12*model_pred_ind$VISITa +
        B22*model_pred_ind$ARM*model_pred_ind$CATAGE1*model_pred_ind$VISITa + B23*model_pred_ind$ARM*model_pred_ind$CATAGE5*model_pred_ind$VISITa + B24*model_pred_ind$ARM*model_pred_ind$CATAGE12*model_pred_ind$VISITa +
        B25*model_pred_ind$GUINEA*model_pred_ind$VISITa + B26*model_pred_ind$SIERRA*model_pred_ind$VISITa + B27*model_pred_ind$LIBERIA*model_pred_ind$VISITa +
        B28*model_pred_ind$ARM*model_pred_ind$GUINEA*model_pred_ind$VISITa + B29*model_pred_ind$ARM*model_pred_ind$SIERRA*model_pred_ind$VISITa + B30*model_pred_ind$ARM*model_pred_ind$LIBERIA*model_pred_ind$VISITa +
        B31*model_pred_ind$ABBn*model_pred_ind$VISITa + B32*model_pred_ind$ARM*model_pred_ind$ABBn*model_pred_ind$VISITa +
        B33*model_pred_ind$ARM*model_pred_ind$VISITb +
        B34*model_pred_ind$CATAGE1*model_pred_ind$VISITb + B35*model_pred_ind$CATAGE5*model_pred_ind$VISITb + B36*model_pred_ind$CATAGE12*model_pred_ind$VISITb +
        B37*model_pred_ind$ARM*model_pred_ind$CATAGE1* model_pred_ind$VISITb + B38*model_pred_ind$ARM*model_pred_ind$CATAGE5* model_pred_ind$VISITb + B39*model_pred_ind$ARM*model_pred_ind$CATAGE12*model_pred_ind$VISITb +
        B40*model_pred_ind$GUINEA*model_pred_ind$VISITb + B41*model_pred_ind$SIERRA*model_pred_ind$VISITb + B42*model_pred_ind$LIBERIA*model_pred_ind$VISITb +
        B43*model_pred_ind$ARM*model_pred_ind$GUINEA*model_pred_ind$VISITb + B44*model_pred_ind$ARM*model_pred_ind$SIERRA*model_pred_ind$VISITb + B45*model_pred_ind$ARM*model_pred_ind$LIBERIA*model_pred_ind$VISITb +
        B46*model_pred_ind$ABBn*model_pred_ind$VISITb + B47*model_pred_ind$ARM*model_pred_ind$ABBn*model_pred_ind$VISITb
      
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
Time_AbTrhreshold_Estimation_Function <- function(fixed_params,random_effects,covariate,sd_Err,Ab_thresh,breaking_point=90){
  
  with(as.list(c(fixed_params)),{
    
    dfsim <- covariate
    dfsim$Gam <- random_effects
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
      # j <- 1
      time_pred_ind <- data.frame(dfsim[j,],Ab_thresh=Ab_thresh)
      estimated_time <- optimize(function(time) abs((B0 + time_pred_ind$Gam + B1*min(time,breaking_point) + B2*max(time-breaking_point,0) +
                                                       B3*time_pred_ind$ARM + B4*time_pred_ind$CATAGE1 + B5*time_pred_ind$CATAGE5 + B6*time_pred_ind$CATAGE12 +
                                                       B7*time_pred_ind$CATAGE1*time_pred_ind$ARM + B8*time_pred_ind$CATAGE5*time_pred_ind$ARM + B9*time_pred_ind$CATAGE12*time_pred_ind$ARM +
                                                       B10*time_pred_ind$GUINEA + B11*time_pred_ind$SIERRA + B12*time_pred_ind$LIBERIA +
                                                       B13*time_pred_ind$GUINEA*time_pred_ind$ARM + B14*time_pred_ind$SIERRA*time_pred_ind$ARM + B15*time_pred_ind$LIBERIA*time_pred_ind$ARM +
                                                       B16*time_pred_ind$ABBn + B17*time_pred_ind$ABBn*time_pred_ind$ARM +B18*time_pred_ind$ARM*min(time,breaking_point) +
                                                       B19*time_pred_ind$CATAGE1*min(time,breaking_point) + B20*time_pred_ind$CATAGE5*min(time,breaking_point) + B21*time_pred_ind$CATAGE12*min(time,breaking_point) +
                                                       B22*time_pred_ind$ARM*time_pred_ind$CATAGE1*min(time,breaking_point) + B23*time_pred_ind$ARM*time_pred_ind$CATAGE5*min(time,breaking_point) + B24*time_pred_ind$ARM*time_pred_ind$CATAGE12*min(time,breaking_point) +
                                                       B25*time_pred_ind$GUINEA*min(time,breaking_point) + B26*time_pred_ind$SIERRA*min(time,breaking_point) + B27*time_pred_ind$LIBERIA*min(time,breaking_point) +
                                                       B28*time_pred_ind$ARM*time_pred_ind$GUINEA*min(time,breaking_point) + B29*time_pred_ind$ARM*time_pred_ind$SIERRA*min(time,breaking_point) + B30*time_pred_ind$ARM*time_pred_ind$LIBERIA*min(time,breaking_point) +
                                                       B31*time_pred_ind$ABBn*min(time,breaking_point) + B32*time_pred_ind$ARM*time_pred_ind$ABBn*min(time,breaking_point) +
                                                       B33*time_pred_ind$ARM*max(time-breaking_point,0) +
                                                       B34*time_pred_ind$CATAGE1*max(time-breaking_point,0) + B35*time_pred_ind$CATAGE5*max(time-breaking_point,0) + B36*time_pred_ind$CATAGE12*max(time-breaking_point,0) +
                                                       B37*time_pred_ind$ARM*time_pred_ind$CATAGE1* max(time-breaking_point,0) + B38*time_pred_ind$ARM*time_pred_ind$CATAGE5* max(time-breaking_point,0) + B39*time_pred_ind$ARM*time_pred_ind$CATAGE12*max(time-breaking_point,0) +
                                                       B40*time_pred_ind$GUINEA*max(time-breaking_point,0) + B41*time_pred_ind$SIERRA*max(time-breaking_point,0) + B42*time_pred_ind$LIBERIA*max(time-breaking_point,0) +
                                                       B43*time_pred_ind$ARM*time_pred_ind$GUINEA*max(time-breaking_point,0) + B44*time_pred_ind$ARM*time_pred_ind$SIERRA*max(time-breaking_point,0) + B45*time_pred_ind$ARM*time_pred_ind$LIBERIA*max(time-breaking_point,0) +
                                                       B46*time_pred_ind$ABBn*max(time-breaking_point,0) + B47*time_pred_ind$ARM*time_pred_ind$ABBn*max(time-breaking_point,0) + rnorm(1,mean=0,sd=sd_Err)) - Ab_thresh),interval = c(0,1e10),tol=1e-6)
      time_pred_ind$time <- estimated_time$minimum
      time_pred_ind$objective_time <- estimated_time$objective
      
      Time_prediction <- rbind(Time_prediction,time_pred_ind)
    }
    return(Time_prediction)
  })
}
# --------------------- #



# ----- Arguments for calculation server ----- #
# See the batch file "Simulation_results_Ad26.sh" for an example of run
args <- commandArgs(trailingOnly = TRUE)             # List of arguments called in the batch file 
slurm_arrayid <- Slurm_env(x='SLURM_ARRAY_TASK_ID')  # calculation was splitted into multiple jobs
arr <- as.numeric(slurm_arrayid)
# --------------------- #





# ****************************************************** ####
#       PART 1: Estimation of the MEM on raw data        ####
# ***************************************************** ####
# -- Load of data ####
load("Data/dfprevacAD26.Rda")
dfprevac$VISITb[which(dfprevac$VISIT == 365 & dfprevac$VISITb == 180)] <- 185
breaking_timepoint <- 90  # time post-peak
# ----- #

# -- Estimation of the Mixed effects model ####
MEM_Ad26 <- lmer(TEST_MED_FINAL_LOG ~ VISITa*ARM*CATAGE + VISITa*ARM*COUNTRYID + VISITa*ARM*ABBn + VISITb*ARM*CATAGE + VISITb*ARM*COUNTRYID + VISITb*ARM*ABBn +  (1 | SUBJID_N) , dfprevac, REML = T)

# Fixed parameters
Fixed_parameters_order <- c("B0", "B1", "B3", "B4", "B6", "B5", "B10", "B12", "B11", "B16", "B2", "B18", "B19", "B21", "B20",
                            "B7", "B9", "B8", "B25", "B27", "B26", "B13", "B15", "B14", "B31", "B17", "B33", "B34", "B36", "B35",
                            "B40", "B42", "B41", "B46", "B22", "B24", "B23", "B28", "B30", "B29", "B32", "B37", "B39", "B38",
                            "B43", "B45", "B44", "B47")
MEM_Fixed_parameters <-  summary(MEM_Ad26)$coefficients[,1]
names(MEM_Fixed_parameters) <- Fixed_parameters_order
MEM_VCoV_Mat <- as.matrix(vcov(MEM_Ad26))
colnames(MEM_VCoV_Mat) <- Fixed_parameters_order ; rownames(MEM_VCoV_Mat) <- Fixed_parameters_order
# Random effects
Sd_Intercept <- attr(VarCorr(MEM_Ad26)$SUBJID_N, "stddev")[1]
Var_Var_Intercept <- 9.024945e-06    # Values extracted from the estimated covariance matrix
# # Error model
Sd_Residual <- attr(VarCorr(MEM_Ad26),"sc")
Res_Var <- 3.638740e-06   # Values extracted from the estimated covariance matrix
# ----- #

# -- Extraction of the database required for the simulations
Data <- dfprevac[,c("SUBJID_N","VISIT","COUNTRYID","ARM","SEX","TEST_MED_FINAL_LOG","CATAGE","ABBn","VISITa","VISITb")]
colnames(Data) <- c("ID","TIME","COUNTRY","ARM","SEX","LOG10Ab","CATAGE","ABBn","TIMEa","TIMEb")
Data$ABBnInf200 <- 1*(Data$ABBn < log10(200))
# ----- #

# -- Estimation of the proportion of each type of profile in the population
Profiles <- ddply(.data=Data,.variables = .(COUNTRY,ARM,CATAGE,SEX,ABBnInf200),summarise,
                  N=length(unique(ID)),Proportion = length(unique(ID))/length(unique(Data$ID)),
                  Mean_Abbn = mean(ABBn,na.rm=T),sd_Abbn = sd(ABBn,na.rm=T))
# ----- #





# ****************************************** ####
#       PART 2: Generation of results        ####
# **************************************** ####
Folder_results <- "Rcode/Results/Ad26_2"
dir.create(Folder_results,recursive = TRUE)

Nb_Individual <- as.numeric(args[1])
Ab_thresholds <- seq(100,5000,by=100)



simus <- ((arr-1)*10+1):(arr*10)

for(n in simus){
  
  
  # -- Data simulation ####
  
  print("-- Simulation of individuals")
  # --Simulation of individuals
  Simulated_Individual <- Simulation_database_Function(obs_data = Data,profiles=Profiles,nb_ind=Nb_Individual)
  
  # - Simulation of model parameters
  # Time of simulation
  Time_simulation <- sort(unique(Data$TIME))
  
  print("-- Simulation of population parameters")
  # Fixed parameters
  Simulated_parameters <- Simulation_Population_Parameters_Function(parameters=MEM_Fixed_parameters,varcov_mat=MEM_VCoV_Mat)[1,]
  
  print("-- Simulation of random effects")
  # Random Effects
  Simulated_RE <- Simulation_Random_Effects_Function(Var_Var_Intercept = Var_Var_Intercept, Sd_Intercept = Sd_Intercept,nb_ind=Nb_Individual)
  
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
    ind_simulated_data <- Model_prediction_Function(fixed_params=Simulated_parameters,random_effects=Simulated_RE[i],sd_Err=Simulated_sd_ErrMod,covariate=Simulated_Individual[i,])
    ind_simulated_data$SEX <- ifelse(ind_simulated_data$SEX == 1,"Woman","Man")
    ind_simulated_data$Time <- ind_simulated_data$Time - 90  # Time post-peak
    ind_simulated_data$breakingPoint <- 90 # Time post-peak
    return(ind_simulated_data)
  }
  stopCluster(Cluster)
  # ----- #

  
  
  # -- Calculation of the probability to reach an Ab threshold at day 365 ####
  
  print("Probability estimation")
  Profile_Vac <- unique(Profiles[which(Profiles$ARM == 1),c("ARM","CATAGE")])
  
  Simulated_data$Ab <- 10^Simulated_data$log10Ab
  Simulated_data_Vac <- subset(Simulated_data,ARM == 1)
  Simulated_data_Vac <- Simulated_data_Vac[,c("ID","CATAGE","Time","log10Ab","Ab")]
  
  NCores <- detectCores()-1
  Cluster <- parallel::makeCluster(NCores,type="SOCK")
  registerDoSNOW(Cluster)
  
  clusterExport(Cluster,c("Profile_Vac","Ab_thresholds","Simulated_data_Vac","Proba_Above_AbThres_Function"))
  
  # i <- 1
  Probability_results <- foreach(i=seq(1,length(Ab_thresholds)),.errorhandling = "remove",.packages = c("parallel","foreach","doSNOW","plyr","mvtnorm"))%dopar%{
    
    Proba_Profile_Vs_Thresh <- Profile_Vac
    Proba_Profile_Vs_Thresh <- cbind(Proba_Profile_Vs_Thresh,Ab_thresh=Ab_thresholds[i],stringsAsFactors=F)
    Proba_Profile_Vs_Thresh$Proba_D365 <- sapply(seq(1,nrow(Proba_Profile_Vs_Thresh)),function(j,profile,data,time){
      tmp_data <- subset(data,CATAGE == as.character(profile$CATAGE[j]))
      proba <- Proba_Above_AbThres_Function(data=tmp_data,time=time-90,thres=profile$Ab_thresh[j])
      return(proba)
    },data=Simulated_data_Vac,profile=Proba_Profile_Vs_Thresh,time=365)
    return(Proba_Profile_Vs_Thresh)
  }
  stopCluster(Cluster)
  
  Probability_results <- cbind(Simu=n,Probability_results)
  write.table(Probability_results,file = paste(Folder_results,paste("Ad26MVA_proba_Simu",n,".csv",sep=""),sep="/"),row.names=FALSE,sep="\t",dec=".")
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
  
  Median_Time <- ddply(.data=Estimated_time,.variables = .(ARM,CATAGE,Ab_thresh),summarise,Time=median(time,na.rm=TRUE))
  Median_Time <- cbind(Simu=n,Median_Time)
  write.table(Median_Time,file = paste(Folder_results,paste("Ad26_Time_Simu",n,".csv",sep=""),sep="/"),row.names=FALSE,sep="\t",dec=".")
  # ----- #
}
  
  
  




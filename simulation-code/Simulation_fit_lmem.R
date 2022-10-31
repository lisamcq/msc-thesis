#Simulation study: estimate the lme model 
library(nlme)
library(dplyr)
library(tidyr)
library(mvtnorm)

#read in simulated values of explanatory variables
sim_pred = read.csv(file = "Simulated_explanatories_v2.csv",
                    header = TRUE)
sim_pred = sim_pred %>% dplyr::select(-X)


#set values of model coefficients
beta_time = 0.479
beta0 = 6.038
beta_cancer = -0.35
r = 0.65
sd_error = 2.5
sd_random = 2.4

num_timepts = 43
num_postdiag = 12
#cancer effect sizes to test: -0.35, -0.1
#sd_error to test: 2.5
#sd_random to test: 2.4 

#function for the variance matrix of the autoregressive error model
Sigma <- function(t, r, var_err){
  corrmat = matrix(nrow=t, ncol = t)
  for(i in 1:t){
    for(j in 1:t){
      corrmat[i, j] = r^(abs(i-j))
    }
  }
  return(var_err*corrmat)
}

#function for simulating from the autoregressive error model
simulation_autoerror <- function(N = samplesize, timepoints = num_timepts, num_pts_postdiag = num_postdiag,
                                 intercept = beta0, cancer_effect = beta_cancer, 
                                 time_effect = beta_time,
                                 rho = r, sd_e = sd_error, sd_r = sd_random){
  #initialize columns for each year with NA values
  columnnames = paste("Year", 0:(timepoints-1), sep = "")
  simdata = sim_pred[1:N,]
  simdata[,columnnames] <- NA
  
  time_vec= 0:42
  cancer_vec =  c(rep(0, timepoints-num_pts_postdiag), rep(1, num_pts_postdiag))
  #for each subject
  for(i in 1:N){
    #vector of errors for all time points including baseline
    sig = Sigma(t=timepoints, r=rho, var_err = sd_e^2)
    error = as.vector(rmvnorm(1,mean=rep(0, timepoints), 
                              sigma =sig))
    #subject-specific effect
    random_int = rnorm(1, 0, sd_r)
    #time vector
    
    beta_vec = c(intercept, cancer_effect, time_effect)
    X = matrix(c(rep(1, timepoints),  #intercept
                 cancer_vec,
                 #cancer indicator causing effect post-diagnosis
                 time_vec), nrow = timepoints, ncol = length(beta_vec))
    
    b_vec = rep(random_int, timepoints)
    simdata[i,columnnames]= X%*%beta_vec + b_vec + error
    #for(j in columnnames){
    #   if(simdata[i, j] < 0) simdata[i,j] =0
    # }
  }
  return(simdata)
}

#function for estimating the lme model
fit_LMEM <- function(simdata, num_timepoints = num_timepts){
  simdata_long <- pivot_longer(data=simdata, cols=paste("Year", 0:(num_timepoints-1), sep = ""), names_to = "Year_since18", 
                               values_to = "transf_income", names_prefix = "Year")
  simdata_long$Year_since18 = as.integer(simdata_long$Year_since18)
  simdata_long = simdata_long %>% mutate(Cancer_timevary = ifelse(Year_since18 < Age_diag-18, 0, 1))
  mod = lme(transf_income ~ Year_since18 + Cancer_timevary, data=simdata_long, 
            random = ~1|SubjectID,
            control = lmeControl(opt = "optim"))
  beta_int = unname(fixed.effects(mod)["(Intercept)"])
  se_beta_int = summary(mod)$tTable["(Intercept)",2]
  beta_c = unname(fixed.effects(mod)["Cancer_timevary"])
  se_beta_c = summary(mod)$tTable["Cancer_timevary",2]
  beta_t = unname(fixed.effects(mod)["Year_since18"])
  se_beta_t = summary(mod)$tTable["Year_since18",2]
  
  sd_random_est = sqrt(getVarCov(mod)[1,1])
  sd_error_est = sqrt(getVarCov(mod, type = "conditional")[[1]][1,1])
  return(c(beta_int, se_beta_int, beta_c, se_beta_c, 
           beta_t, se_beta_t,
           sd_random_est, sd_error_est))
  
}


#function that simulates data from autoregressive error model and fits lme model
model_resultsLMEM <- function(B, N, timepoints, intercept = beta0,
                          cancer_effect = beta_cancer, time_effect = beta_time,
                          r_errormod = r, sd_e = sd_error, sd_r = sd_random){
  model_outputs = matrix(nrow = B, ncol = 8)
  colnames(model_outputs) = c("Beta0_hat", "Beta0 SE", "Cancer effect", 
                              "Cancer SE","Time effect", "Time SE",
                              "sigma_u_hat", "sigma_hat")
  
  for(i in 1:B){
    simdata = simulation_autoerror(N = N, timepoints = timepoints,
                                   intercept = intercept, cancer_effect = cancer_effect,
                                   time_effect = time_effect,
                                   rho=r_errormod, sd_e = sd_e, sd_r = sd_r)
    model_outputs[i,] = fit_LMEM(simdata = simdata, num_timepoints = timepoints)
  }
  return(model_outputs)
}
#model_resultsLMEM(B=5, N = 1000, timepoints=num_timepts)

#values for simulation
B=100
#B=50
#samplesizes = c(50, 250, 1000)
samplesizes =c(50, 1000)
#effectsizes = c(log(0.9), log(0.7))
effectsizes =c(-0.1, -0.35)
#errorsds = c(0.1, 0.25)
errorsds = c(1, 2.5)
#randomsds = c(0.3, 0.6)
randomsds = c(1, 2.5)

list_model_results=list(length=length(samplesizes)*length(effectsizes)*length(errorsds)*length(randomsds))
scenario_names = vector(length=length(samplesizes)*length(effectsizes)*length(errorsds)*length(randomsds))
set.seed(10000)
i=1
#for each value, run simulation
for (N in samplesizes) {
  for(effsize in effectsizes){
    for(errorsd in errorsds){
      for(randomsd in randomsds){
        scenario_names[i] = paste(N,round(effsize,2), errorsd, randomsd)
        list_model_results[[i]] = model_resultsLMEM(B = B, N = N, timepoints = num_timepts, 
                                                  intercept = beta0, cancer_effect = effsize,
                                                  time_effect = beta_time, 
                                                  r_errormod = r,
                                                  sd_e = errorsd, sd_r = randomsd)
        
        i=i+1
      }
    }
    
  }
}

names(list_model_results) = scenario_names
#output results
saveRDS(list_model_results, file = "LMEM_sim2_results_v6.rds")


library(dplyr)
library(tidyr)
library(mvtnorm)
library(numDeriv)
library(snowfall)
library(MASS)

#function for the matrix B in the autoregressive response model formula
Bmatrix <- function(t, f_rho){
  #rho = coefficient
  #t = number of time points including baseline
  Back = matrix(nrow = t, ncol = t)
  rho_vec = vector()
  for (i in 0:(t-1)) {
    rho_vec = c(rho_vec,((exp(f_rho)-1)/(1+exp(f_rho)))^i)
  }
  Back[,1] = rho_vec
  for (j in 2:t) {
    Back[,j] = c(0,Back[1:t-1,j-1])
  }
  return(Back)
}

#function for the variance matrix of the autoregressive response model
variance_matrix <- function(n, log_var_random, log_var_error){
  var_matrix = matrix(nrow=n, ncol=n)
  var_matrix[upper.tri(var_matrix, diag = F)] = exp(log_var_random)
  var_matrix[lower.tri(var_matrix, diag = F)] = exp(log_var_random)
  diag(var_matrix) = exp(log_var_random) + exp(log_var_error)
  return(var_matrix)
}

#function for the log-likelihood of the autoregressive response model
log_like_full <- function(f_beta_rho_var, X, n, Y, p){
  beta = f_beta_rho_var[1:p]
  f_rho = f_beta_rho_var[p+1]
  log_sigma1 = f_beta_rho_var[p+2]
  log_sigma2 = f_beta_rho_var[p+3]
  N = length(X)
  Back = lapply(n, Bmatrix, f_rho=f_rho) #calculate B_i for each n_i
  var_list = lapply(n, variance_matrix, log_var_random=log_sigma1, log_var_error=log_sigma2) 
  #calculate V_i for each n_i
  Sig = lapply(1:N, function(i) Back[[i]]%*%var_list[[i]]%*%t(Back[[i]])) #calculate \Sigma_i for each i
  log_density = sapply(1:N, function(i) mvtnorm::dmvnorm(Y[[i]], mean = Back[[i]]%*%X[[i]]%*%beta, sigma = Sig[[i]], log = TRUE, checkSymmetry = FALSE))
  return(-2*sum(log_density))
}

#pull in values of the explanatory variables for simulation
sim_pred = read.csv(file = "Simulated_explanatories_v2.csv",
                    header = TRUE)
sim_pred = sim_pred %>% dplyr::select(-X)
#set true values for the model coefficients
beta_time = 0.479
beta0 = 6.038
beta_cancer = -0.35
r = 0.65
sd_error = 2.5
sd_random = 2.5

#set number of time points (total and post-diagnosis)
num_timepts = 43
num_postdiag = 12

#set sample size
samplesize = 1000

#function of the variance of the autoregressive error model
Sigma <- function(t, r, var_err){
  corrmat = matrix(nrow=t, ncol = t)
  for(i in 1:t){
    for(j in 1:t){
      corrmat[i, j] = r^(abs(i-j))
    }
  }
  return(var_err*corrmat)
}

#function for simulating data from the autoregressive error model
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


#function for estimating the autoregressive response model
fit_autoresponse <- function(simdata, timepoints = num_timepts,
                             startpts_beta,
                             startpt_rho, startpts_var){
  simdata_long <- pivot_longer(data=simdata, cols=paste("Year", 0:(timepoints-1), sep = ""), names_to = "Year_since18", 
                               values_to = "transf_income", names_prefix = "Year")
  simdata_long$Year_since18 = as.integer(simdata_long$Year_since18)
  simdata_long = simdata_long %>% mutate(Cancer_timevary = ifelse(Year_since18 < Age_diag-18, 0, 1))
  
  simdata_long$Intercept = 1
  X_sim = list()
  Y_sim = list()
  subjectIDs = unique(simdata_long$SubjectID)
  n = vector()
  for(j in 1:length(subjectIDs)){
    X_sim[[j]] = as.matrix(simdata_long %>% filter(SubjectID==subjectIDs[j]) %>% dplyr::select(Intercept,
                                                                                               Cancer_timevary, 
                                                                                               Year_since18))
    Y_sim[[j]] = as.vector(as.matrix(simdata_long %>% filter(SubjectID==subjectIDs[j]) %>% dplyr::select(transf_income)))
    n[j] = sum(simdata_long$SubjectID == subjectIDs[j])
    #dim(X_sim[[j]]) = c(n[j], num_explanatory)
  }
  #optimize the log-likelihood
  op = optim(fn = log_like_full, par = c(startpts_beta, startpt_rho, startpts_var), X=X_sim, 
             n= n, Y=Y_sim, p =3,
              method ="BFGS",
              hessian = TRUE)
    return(op)
  
}
#fit_autoresponse(simdata,
# startpts_beta = c(8, -0.2, 0.2), startpt_rho = 0.3, startpts_var = c(3,3))

#vectors of parameter values for simulation
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

#function for running simulation in parallel
paraLapply = function(arr,fn,cores)
{
  sfInit(parallel=TRUE,cpus=cores) ## change to 4,8 w/e
  
  fun = function(){x = sfLapply( x=arr, fun=fn ) ; x}
  sfExportAll()
  
  res = fun() ;
  
  sfRemoveAll( ) 	## If you break the function (using control+C), you should
  sfStop() 		## run these commands unless you really like having zero RAM (hardcore memory leakzz)
  res ;
}

#input to the parallel function - simulates data and then fits autoregressive response model
parallel_replicates <- function(b){
  library(mvtnorm)
  library(dplyr)
  library(tidyr)
  simdata = simulation_autoerror(N = ss, timepoints = num_timepts, intercept = beta0, 
                                 cancer_effect = effsize,
                                 time_effect = beta_time,
                                 rho=r, sd_e = errorsd, sd_r = randomsd)

  return(fit_autoresponse(simdata = simdata,
                          timepoints = num_timepts, startpts_beta = initialbeta,
                          startpt_rho = initialrho, startpts_var = initialvar))
}


#test = paraLapply(1:2, parallel_replicates,2)
set.seed(123)
#initial values for the model estimation algorithm
initialrho = r + runif(1, -0.3, 0.3)
initialvarrandom = c(1 + runif(1,-0.5, 0.5), 2.5 + runif(1,-1, 1))
initialvarerror = c(1 + runif(1,-0.5, 0.5), 2.5 + runif(1,-1, 1))
initialbeta2 = c(beta0 + runif(1,-2, 2), -0.1 + runif(1, -0.05, 0.05), -0.35 + runif(1, -0.1, 0.1), beta_time + runif(1, -0.2, 0.2))

set.seed(10000)
#set up lists to save model results in
list_model_est=list(length=length(samplesizes)*length(effectsizes)*length(errorsds)*length(randomsds))
list_model_vars=list(length=length(samplesizes)*length(effectsizes)*length(errorsds)*length(randomsds))
scenario_names = vector(length=length(samplesizes)*length(effectsizes)*length(errorsds)*length(randomsds))
i=1
#for each value of each parameter, run simulation and save results
for(ss in samplesizes){
  for(effsize in effectsizes){
    for(errorsd in errorsds){
      for(randomsd in randomsds){
        scenario_names[i] = paste(ss, round(effsize,2), errorsd, randomsd)
        print(scenario_names[i])
        if(effsize == -0.1) initialbeta = initialbeta2[-3]
        if(effsize == -0.35) initialbeta = initialbeta2[-2]
        if(errorsd == 1 & randomsd == 1) initialvar = log(c(initialvarrandom[1], initialvarerror[1]))
        if(errorsd == 2.5 & randomsd == 1) initialvar = log(c(initialvarrandom[1], initialvarerror[2]))
        if(errorsd == 1 & randomsd==2.5) initialvar = log(c(initialvarrandom[2], initialvarerror[1]))
        if(errorsd == 2.5 & randomsd == 2.5) initialvar = log(c(initialvarrandom[2], initialvarerror[2]))
        temp = paraLapply(1:B, parallel_replicates, 10)
        hessian_est = lapply(temp, function(i) i$hessian)
        params_est = lapply(temp, function(i) i$par)
        temp2 = matrix(unlist(params_est), nrow = B, byrow = TRUE)
        colnames(temp2) = c("Beta0_hat", "Cancer effect", 
                            "Time effect", "f_rho", "log_sigma_u^2", "log_sigma^2")
        list_model_est[[i]] = temp2
        list_model_vars[[i]] = hessian_est
        
        i=i+1
      }
    }
    
    
  }
}

names(list_model_est) = scenario_names
names(list_model_vars) = scenario_names
#output results
saveRDS(list_model_est, file = "Autoresponse_exp3_sim2_estimates_v6.rds")
saveRDS(list_model_vars, file = "Autoresponse_exp3_sim2_variance_v6.rds")


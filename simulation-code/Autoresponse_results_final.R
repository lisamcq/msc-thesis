#Calculate bias and coverage probability of estimates from autoregressive response error model

#read in results
est = readRDS(file = "Autoresponse_exp3_sim2_estimates_v6.rds")
hess =readRDS(file = "Autoresponse_exp3_sim2_variance_v6.rds")
#read in simulated values of explanatory variables
sim_pred = read.csv("Simulated_explanatories_v2.csv", header = T)

B=nrow(est[[1]])
est2 = lapply(est, function(i) cbind(i[,1:3], (exp(i[,4]) - 1)/(exp(i[,4])+1), sqrt(exp(i[,5:6]))))

#function for calculating the marginal effect estimate
marginaleffect <- function(beta_rho, t, t_d){
  beta = beta_rho[1]
  rho = beta_rho[2]
  return((1-rho^(t-t_d+1))/(1-rho)*beta)
}
marginaleffect(c(est2[[1]][1,2], est2[[1]][1,4]), t = 31, t_d=30)

#function for calculating the standard error of the marginal effect estimate
marginalse <- function(beta_rho, t,  t_d, var){
  g = grad(marginaleffect, x = beta_rho, t =t, t_d = t_d)
  var_sub = var[c(2,4), c(2,4)]
  return(sqrt(t(g)%*%(solve(0.5*var_sub))%*%g))
}
marginalse(c(est2[[1]][1,2], est2[[1]][1,4]), t=1, t_d=0, var = hess[[1]][[1]])

#coverage probability of conditional effect
beta_cancer_vec = rep(c(rep(-0.1, 4), rep(-0.35,4)),2)
cov_prob_cond = vector(length = length(est2))
for(i in 1:length(est2)){
  temp = vector()
  for(j in 1:B){
    temp[j] = est2[[i]][j, "Cancer effect"] - qnorm(0.975)*sqrt(solve(0.5*hess[[i]][[j]])[2,2]) <= beta_cancer_vec[i] &
      est2[[i]][j, "Cancer effect"] + qnorm(0.975)*sqrt(solve(0.5*hess[[i]][[j]])[2,2]) >= beta_cancer_vec[i]
  } 
  cov_prob_cond[i] = sum(temp)/B
}
names(cov_prob_cond) = names(est2)

est2[[16]][, "Cancer effect"] - qnorm(0.975)*sapply(hess[[16]], function(i) sqrt(solve(0.5*i)[2,2]))

#bias of conditional effect
avgbeta = sapply(est2, function(i) mean(i[,"Cancer effect"]))
avgbias = avgbeta - beta_cancer_vec
saveRDS(avgbias, "Bias_conditionalARR.rds")

#sd of bias
sdbiasARR_cond = sapply(1:16, function(i) sd(est2[[i]][,"Cancer effect"] - beta_cancer_vec[i]))
saveRDS(sdbiasARR_cond, file = "SD_Bias_ARR_cond.rds")

#empirical sd of beta^hat
sdbetahatARR_cond = sapply(1:16, function(i) sd(est2[[i]][,"Cancer effect"]))
#average se of betahat
#se betahat
sebetahat = list(length=16)
for(i in 1:16){
  sebetahat[[i]] = vector(length=B)
  for(j in 1:B)
    sebetahat[[i]][j] = sqrt(solve(0.5*hess[[i]][[j]])[2,2])
}
avgsebetahatARR_cond = sapply(1:16, function(i) mean(sebetahat[[i]]))

saveRDS(sdbetahatARR_cond, file = "SD_betahat_ARR_cond.rds")
saveRDS(avgsebetahatARR_cond, file = "Avgse_betahat_ARR_cond.rds")

#CI for coverage prob
successARR_cond = cov_prob_cond*100
CI_ARR_cond_L = binom.confint(x = successARR_cond, n=100, methods = "ac")$lower
CI_ARR_cond_U = binom.confint(x = successARR_cond, n=100, methods = "ac")$upper

################
#Marginal effect
#one year post-diagnosis

cov_prob_marg1 = vector(length = length(est2))
avgbetamarg = vector(length = length(est2))
margeff = list()
margse = list()
for(i in 1:length(est2)){
  temp = vector(length=B)
  efftemp = vector(length=B)
  margeff[[i]] = vector(length=B)
  margse[[i]] = vector(length=B)
  for(j in 1:B){
    margeff[[i]][j] = marginaleffect(c(est2[[i]][j, "Cancer effect"], est2[[i]][j,4]), t=32,t_d=31)
    margse[[i]][j] = marginalse(c(est2[[i]][j, "Cancer effect"], est2[[i]][j,4]), t=32, t_d=31, var = hess[[i]][[j]])
    temp[j] = margeff[[i]][j] - qnorm(0.975)*margse[[i]][j] <= beta_cancer_vec[i] &
      margeff[[i]][j] + qnorm(0.975)*margse[[i]][j] >= beta_cancer_vec[i]
   
  } 
  avgbetamarg[i] = mean(margeff[[i]])
  cov_prob_marg1[i] = sum(temp)/B
}
names(cov_prob_marg1) = names(est2)
names(avgbetamarg) = names(est2)
#bias in marginal effect

avgbiasmarg = avgbetamarg - beta_cancer_vec
saveRDS(avgbiasmarg, file = "Bias_marginalARR.rds")

#sd of bias
sdbiasARR_marg = sapply(1:16, function(i) sd(margeff[[i]] - beta_cancer_vec[i]))
saveRDS(sdbiasARR_marg, file = "SD_Bias_ARR_marg.rds")

#empirical sd of beta^hat
sdbetahatARR_marg = sapply(1:16, function(i) sd(margeff[[i]]))
#average se of betahat
avgsebetahatARR_marg = sapply(1:16, function(i) mean(margse[[i]]))

saveRDS(sdbetahatARR_marg, file = "SD_betahat_ARR_marg.rds")
saveRDS(avgsebetahatARR_marg, file = "Avgse_betahat_ARR_marg.rds")


#Agresti confidence interval for coverage probability
library(binom)
successARR_marg = cov_prob_marg1*100
CI_ARR_marg_L = binom.confint(x = successARR_marg, n=100, methods = "ac")$lower
CI_ARR_marg_U = binom.confint(x = successARR_marg, n=100, methods = "ac")$upper
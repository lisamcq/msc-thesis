#Code for calculating the bias and coverage probability of estimates from the LME model

#read in results
model_resultsLME = readRDS(file = "LMEM_sim2_results_v6.rds")

B = 100

#calculate coverage probability
beta_cancer_vec = rep(c(rep(-0.1, 4), rep(-0.35,4)),2)
cov_prob_LME = vector(length = length(model_resultsLME))
for(i in 1:length(model_resultsLME)){
  cov_prob_LME[i] = sum(apply(as.matrix(model_resultsLME[[i]]), MARGIN = 1, FUN = function(j) 
    j["Cancer effect"] - qnorm(0.975)*j["Cancer SE"] <= beta_cancer_vec[i] & 
      j["Cancer effect"] + qnorm(0.975)*j["Cancer SE"] >= beta_cancer_vec[i]))/B
}
names(cov_prob_LME) = names(model_resultsLME)

saveRDS(cov_prob_LME, file ="LMEcov_prob_results.rds")

#average over N=50
mean(cov_prob[1:8]) #= 0.71875
##average over N = 1000
mean(cov_prob[9:16]) #= 0.7
#average over effect size of -0.1
mean(cov_prob[c(1:4,9:12)]) #= 0.71625
#average over effect size of -0.35
mean(cov_prob[c(5:8, 13:16)]) #= 0.7025
#average over sd_random = 1
mean(cov_prob[c(1,3,5,7,9,11,13,15)]) #= 0.72375
#average over sd_random=2.5
mean(cov_prob[c(2,4,6,8,10,12,14,16)]) #=0.695
#average over sd_error=1
mean(cov_prob[c(1,2,5,6,9,10,13,14)]) #= 0.7075
#average over sd_error = 2.5
mean(cov_prob[c(3,4,7,8,11,12,15,16)]) #=0.71125

###########################
#determining the average bias
#calculate average \beta\hat for each run
#bias is average betahat - true beta

avgbeta = sapply(model_resultsLME, function(i) mean(i[,"Cancer effect"]))
#names(avgbeta) = names(model_resultsLME)

avgbias = avgbeta - beta_cancer_vec
saveRDS(avgbias, file = "LMEbias_results.rds")
#average over N=50
mean(avgbias[1:8]) #=0.009513165
##average over N = 1000
mean(avgbias[9:16]) #= -0.0003121726
#average over effect size of -0.1
mean(avgbias[c(1:4,9:12)]) #= 0.009646276
#average over effect size of -0.35
mean(avgbias[c(5:8, 13:16)]) #= -0.0004452831
#average over sd_random = 1
mean(avgbias[c(1,3,5,7,9,11,13,15)]) #= 0.002540541
#average over sd_random=2.5
mean(avgbias[c(2,4,6,8,10,12,14,16)]) #= 0.006660452
#average over sd_error=1
mean(avgbias[c(1,2,5,6,9,10,13,14)]) #= -0.0001847093
#average over sd_error = 2.5
mean(avgbias[c(3,4,7,8,11,12,15,16)]) #=0.009385702

#average standard error of betahat
avgsebeta = sapply(model_resultsLME, function(i) mean(i[,"Cancer SE"]))
#calculate true SE of beta?

#sd of bias
sdbiasLME = sapply(1:16, function(i) sd(model_resultsLME[[i]][,"Cancer effect"] - beta_cancer_vec[i]))
saveRDS(sdbiasLME, file = "SD_Bias_LME.rds")

#SE plot
#sample sd of betahat
sdbetahatLME = sapply(1:16, function(i) sd(model_resultsLME[[i]][,"Cancer effect"]))
avgsebetahatLME = sapply(1:16, function(i) mean(model_resultsLME[[i]][,"Cancer SE"]))
saveRDS(sdbetahatLME, file = "SD_betahat_LME.rds")
saveRDS(avgsebetahatLME, file = "Avgse_betahat_LME.rds")


#CI for coverage prob
successLME = cov_prob_LME*100
CI_LME_L = binom.confint(x = successLME, n=100, methods = "ac")$lower
CI_LME_U = binom.confint(x = successLME, n=100, methods = "ac")$upper

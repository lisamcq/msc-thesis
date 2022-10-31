#Code for calculating the bias and coverage probability of estimates from the autoregressive
## error model

#read in results
list_model_results_3 = readRDS(file = "Autoerror_exp3_sim2_results_v6.rds")

B = 100

#coverage probability
beta_cancer_vec = rep(c(rep(-0.1, 4), rep(-0.35,4)),2)
cov_prob_ARE = vector(length = length(list_model_results_3))
for(i in 1:length(list_model_results_3)){
  cov_prob_ARE[i] = sum(apply(as.matrix(list_model_results_3[[i]]), MARGIN = 1, FUN = function(j) 
    j["Cancer effect"] - qnorm(0.975)*j["Cancer SE"] <= beta_cancer_vec[i] & 
      j["Cancer effect"] + qnorm(0.975)*j["Cancer SE"] >= beta_cancer_vec[i]))/B
}
names(cov_prob_ARE) = names(list_model_results_3)

saveRDS(cov_prob_ARE, file ="AREcov_prob_results.rds")

#average over N=50
mean(cov_prob[1:8]) #= 0.96
##average over N = 1000
mean(cov_prob[9:16]) #= 0.95125
#average over effect size of -0.1
mean(cov_prob[c(1:4,9:12)]) #= 0.95375
#average over effect size of -0.35
mean(cov_prob[c(5:8, 13:16)]) #= 0.9575
#average over sd_random = 1
mean(cov_prob[c(1,3,5,7,9,11,13,15)]) #= 0.95875
#average over sd_random=2.5
mean(cov_prob[c(2,4,6,8,10,12,14,16)]) #=0.9525
#average over sd_error=1
mean(cov_prob[c(1,2,5,6,9,10,13,14)]) #= 0.95625
#average over sd_error = 2.5
mean(cov_prob[c(3,4,7,8,11,12,15,16)]) #=0.955


#bias
avgbetaARE = sapply(list_model_results_3, function(i) mean(i[,"Cancer effect"]))
#names(avgbeta) = names(model_resultsLME)

biasARE = avgbetaARE - beta_cancer_vec
saveRDS(biasARE, file = "Bias_ARE.rds")

#sample sd of bias
sd(list_model_results_3[[1]][,"Cancer effect"] - beta_cancer_vec[1])
sdbiasARE = sapply(1:16, function(i) sd(list_model_results_3[[i]][,"Cancer effect"] - beta_cancer_vec[i]))
saveRDS(sdbiasARE, file = "SD_Bias_ARE.rds")

#sample sd of betahat
sdbetahatARE = sapply(1:16, function(i) sd(list_model_results_3[[i]][,"Cancer effect"]))
avgsebetahatARE = sapply(1:16, function(i) mean(list_model_results_3[[i]][,"Cancer SE"]))
saveRDS(sdbetahatARE, file = "SD_betahat_ARE.rds")
saveRDS(avgsebetahatARE, file = "Avgse_betahat_ARE.rds")


#CI for coverage prob
successARE = cov_prob_ARE*100
CI_ARE_L = binom.confint(x = successARE, n=100, methods = "ac")$lower
CI_ARE_U = binom.confint(x = successARE, n=100, methods = "ac")$upper

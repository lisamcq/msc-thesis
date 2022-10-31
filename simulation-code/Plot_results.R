#Code for visualizing bias and coverage probability of estimates from all 3 models

library(ggplot2)
library(dplyr)
library(tidyr)

#bias
biasARRmarg = readRDS(file = "Bias_marginalARR.rds")
biasARRcond = readRDS(file = "Bias_conditionalARR.rds")

biasLME = readRDS(file = "LMEbias_results.rds")
biasARE = readRDS(file = "Bias_ARE.rds")
sdbiasARR_marg = readRDS(file = "SD_Bias_ARR_marg.rds")
sdbiasARR_cond = readRDS(file = "SD_Bias_ARR_cond.rds")
sdbiasARE = readRDS(file = "SD_Bias_ARE.rds")
sdbiasLME = readRDS(file = "SD_Bias_LME.rds")

bias_df = data.frame(Run = 1:16, ARR_marg = biasARRmarg, ARR_cond = biasARRcond, ARE = biasARE, LME = biasLME)
bias_sd_df = data.frame(Run = 1:16, ARR_marg = sdbiasARR_marg, ARR_cond = sdbiasARR_cond, ARE = sdbiasARE, LME = sdbiasLME)
bias_df_long = pivot_longer(data = bias_df, cols = ARR_marg:LME, names_to = "Method", values_to = "Bias")
bias_sd_df_long = pivot_longer(data = bias_sd_df, cols = ARR_marg:LME, names_to = "Method", values_to = "Bias_SD")
bias_plot_data = merge(bias_df_long, bias_sd_df_long, all = TRUE)

bias_plot = ggplot(data = bias_plot_data, aes(x = Run, y= Bias, group = Method)) +
  geom_point(aes(colour = Method), size =5) + 
  geom_hline(yintercept = 0, size = 2) + 
  labs(y = "Estimated bias", x = "Run") +
  scale_x_continuous(minor_breaks = 1:16) +
  scale_colour_discrete(breaks = c("ARR_marg", "ARR_cond", "LME", "ARE"), 
                        labels = c("ARR Marginal", "ARR Conditional", "LME", "ARE"),
                        name = "Estimator") +
  geom_errorbar(aes(ymin = Bias-qnorm(0.975)*Bias_SD, ymax = Bias+qnorm(0.975)*Bias_SD, width=0.5, colour = Method))
bias_plot

#SE plot
sdbetahatARR_marg = readRDS(file = "SD_betahat_ARR_marg.rds")
sdbetahatARR_cond = readRDS(file = "SD_betahat_ARR_cond.rds")
sdbetahatLME = readRDS(file = "SD_betahat_LME.rds")
sdbetahatARE = readRDS(file = "SD_betahat_ARE.rds")
avgsebetahatARR_marg = readRDS(file = "Avgse_betahat_ARR_marg.rds")
avgsebetahatARR_cond = readRDS(file = "Avgse_betahat_ARR_cond.rds")
avgsebetahatARE = readRDS(file = "Avgse_betahat_ARE.rds")
avgsebetahatLME = readRDS(file = "Avgse_betahat_LME.rds")


SE_df = data.frame(Run = 1:16, ESDARR_marg = sdbetahatARR_marg, ESDARR_cond = sdbetahatARR_cond,
                   ESDARE = sdbetahatARE, ESDLME = sdbetahatLME, 
                   AvgSEARR_marg = avgsebetahatARR_marg, AvgSEARR_cond = avgsebetahatARR_cond,
                   AvgSEARE = avgsebetahatARE, AvgSELME = avgsebetahatLME)
ESE_df_long = pivot_longer(data = SE_df[,1:5], cols = ESDARR_marg:ESDLME, names_to = "Method", values_to = "Empirical_SD",
                           names_prefix = "ESD")
AvgSE_df_long = pivot_longer(data = SE_df[,c(1,6:9)], cols = AvgSEARR_marg:AvgSELME, names_to = "Method", values_to = "Average_SE",
                             names_prefix = "AvgSE")
SE_plot_data = merge(ESE_df_long, AvgSE_df_long, all = TRUE)

SE_plot = ggplot(data = SE_plot_data, aes(x = Empirical_SD, y = Average_SE, group = Method)) +
  geom_point(aes(colour=Method), size = 5)+
  geom_abline(slope = 1, intercept = 0, size = 2)+
  labs(x = "Empirical standard error", y ="Average model-based standard error") +
  scale_colour_discrete(breaks = c("ARR_marg", "ARR_cond", "LME", "ARE"), 
                        labels = c("ARR Marginal", "ARR Conditional", "LME", "ARE"),
                        name = "Estimator")
SE_plot

#Coverage probability plot
covprob_df = data.frame(Run = 1:16, ARR_marg = cov_prob_marg1, ARR_cond =cov_prob_cond, ARE = cov_prob_ARE, 
                        LME = cov_prob_LME)
covprob_CI_df = data.frame(Run = 1:16, L_ARR_marg = CI_ARR_marg_L, U_ARR_marg = CI_ARR_marg_U,
                           L_ARR_cond = CI_ARR_cond_L, U_ARR_cond = CI_ARR_cond_U,
                           L_ARE = CI_ARE_L, U_ARE = CI_ARE_U, L_LME = CI_LME_L, U_LME = CI_LME_U
                           )
covprob_df_long = pivot_longer(data = covprob_df, cols = ARR_marg:LME, names_to = "Method", values_to = "Covprob")
CIL_df_long = pivot_longer(data = covprob_CI_df[,c(1,2,4,6,8)], cols = L_ARR_marg:L_LME, names_to = "Method", values_to = "CI_L",
                           names_prefix = "L_")
CIU_df_long = pivot_longer(data = covprob_CI_df[,c(1,3,5,7,9)], cols = U_ARR_marg:U_LME, names_to = "Method", values_to = "CI_U",
                           names_prefix = "U_")
covprob_plot_data = merge(covprob_df_long, CIL_df_long, all = TRUE)
covprob_plot_data = merge(covprob_plot_data, CIU_df_long, all=TRUE)

covprob_plot = ggplot(data = covprob_plot_data, aes(x = Run, y= Covprob, group = Method)) +
  geom_point(aes(colour = Method), size =5) + 
  geom_hline(yintercept = 0.95, size = 2) + 
  scale_x_continuous(minor_breaks = 1:16) +
  scale_colour_discrete(breaks = c("ARR_marg", "ARR_cond", "LME", "ARE"), 
                        labels = c("ARR Marginal", "ARR Conditional", "LME", "ARE"),
                        name = "Estimator") +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_U, width=0.5, colour = Method)) +
  labs(y = "Estimated coverage probability")
covprob_plot



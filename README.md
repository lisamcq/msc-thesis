# msc-thesis
My master's thesis and code from the simulation study.

My master's thesis is posted in this repository and available for reading if you think you would enjoy 63 pages of in-depth discussion on statistical models and their mean/variance-covariance structures. I have included the abstract here for you to gauge your interest.

Title: Autoregressive mixed effects models and an application to annual income of cancer survivors
Abstract: Longitudinal observations of income are often strongly autocorrelated, even after adjusting for independent variables. We explore two common longitudinal models that allow for residual autocorrelation: 1. the autoregressive error model (a linear mixed effects model with an AR(1) covariance structure), and 2. the autoregressive response model (a linear mixed effects model that includes the first lag of the response variable as an independent variable). We explore the theoretical properties of these models and illustrate the behaviour of parameter estimates using a simulation study. Additionally, we apply the models to a data set containing repeated (annual) observations of income and sociodemographic variables on a sample of breast cancer survivors. Our preliminary results suggest that the autoregressive
response model may severely overestimate the magnitude of the effect of cancer. Our findings will guide future, comprehensive study of the short- and long-term effects of a breast cancer diagnosis on a survivor’s annual net income.

The repository also contains the code for the simulation study. The simulation study investigated the performance of the three models discussed (linear mixed effects model with a random intercept and independent errors, autoregressive error model, and autoregressive response model) when applied to data simulated from the autoregressive response model. To understand the simulation study and code, I suggest reading Chapter 3 (Methods) and Chapter 5 (Simulation Study) of my thesis at minimum.

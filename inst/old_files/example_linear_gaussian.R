#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file provides an example of univariate Linear Gaussian state-space model
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# set working directory to directory with all the required R files
setwd("C:/Users/jianfushao/Documents/Harvard/MyHarvard/Research/Model Selection for SSM (Pierre and Jie)")

# load the different functions
source("functions_model.R")
source("functions_hscore_continuous.R")

# load the different models (NOTE: the priors are specified inside the model file)
source("model_linear_gaussian_univariate.R")

#------------ Simulate some data from a Linear Gaussian SSM ----------#
T = 50
sim = simulateData(model.lineargaussian,theta = c(0.8,1,1,1),T)
X = sim$X
Y = sim$Y

#---------- Compute the prequential Hyvarinen score for the given model ----------#
Nx = 100
Ntheta = 100
hscore = hyvarinenContinuous(Y,T,model.lineargaussian,Ntheta,Nx)
plot(hscore$Hscore,type = 'l',xlab = "Time horizon T",ylab = "SH(T)")



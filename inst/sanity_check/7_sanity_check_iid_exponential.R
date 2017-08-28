##################################################################################################
# This implements example 3.3. in O'Hagan (1995).
# Model 1 = iid expo(theta) with reference prior
# Model 2 = iid Log-Normal(mu,sigma2) with reference prior
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(1)
#--------------------------------------------------------------------------------------------
# Generate some data
nobservations = 30
Y = rexp(nobservations,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$ess_threshold = 0.5
#--------------------------------------------------------------------------------------------
a0 = 1
b0 = 1
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
results = data.frame()
post_all = data.frame()
result = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  hscore(observations, get_model_iid_exponential(a0,b0), algorithmic_parameters)
}
for (r in 1:repl){
  results = rbind(results,data.frame(logevidence = result[[r]]$logevidence,
                                     hscore = result[[r]]$hscore,
                                     time = 1:nobservations,
                                     model = 1,
                                     repl = r))
}
#--------------------------------------------------------------------------------------------
a_post = rep(NA,nobservations)
b_post = rep(NA,nobservations)
for (t in 1:nobservations){
  a_post[t] = a0 + t
  b_post[t] = b0 + sum(observations[,1:t])
}
M = 2000
post_exact = data.frame()
for (r in 1:repl){
  for (t in 1:nobservations){
    post_exact = rbind(post_exact,data.frame(lambda = rgamma(M,shape = a_post[t],rate = b_post[t]),
                                             W = 1/M,
                                             time = t,
                                             model = 1,
                                             repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# Compute the "exact" h-score for model 1 (exponential)
results$hscore_analytical = NA
for (r in 1:repl){
  for (t in 1:nobservations){
    yt = observations[,t]
    post = subset(post_exact,time==t & repl==r)
    theta1 = -post$lambda
    post_mean = rbind(sum(theta1*post$W))
    post_var = cov.wt(cbind(theta1),wt = post$W)$cov
    results$hscore_analytical[results$time==t&results$repl==r] = (post_mean)^2 + 2*post_var
  }
  results$hscore_analytical[results$repl==r] = cumsum(results$hscore_analytical[results$repl==r])
}
# Check the h-score
ggplot(results) +
  geom_line(aes(time, hscore/time, color = factor(model), group=interaction(model,repl))) +
  geom_point(aes(time, hscore_analytical/time, group=interaction(model,repl)), col ="blue")

#--------------------------------------------------------------------------------------------



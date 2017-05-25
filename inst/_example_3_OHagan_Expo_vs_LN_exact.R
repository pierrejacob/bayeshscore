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
algorithmic_parameters$ess_threshold = 0.5
#--------------------------------------------------------------------------------------------
# hyperparameters for model 1
a0 = 1
b0 = 1
# hyperparameters for model 2
mu0 = 0
kappa0 = 1
nu0 = 1
sigma02 = 1
#--------------------------------------------------------------------------------------------
# define models
# Note: set a, b, nu0 and kappa0 close to 0 close to infinity for improper reference priors
model = function(i){
  if (i==1) {return(get_model_iid_exponential(a0, b0))}
  if (i==2) {return(get_model_iid_lognormal(mu0, kappa0, nu0, sigma02))}
}
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
results = data.frame()
for (m in 1:2){
  result = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results = rbind(results,data.frame(logevidence = result[[r]]$logevidence,
                                       hscore = result[[r]]$hscore,
                                       time = 1:nobservations,
                                       model = m,
                                       repl = r))
  }
}
#--------------------------------------------------------------------------------------------
M = 2^10
# Samples exact posterior from model 1
post_exact1 = data.frame()
a_post = rep(NA,nobservations)
b_post = rep(NA,nobservations)
for (t in 1:nobservations){
  a_post[t] = a0 + t
  b_post[t] = b0 + sum(observations[,1:t])
}
for (r in 1:repl){
  for (t in 1:nobservations){
    post_exact1 = rbind(post_exact1,data.frame(lambda = rgamma(M,shape = a_post[t],rate = b_post[t]),
                                               W = 1/M,
                                               time = t,
                                               model = 1,
                                               repl = r))
  }
}
# Samples exact posterior from model 2
post_exact2 = data.frame()
mu_post = rep(NA,nobservations)
kappa_post = rep(NA,nobservations)
nu_post = rep(NA,nobservations)
sigma2_post = rep(NA,nobservations)
for (t in 1:nobservations){
  logy_1_t = log(observations[,1:t])
  mt = mean(logy_1_t)
  st2 = var(logy_1_t)
  mu_post[t] = mu0*kappa0/(kappa0+t) + mt*t/(kappa0+t)
  kappa_post[t] = kappa0 + t
  nu_post[t] = nu0 + t
  sigma2_post[t] = (1/nu_post[t]) * (nu0*sigma02 + sum((logy_1_t-mt)^2) + (kappa0*t/(kappa0+t))*(mt-mu0)^2)
}
for (r in 1:repl){
  for (t in 1:nobservations){
    exact_draw = rnorminvchisq(M, mu_post[t], kappa_post[t], nu_post[t], sigma2_post[t])
    post_exact2 = rbind(post_exact2,data.frame(mu = exact_draw[1,],
                                               sigma2 = exact_draw[2,],
                                               W = 1/M,
                                               time = t,
                                               model = 2,
                                               repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# Compute the "exact" h-score for model 1 (exponential)
results$hscore_analytical = NA
for (r in 1:repl){
  for (t in 1:nobservations){
    yt = observations[,t]
    post = subset(post_exact1,time==t & repl==r)
    theta1 = -post$lambda
    post_mean = rbind(sum(theta1*post$W))
    post_var = cov.wt(cbind(theta1),wt = post$W)$cov
    results$hscore_analytical[results$time==t&results$repl==r&results$model==1] = (post_mean)^2 + 2*post_var
  }
  results$hscore_analytical[results$repl==r&results$model==1] = cumsum(results$hscore_analytical[results$repl==r&results$model==1])
}
# Compute the "exact" h-score for model 2 (log-normal)
for (r in 1:repl){
  for (t in 1:nobservations){
    yt = observations[,t]
    post = subset(post_exact2,time==t & repl==r)
    theta1 = post$mu/post$sigma2
    theta2 = -1/(2*post$sigma2)
    post_mean = rbind(sum(theta1*post$W),sum(theta2*post$W))
    post_var = cov.wt(cbind(theta1,theta2),wt = post$W)$cov
    d1a = -1/yt
    d2a = 1/(yt^2)
    d = rbind(-1/(yt^2),-(2/(yt^2))*(log(yt)-1))
    Jt = cbind(1/yt, 2*log(yt)/yt)
    results$hscore_analytical[results$time==t&results$repl==r&results$model==2] = 2*d2a + 2*t(d)%*%post_mean + (d1a + Jt%*%post_mean)^2 + 2*Jt%*%post_var%*%t(Jt)
  }
  results$hscore_analytical[results$repl==r&results$model==2] = cumsum(results$hscore_analytical[results$repl==r&results$model==2])
}
# Check the h-score
ggplot(results) +
  geom_line(aes(time, hscore/time, color = factor(model), group=interaction(model,repl))) +
  geom_point(aes(time, hscore_analytical/time, color = factor(model), group=interaction(model,repl)))

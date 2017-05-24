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
nobservations = 50
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
mu0 = 0
kappa0 = 1
nu0 = 1
sigma02 = 1
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
results = data.frame()
post_all = data.frame()
result = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  hscore(observations, get_model_iid_lognormal(mu0, kappa0,nu0,sigma02), algorithmic_parameters)
}
for (r in 1:repl){
  results = rbind(results,data.frame(logevidence = result[[r]]$logevidence,
                                     hscore = result[[r]]$hscore,
                                     time = 1:nobservations,
                                     model = 2,
                                     repl = r))
  for (t in 1:nobservations){
    post_all = rbind(post_all,data.frame(mu = c(result[[r]]$thetas_history[[t+1]][1,]),
                                         sigma2 = c(result[[r]]$thetas_history[[t+1]][2,]),
                                         W = result[[r]]$normw_history[[t+1]],
                                         time = t,
                                         model = 2,
                                         repl = r))
  }
}
#--------------------------------------------------------------------------------------------
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
M = 5000
post_exact = data.frame()
for (r in 1:repl){
  for (t in 1:nobservations){
    exact_draw = rnorminvchisq(M, mu_post[t], kappa_post[t], nu_post[t], sigma2_post[t])
    post_exact = rbind(post_exact,data.frame(mu = exact_draw[1,],
                                             sigma2 = exact_draw[2,],
                                             W = 1/M,
                                             time = t,
                                             model = 2,
                                             repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# Compute the "exact" h-score for model 2 (log-normal)
results$hscore_analytical = NA
for (r in 1:repl){
  for (t in 1:nobservations){
    yt = observations[,t]
    post = subset(post_exact,time==t & repl==r)
    theta1 = post$mu/post$sigma2
    theta2 = -1/(2*post$sigma2)
    post_mean = rbind(sum(theta1*post$W),sum(theta2*post$W))
    post_var = cov.wt(cbind(theta1,theta2),wt = post$W)$cov
    d1a = -1/yt
    d2a = 1/(yt^2)
    d = rbind(-1/(yt^2),-(2/(yt^2))*(log(yt)-1))
    Jt = cbind(1/yt, 2*log(yt)/yt)
    results$hscore_analytical[results$time==t&results$repl==r] = 2*d2a + 2*t(d)%*%post_mean + (d1a + Jt%*%post_mean)^2 + 2*Jt%*%post_var%*%t(Jt)
  }
  results$hscore_analytical[results$repl==r] = cumsum(results$hscore_analytical[results$repl==r])
}
# Check the h-score
ggplot(results) +
  geom_line(aes(time, hscore/time, color = factor(model), group=interaction(model,repl))) +
  geom_point(aes(time, hscore_analytical/time, group=interaction(model,repl)), col ="blue")

#--------------------------------------------------------------------------------------------



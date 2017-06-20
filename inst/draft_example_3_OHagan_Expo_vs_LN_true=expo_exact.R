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
nobservations = 50000
Y = rexp(nobservations,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
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
# repl = 5 #number of replications
# registerDoParallel(cores=5) #number of workers in parallel
# #--------------------------------------------------------------------------------------------
# # Monitor progress in parallel via log file
# logfilename = "results.log"
# writeLines(c(""), logfilename)
# sink(logfilename, append = TRUE)
# #--------------------------------------------------------------------------------------------
# results = data.frame()
# for (m in 1:2){
#   gc() # attempt to limit RAM usage
#   cat("Model ",toString(m)," started at:", toString(Sys.time()))
#   result = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
#     gc() # attempt to limit RAM usage
#     sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
#     hscore(observations, model(m), algorithmic_parameters)
#   }
#   for (r in 1:repl){
#     results = rbind(results,data.frame(logevidence = result[[r]]$logevidence,
#                                        hscore = result[[r]]$hscore,
#                                        time = 1:nobservations,
#                                        model = m,
#                                        repl = r))
#   }
# }
# sink()
#--------------------------------------------------------------------------------------------
M = 2^10
# Compute the "exact" h-score for model 1 (exponential)
hscore_analytical1 = rep(NA,nobservations)
a_post = rep(NA,nobservations)
b_post = rep(NA,nobservations)
progressbar = txtProgressBar(0,nobservations,style = 3)
for (t in 1:nobservations){
  a_post[t] = a0 + t
  b_post[t] = b0 + sum(observations[,1:t])
  theta1_post = rgamma(M,shape = a_post[t],rate = b_post[t])
  post_mean = mean(theta1_post)
  post_var = var(theta1_post)
  hscore_analytical1[t] = (post_mean)^2 + 2*post_var
  setTxtProgressBar(progressbar,t)
}
hscore_analytical1 = cumsum(hscore_analytical1)

# Compute the "exact" h-score for model 2 (log-normal)
hscore_analytical2 = rep(NA,nobservations)
post_exact2 = data.frame()
mu_post = rep(NA,nobservations)
kappa_post = rep(NA,nobservations)
nu_post = rep(NA,nobservations)
sigma2_post = rep(NA,nobservations)
progressbar = txtProgressBar(0,nobservations,style = 3)
for (t in 1:nobservations){
  logy_1_t = log(observations[,1:t])
  mt = mean(logy_1_t)
  st2 = var(logy_1_t)
  mu_post[t] = mu0*kappa0/(kappa0+t) + mt*t/(kappa0+t)
  kappa_post[t] = kappa0 + t
  nu_post[t] = nu0 + t
  sigma2_post[t] = (1/nu_post[t]) * (nu0*sigma02 + sum((logy_1_t-mt)^2) + (kappa0*t/(kappa0+t))*(mt-mu0)^2)
  draw_post = rnorminvchisq(M, mu_post[t], kappa_post[t], nu_post[t], sigma2_post[t])
  mu_post = draw_post[1,]
  sigma2_post = draw_post[2,]
  theta1_post = mu_post/sigma2_post
  theta2_post = -1/(2*sigma2_post)
  post_mean = rbind(mean(theta1_post),mean(theta2_post))
  post_var = var(cbind(theta1_post,theta2_post))
  yt = observations[,t]
  d1a = -1/yt
  d2a = 1/(yt^2)
  d = rbind(-1/(yt^2),-(2/(yt^2))*(log(yt)-1))
  Jt = cbind(1/yt, 2*log(yt)/yt)
  hscore_analytical2[t] = 2*d2a + 2*t(d)%*%post_mean + (d1a + Jt%*%post_mean)^2 + 2*Jt%*%post_var%*%t(Jt)
  setTxtProgressBar(progressbar,t)
}
hscore_analytical2 = cumsum(hscore_analytical2)

# Check the h-score
ggplot() +
  # geom_line(data=results,aes(time, hscore/time, color = factor(model), group=interaction(model,repl))) +
  # geom_point(aes(1:nobservations, hscore_analytical1/(1:nobservations)),col="red") +
  # geom_point(aes(1:nobservations, hscore_analytical2/(1:nobservations)),col="blue") +
  geom_line(aes(1:nobservations, hscore_analytical1/(1:nobservations)),col="red") +
  geom_line(aes(1:nobservations, hscore_analytical2/(1:nobservations)),col="blue")



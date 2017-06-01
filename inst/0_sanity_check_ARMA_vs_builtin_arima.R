##################################################################################################
# Sanity check ARMA vs arima built-in package
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
# library(doParallel)
# library(foreach)
set.seed(1)

nu0 = 1
sigma02 = 1
p = 2
q = 1
# # set algorithmic parameters
# algorithmic_parameters = list()
# algorithmic_parameters$Ntheta = 2^10
# algorithmic_parameters$verbose = TRUE
# algorithmic_parameters$store_theta = TRUE
# algorithmic_parameters$store_byproducts= TRUE
# #--------------------------------------------------------------------------------------------
# repl = 5 #number of replications
# registerDoParallel(cores=5) #number of workers in parallel
# #--------------------------------------------------------------------------------------------
nobservations = 20
model = get_model_ARMA(p,q,nu0,sigma02)

get_stationaryparameters = function(p,q){
  AR_coeffs = runif(p,-1,1)
  MA_coeffs = runif(q,-1,1)
  AR_roots = polyroot(c(1,-AR_coeffs))
  MA_roots = polyroot(c(1,-MA_coeffs))
  accept = ((sum(abs(AR_roots)<=1))==0)&&((sum(abs(MA_roots)<=1))==0)&&(length(intersect(AR_roots,MA_roots))==0)
  while (!accept){
    AR_coeffs = runif(p,-1,1)
    MA_coeffs = runif(q,-1,1)
    accept = ((sum(abs(polyroot(c(1,-AR_coeffs)))<=1))==0)&&((sum(abs(polyroot(c(1,-MA_coeffs)))<=1))==0)
  }
  return (c(AR_coeffs,MA_coeffs,1))
}
theta = get_stationaryparameters(p,q)
observations = simulateData(model,theta,nobservations)$Y


# results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
#   hscore(observations, get_model_ARMA(p,q,nu0,sigma02), algorithmic_parameters)
# }

ll = function(theta,t){
  KF = model$initialize_byproducts(theta,observations[,1:t,drop=FALSE])
  for (i in 1:t){
    KF = model$update_byproduct(KF,i,theta,observations[,1:t,drop=FALSE])
  }
  return (model$likelihood(observations[,1:t,drop=FALSE],t,theta,KF))
}
# approxMLE = function(r,t){
#   result = results[[r]]
#   iMLE = which.max(result$normw_history[[t+1]])
#   approxMLE = result$thetas_history[[t+1]][,iMLE]
#   KF = result$byproducts_history[[t+1]][[iMLE]]
#   return (optim(approxMLE,function(theta)ll(theta,t), control = list(fnscale = -1))$value)
# }
llMLE = function(thetainit,t){
  # cat("\n Step",t,thetainit,"\n")
  return (optim(thetainit,function(theta)ll(theta,t), control = list(fnscale = -1))$value)
}

nmin = 10
# approxMLE.df = data.frame()
# for (r in 1:repl){
#   approxMLE.df = rbind(approxMLE.df, data.frame(time = nmin:nobservations,
#                                                 r = r,
#                                                 MaxLL = sapply(nmin:nobservations,function(t)approxMLE(r,t))))
# }
approxLLMLE = rep(NA,nobservations-nmin+1)
arimaLLMLE = rep(NA,nobservations-nmin+1)
pb = txtProgressBar(0,length(approxLLMLE),style = 3)
for (t in nmin:nobservations) {
  arimafit = arima(c(observations[,1:t]),c(p,0,q),include.mean = FALSE)
  i = t-nmin+1
  arimaLLMLE[i] = arimafit$loglik
  approxLLMLE[i] = llMLE(c(arimafit$coef,arimafit$sigma2),t)
  setTxtProgressBar(pb,i)
}
# arimaLLMLE = sapply(nmin:nobservations,function(t)arima(c(observations[,1:t]),c(p,0,q),include.mean = FALSE)$loglik)
# ggplot() +
#   geom_line(data = approxMLE.df, aes(time,MaxLL,group = r)) +
#   geom_line(aes(nmin:nobservations,arimaLLMLE),col="blue",size=2)
ggplot() +
  geom_line(aes(nmin:nobservations,approxLLMLE), size = 1) +
  geom_line(aes(nmin:nobservations,arimaLLMLE),col="red",linetype="dashed",size=2)

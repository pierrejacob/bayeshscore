##################################################################################################
# TOy example: Kernel density estimation for log-derivatives of the predictive distributions
##################################################################################################
rm(list = ls())
library(ks)
library(HyvarinenSSM)
library(numDeriv)
library(ggplot2)
library(gridExtra)
set.seed(29)

# arbitrary observation at which the derivatives are to be estimated
mu = 0
sigma2 = 1

Yt = rnorm(1,mu,sqrt(sigma2))
sigmastar = 2
sigma2star = (sigmastar)^2

# # Define hyperparameters
Ny = 500
####################################### NORMAL ######################################
Y = rnorm(Ny)
# lower <- Yt - 0.5
# upper <- Yt + 0.5
lower <- -4
upper <- 4
# Y = rep(NA,Ny)
# for (i in 1:Ny){
#   draw = rnorm(1)
#   while ((draw < lower)||(draw > upper)){
#     draw = rnorm(1)
#   }
#   Y[i] = draw
# }


####################################### CAUCHY ######################################
# epsilon = 2
# upper <- Yt + epsilon
# lower <- Yt - epsilon
# Y = rep(NA,Ny)
# for (i in 1:Ny){
#   draw = rcauchy(1)
#   while ((draw < lower)||(draw > upper)){
#     draw = rcauchy(1)
#   }
#   Y[i] = draw
# }

####################################### NORMAL ######################################
##### compute density estimators
xgrid = seq(lower,upper,0.001)
# exact log-density
ll_exact = dnorm(xgrid,log=TRUE)
l_exact = exp(ll_exact)
# via KDE
l_kde = kdde(Y,deriv.order = 0,eval.points = xgrid)$estimate
ll_kde = log(l_kde)
# # first derivative of density
d1_kde = kdde(Y,deriv.order = 1,eval.points = xgrid)$estimate
d1_exact = grad(function(x)dnorm(x),xgrid)
# first derivative of log-density
d1ll_kde = d1_kde/l_kde
d1ll_exact = grad(func = function(x) dnorm(x, log = TRUE), xgrid)
# # second derivative of density
d2_kde = kdde(Y,deriv.order = 2,eval.points = xgrid)$estimate
d2_exact = sapply(xgrid, function(y) hessian(function(x) dnorm(x), y))
# second derivative of log-density
d2ll_kde = d2_kde/l_kde - (d1_kde/l_kde)^2
d2ll_exact = sapply(xgrid, function(y) hessian(func = function(x) dnorm(x, log = TRUE), y))


# ####################################### CAUCHY ######################################
# ##### compute density estimators
# xgrid = seq(lower,upper,0.001)
# # exact log-density
# ll_exact = dcauchy(xgrid,log=TRUE) - log(pcauchy(upper) - pcauchy(lower))
# l_exact = exp(ll_exact)
# # via KDE
# l_kde = kdde(Y,deriv.order = 0,eval.points = xgrid)$estimate
# ll_kde = log(l_kde)
# # # first derivative of density
# d1_kde = kdde(Y,deriv.order = 1,eval.points = xgrid)$estimate
# d1_exact = grad(function(x)dcauchy(x),xgrid)
# # first derivative of log-density
# d1ll_kde = d1_kde/l_kde
# d1ll_exact = grad(func = function(x) dcauchy(x, log = TRUE), xgrid)
# # # second derivative of density
# d2_kde = kdde(Y,deriv.order = 2,eval.points = xgrid)$estimate
# d2_exact = sapply(xgrid, function(y) hessian(function(x) dcauchy(x), y))
# # second derivative of log-density
# d2ll_kde = d2_kde/l_kde - (d1_kde/l_kde)^2
# d2ll_exact = sapply(xgrid, function(y) hessian(func = function(x) dcauchy(x, log = TRUE), y))



# density + error log-density
g0 = ggplot() +
  geom_histogram(data = data.frame(Y = Y), aes(Y,..density..),alpha=0.6) +
  geom_line(aes(xgrid,l_kde),col="blue") +
  geom_line(aes(xgrid,l_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)
g0bis = ggplot() +
  geom_line(aes(xgrid,(ll_kde - ll_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed")

# first derivative + absolute error
g1l = ggplot() +
  geom_line(aes(xgrid,d1_kde),col="blue") +
  geom_line(aes(xgrid,d1_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)
g1lbis = ggplot() +
  geom_line(aes(xgrid,abs(d1_kde - d1_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()

# first log-derivative + absolute error
g1ll = ggplot() +
  geom_line(aes(xgrid,d1ll_kde),col="blue") +
  geom_line(aes(xgrid,d1ll_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)
g1llbis = ggplot() +
  geom_line(aes(xgrid,abs(d1ll_kde - d1ll_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()

# second derivative + absolute error
g2l = ggplot() +
  geom_line(aes(xgrid,d2_kde),col="blue") +
  geom_line(aes(xgrid,d2_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)
g2lbis = ggplot() +
  geom_line(aes(xgrid,abs(d2_kde - d2_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()

# second log-derivative + absolute error
g2ll = ggplot() +
  geom_line(aes(xgrid,d2ll_kde),col="blue") +
  geom_line(aes(xgrid,d2ll_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)
g2llbis = ggplot() +
  geom_line(aes(xgrid,abs(d2ll_kde - d2ll_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()


grid.arrange(g0,g0bis,g1l,g1lbis,g2l,g2lbis,nrow=3)
# grid.arrange(g0,g0bis,g1ll,g1llbis,g2ll,g2llbis,nrow=3)


## direct derivative estimation method
range_lambda = 10^seq(-3,3,0.25)
range_sigma2 = 10^seq(-4,1,0.1625)
kfold = 10

weighted = TRUE
Ysmaller = Y[sample(1:length(Y),min(length(Y),500),replace=FALSE)]

if(weighted){
  # cvfit1 = get_CV_RBF_weighted(Ysmaller, kfold, range_lambda, range_sigma2, order = 1, Yt, sigma2star)
  # sigma21 = cvfit1$sigma2
  # lambda1 = cvfit1$lambda
  # fit1 = get_derivative_RBF_weighted(Y,lambda1,sigma21,1,Yt,sigma2star,xgrid)
  # theta1 = fit1$theta
  cvfit1 = get_CV_RBF_weighted_local(Y, kfold, range_lambda, range_sigma2, order = 1, Yt, sigma2star)
  sigma21 = cvfit1$sigma2
  lambda1 = cvfit1$lambda
  fit1 = get_derivative_RBF_weighted_local(Y,lambda1,sigma21,1,Yt,sigma2star,xgrid)
  theta1 = fit1$theta
} else {
  cvfit1 = get_CV_RBF(Ysmaller, kfold, range_lambda, range_sigma2, order = 1)
  sigma21 = cvfit1$sigma2
  lambda1 = cvfit1$lambda
  fit1 = get_derivative_RBF(Y,lambda1,sigma21,1,xgrid)
  theta1 = fit1$theta
}

cat("sigma2 =",sigma21,"\n")
cat("lambda =",lambda1,"\n")

other_d1s = fit1$estimate
d1est = get_estimate_RBF(Yt,Y,theta1,sigma21)

ggplot() +
  geom_line(aes(xgrid,d1_kde),col="blue",size=1) +
  geom_line(aes(xgrid,d1_exact),col="red",size=2,linetype="dashed") +
  geom_line(aes(xgrid,other_d1s),col="orange",size=1) +
  geom_vline(xintercept = Yt)

#####################################################

if (weighted){
  # cvfit2 = get_CV_RBF_weighted(Ysmaller, kfold, range_lambda, range_sigma2, order = 2, Yt, sigma2star)
  # sigma22 = cvfit1$sigma2
  # lambda2 = cvfit1$lambda
  # fit2 = get_derivative_RBF_weighted(Y,lambda2,sigma22,2,Yt,sigma2star,xgrid)
  # theta2 = fit2$theta
  cvfit2 = get_CV_RBF_weighted_local(Y, kfold, range_lambda, range_sigma2, order = 2, Yt, sigma2star)
  sigma22 = cvfit2$sigma2
  lambda2 = cvfit2$lambda
  fit2 = get_derivative_RBF_weighted_local(Y,lambda2,sigma22,2,Yt,sigma2star,xgrid)
  theta2 = fit2$theta
} else {
  cvfit2 = get_CV_RBF(Ysmaller, kfold, range_lambda, range_sigma2, order = 2)
  sigma22 = cvfit2$sigma2
  lambda2 = cvfit2$lambda
  fit2 = get_derivative_RBF(Y,lambda2,sigma22,2,xgrid)
  theta2 = fit2$theta
}

cat("sigma2 =",sigma22,"\n")
cat("lambda =",lambda2,"\n")

other_d2s = fit2$estimate
d2est = get_estimate_RBF(Yt,Y,theta2,sigma22)

ggplot() +
  geom_line(aes(xgrid,d2_kde),col="blue",size=1) +
  geom_line(aes(xgrid,d2_exact),col="red",size=2,linetype="dashed") +
  geom_line(aes(xgrid,other_d2s),col="orange",size=1) +
  geom_vline(xintercept = Yt)


#####################################################
### WARNING: if conditional sampling, need to rescale the derivatives appropriately before comparing !!
cat("d1 exact:",grad(function(x)dnorm(x), Yt),"\n")
cat("d1 kde:",kdde(Y,deriv.order = 1,eval.points = Yt)$estimate,"\n")
cat("d1 est:",d1est,"\n")

cat("d2 exact:",hessian(function(x)dnorm(x), Yt),"\n")
cat("d2 kde:",kdde(Y,deriv.order = 2,eval.points = Yt)$estimate,"\n")
cat("d2 est:",d2est,"\n")

cat("d1/d0 exact:",grad(function(x)dnorm(x),Yt)/dnorm(Yt),"\n")
cat("d1/d0 est:",d1est/kdde(Y,deriv.order = 0,eval.points = Yt)$estimate,"\n")

cat("d2/d0 exact:",hessian(function(x)dnorm(x),Yt)/dnorm(Yt),"\n")
cat("d2/d0 est:",d2est/kdde(Y,deriv.order = 0,eval.points = Yt)$estimate,"\n")


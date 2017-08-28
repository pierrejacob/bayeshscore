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

# rtest = rnorm
# dtest = dnorm
# ptest = pnorm
rtest = rcauchy
dtest = dcauchy
ptest = pcauchy

# # Define hyperparameters
Ny = 5000
####################################### Draws ######################################
epsilon = 10
upper <- Yt + epsilon
lower <- Yt - epsilon
Y = rep(NA,Ny)
for (i in 1:Ny){
  draw = rtest(1)
  while ((draw < lower)||(draw > upper)){
    draw = rtest(1)
  }
  Y[i] = draw
}
xgrid = seq(lower,upper,0.01)
normalizing = ptest(upper) - ptest(lower)
# Y = rtest(Ny)
# xgrid = seq(-4,4,0.01)
# normalizing = 1

####################################### Estimation ######################################
##### compute density estimators
# exact log-density
ll_exact = dtest(xgrid,log=TRUE) - log(normalizing)
l_exact = exp(ll_exact)
# via KDE
l_kde = kdde(Y,deriv.order = 0,eval.points = xgrid)$estimate
ll_kde = log(l_kde)
# # first derivative of density
d1_kde = kdde(Y,deriv.order = 1,eval.points = xgrid)$estimate
d1_exact = grad(function(x)dtest(x),xgrid)
# first derivative of log-density
d1ll_kde = d1_kde/l_kde
d1ll_exact = grad(func = function(x) dtest(x, log = TRUE), xgrid)
# # second derivative of density
d2_kde = kdde(Y,deriv.order = 2,eval.points = xgrid)$estimate
d2_exact = sapply(xgrid, function(y) hessian(function(x) dtest(x), y))
# second derivative of log-density
d2ll_kde = d2_kde/l_kde - (d1_kde/l_kde)^2
d2ll_exact = sapply(xgrid, function(y) hessian(func = function(x) dtest(x, log = TRUE), y))


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
range_lambda = 10^seq(-1,1,0.25)
range_sigma2 = 10^seq(-0.3,1,0.1625)
kfold = 10

Ysmaller = Y[sample(1:length(Y),min(length(Y),500),replace=FALSE)]

cvfit1 = get_CV_RBF(Ysmaller, kfold, range_lambda, range_sigma2, order = 1)
sigma21 = cvfit1$sigma2
lambda1 = cvfit1$lambda

fit1 = get_derivative_RBF(Y,lambda1,sigma21,1,xgrid)
theta1 = fit1$theta

other_d1s = fit1$estimate

ggplot() +
  geom_line(aes(xgrid,d1_kde),col="blue",size=1) +
  geom_line(aes(xgrid,d1_exact),col="red",size=2,linetype="dashed") +
  geom_line(aes(xgrid,other_d1s),col="orange",size=1) +
  geom_vline(xintercept = Yt)

#####################################################

cvfit2 = get_CV_RBF(Ysmaller, kfold, range_lambda, range_sigma2, order = 2)
sigma22 = cvfit2$sigma2
lambda2 = cvfit2$lambda


fit2 = get_derivative_RBF(Y,lambda2,sigma22,2,xgrid)
theta2 = fit2$theta

other_d2s = fit2$estimate

ggplot() +
  geom_line(aes(xgrid,d2_kde),col="blue",size=1) +
  geom_line(aes(xgrid,d2_exact),col="red",size=2,linetype="dashed") +
  geom_line(aes(xgrid,other_d2s),col="orange",size=1) +
  geom_vline(xintercept = Yt)



### WARNING: if conditional sampling, need to rescale the derivatives appropriately before comparing !!
cat("d1 exact:",grad(function(x)dtest(x), Yt),"\n")
cat("d1 kde:",kdde(Y,deriv.order = 1,eval.points = Yt)$estimate,"\n")
cat("d1 est:",get_estimate_RBF(Yt,Y,theta1,sigma21),"\n")

cat("d2 exact:",hessian(function(x)dtest(x), Yt),"\n")
cat("d2 kde:",kdde(Y,deriv.order = 2,eval.points = Yt)$estimate,"\n")
cat("d2 est:",get_estimate_RBF(Yt,Y,theta2,sigma22),"\n")

cat("d1/d0 exact:",grad(function(x)dtest(x),Yt)/dtest(Yt),"\n")
cat("d1/d0 est:",get_estimate_RBF(Yt,Y,theta1,sigma21)/kdde(Y,deriv.order = 0,eval.points = Yt)$estimate,"\n")

cat("d2/d0 exact:",hessian(function(x)dtest(x),Yt)/dtest(Yt),"\n")
cat("d2/d0 est:",get_estimate_RBF(Yt,Y,theta2,sigma22)/kdde(Y,deriv.order = 0,eval.points = Yt)$estimate,"\n")



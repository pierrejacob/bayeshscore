##################################################################################################
# TOy example: Kernel density estimation for log-derivatives of the predictive distributions
##################################################################################################
rm(list = ls())
library(kedd)
library(numDeriv)
library(ggplot2)
library(gridExtra)
set.seed(29)

# arbitrary observation at which the derivatives are to be estimated
mu = 0
sigma2 = 1

Yt = rnorm(1,mu,sqrt(sigma2))

# Define hyperparameters
Ny = 10^3
# tol = 3
# xgrid = seq(Yt-tol,Yt+tol,length.out = 100)
# Y = rep(NA,Ny)
# for (i in 1:Ny){
#   draw = rnorm(1,mu,sqrt(sigma2))
#   while ((draw<Yt-tol)||(draw>Yt+tol)){
#     draw = rnorm(1,mu,sqrt(sigma2))
#   }
#   Y[i] = draw
# }
Y = rnorm(Ny,mu,sqrt(sigma2))
xgrid = seq(-3,3,length.out = 100)

##### compute density estimators
kernel = "tricube"
# log-density
l_exact = dnorm(xgrid,mu,sqrt(sigma2),log=FALSE)
# l_exact = dnorm(xgrid,mu,sqrt(sigma2),log=FALSE)/(pnorm(Yt+tol)-pnorm(Yt-tol))
l_kde = dkde(x = Y,deriv.order = 0,y = xgrid, kernel = kernel)$est.fx
ll_exact = dnorm(xgrid,mu,sqrt(sigma2),log=TRUE)
# ll_exact = dnorm(xgrid,mu,sqrt(sigma2),log=TRUE) - log((pnorm(Yt+tol)-pnorm(Yt-tol)))
ll_kde = log(l_kde)
# # first derivative of density
# d1_kde = kdde(Y,deriv.order = 1,eval.points = xgrid)$estimate
# d1_exact = grad(function(x)dnorm(x,mu,sqrt(sigma2)),xgrid)
# first derivative of log-density
d1_kde = dkde(x = Y,deriv.order = 1,y = xgrid, kernel = kernel)$est.fx/l_kde
d1_exact = -(xgrid-mu)/sigma2

# # second derivative of density
# d2_kde = kdde(Y,deriv.order = 2,eval.points = xgrid)$estimate
# d2_exact = sapply(xgrid,function(y)hessian(function(x)dnorm(x,mu,sqrt(sigma2)),y))
# second derivative of log-density
d2_kde = dkde(x = Y,deriv.order = 2,y = xgrid, kernel = kernel)$est.fx/l_kde - (d1_kde/l_kde)^2
d2_exact = -1/sigma2


g0 = ggplot() +
  geom_histogram(data = data.frame(Y = Y), aes(Y,..density..),alpha=0.6) +
  geom_line(aes(xgrid,l_kde),col="blue") +
  geom_line(aes(xgrid,l_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)

g0bis = ggplot() +
  geom_line(aes(xgrid,(ll_kde - ll_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed")

g1 = ggplot() +
  geom_line(aes(xgrid,d1_kde),col="blue") +
  geom_line(aes(xgrid,d1_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)

g1bis = ggplot() +
  geom_line(aes(xgrid,abs(d1_kde - d1_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()

g2 = ggplot() +
  geom_line(aes(xgrid,d2_kde),col="blue") +
  geom_line(aes(xgrid,d2_exact),col="red",size=2,linetype="dashed") +
  geom_vline(xintercept = Yt)

g2bis = ggplot() +
  geom_line(aes(xgrid,abs(d2_kde - d2_exact)),col="blue") +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_y_log10()

grid.arrange(g0,g0bis,g1,g1bis,g2,g2bis,nrow=3)

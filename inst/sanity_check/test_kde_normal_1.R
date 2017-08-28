##################################################################################################
# TOy example: Kernel density estimation for log-derivatives of the predictive distributions
##################################################################################################
rm(list = ls())
library(ks)
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
Y = rnorm(Ny,mu,sqrt(sigma2))


##### compute density estimators
xgrid = seq(-3,3,0.1)
# exact log-density
ll_exact = dnorm(xgrid,mu,sqrt(sigma2),log=TRUE)
l_exact = exp(ll_exact)
# via KDE
l_kde = kdde(Y,deriv.order = 0,eval.points = xgrid)$estimate
ll_kde = log(l_kde)
# # first derivative of density
d1_kde = kdde(Y,deriv.order = 1,eval.points = xgrid)$estimate
d1_exact = grad(function(x)dnorm(x,mu,sqrt(sigma2)),xgrid)
# first derivative of log-density
d1ll_kde = d1_kde/l_kde
d1ll_exact = -(xgrid-mu)/sigma2

# # second derivative of density
d2_kde = kdde(Y,deriv.order = 2,eval.points = xgrid)$estimate
d2_exact = sapply(xgrid,function(y)hessian(function(x)dnorm(x,mu,sqrt(sigma2)),y))
# second derivative of log-density
d2ll_kde = d2_kde/l_kde - (d1_kde/l_kde)^2
d2ll_exact = -1/sigma2

# ll_ = function(x)log(kdde(Y,deriv.order = 0,eval.points = x)$estimate)
# d1_kde = grad(ll_,xgrid)
# d2_kde = sapply(xgrid,function(y)hessian(ll_,y))


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
grid.arrange(g0,g0bis,g1ll,g1llbis,g2ll,g2llbis,nrow=3)


## direct derivative estimation method

Ny = 10^3
Y = rnorm(Ny,mu,sqrt(sigma2))

regularizer <- 1
sigma <- 1
psi <- function(ytilde, y) (ytilde - y) / (sigma^2) * exp(-(ytilde - y)^2 / (2*sigma^2))
dpsi <- function(ytilde, y) (-1/(sigma^2) + (ytilde - y)^2/(sigma^4)) * exp(-(ytilde - y)^2 / (2*sigma^2))


G <- matrix(0, Ny, Ny)
for (i in 1:Ny){
  psi_Yi <- sapply(X = 1:Ny, function(k) psi(Y[k], Y[i]))
  G <- G + matrix(psi_Yi, ncol = 1) %*% matrix(psi_Yi, nrow = 1)
}
G <- G / Ny
h <- rep(0, Ny)
for (i in 1:Ny){
  phi_Yi <- sapply(X = 1:Ny, function(k) dpsi(Y[k], Y[i]))
  h <- h + phi_Yi
}
h <- h / Ny

theta_hat <- -solve(G + regularizer * diag(1, Ny, Ny)) %*% h

other_d1ll <- function(x) t(theta_hat) %*% matrix(sapply(X = 1:Ny, function(k) psi(Y[k], x)), ncol = 1)
other_d1lls <- sapply(xgrid, other_d1ll)
ggplot() + geom_line(aes(xgrid,d1ll_kde),col="blue") +
  geom_line(aes(xgrid,d1ll_exact),col="red",size=2,linetype="dashed") +
  geom_line(aes(xgrid,other_d1lls),col="orange") +
  geom_vline(xintercept = Yt)


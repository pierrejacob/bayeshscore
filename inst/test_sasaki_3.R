
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

Yt = rnorm(1)


Ny = 10^7

# rtest = rnorm
# dtest = dnorm
# ptest = pnorm
rtest = rcauchy
dtest = dcauchy
ptest = pcauchy

# epsilon = 100
# upper <- Yt + epsilon
# lower <- Yt - epsilon
# Y = rep(NA,Ny)
# for (i in 1:Ny){
#   draw = rtest(1)
#   while ((draw < lower)||(draw > upper)){
#     draw = rtest(1)
#   }
#   Y[i] = draw
# }
# normalizing = ptest(upper)-ptest(lower)
Y = rtest(Ny)
normalizing = 1


d0_est = get_derivative_RBFlocal(Y,Yt,0.001,order = 0)
d1_est = get_derivative_RBFlocal(Y,Yt,0.002,order = 1)
d2_est = get_derivative_RBFlocal(Y,Yt,0.01,order = 2)


xgrid = seq(-4,4,0.01)


ggplot() +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = d0_est) +
  geom_line(aes(xgrid,dtest(xgrid)/normalizing),linetype="dashed",col="red",size=2)
ggplot() +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = d1_est) +
  geom_line(aes(xgrid,grad(function(y)dtest(y),xgrid)/normalizing),linetype="dashed",col="red",size=2)
ggplot() +
  geom_vline(xintercept = Yt) +
  geom_hline(yintercept = d2_est) +
  geom_line(aes(xgrid,sapply(xgrid,function(x)hessian(function(y)dtest(y),x))/normalizing),linetype="dashed",col="red",size=2)

# density
cat("0th d estim =",d0_est,"\n")
cat("0th d exact =",dtest(Yt),"\n")
# first derivative
cat("1st d estim =",d1_est,"\n")
cat("1st d exact =",grad(function(y)dtest(y),Yt),"\n")
# second derivative
cat("2nd d estim =",d2_est,"\n")
cat("2nd d exact =",hessian(function(y)dtest(y),Yt),"\n")
# first log-derivative
cat("1st log-d estim =",d1_est/d0_est,"\n")
cat("1st log-d exact =",grad(function(y)dtest(y),Yt)/dtest(Yt),"\n")
# second log-derivative
cat("2nd log-d estim =",d2_est/d0_est,"\n")
cat("2nd log-d exact =",hessian(function(y)dtest(y),Yt)/dtest(Yt),"\n")
# hscore
cat("hscore estim =",2*d2_est/d0_est - (d1_est/d0_est)^2,"\n")
cat("hscore exact =",2*hessian(function(y)dtest(y),Yt)/dtest(Yt) - (grad(function(y)dtest(y),Yt)/dtest(Yt))^2,"\n")


# # DIAGNOSTIC: Monte Carlo estimate is not so good because the importance region of the function
# # is much narrower than the sampling region ....
# d2K_values = d2K(Y,Yt,0.01)
# sum(d2K_values>10^-4)/Ny*100
# hist(d2K_values,breaks=100)
# # Note: the scaling factor is arbitrary and for visualization purpose only
# qplot(Y, geom = "histogram", bins = 1000) +
#   geom_area(aes(xgrid,15*abs(d2K(xgrid,Yt,0.01))),fill="red",alpha=0.6)

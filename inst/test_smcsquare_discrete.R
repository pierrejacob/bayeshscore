rm(list = ls())
library(HyvarinenSSM)
library(gridExtra)
set.seed(17)

# Define model and data
observations <- data_kangaroo[c("y1","y2"),]
rangeprior = 10
model1 <- get_model_kangarooLogistic(rangeprior)
model2 <- get_model_kangarooExponential(rangeprior)
model3 <- get_model_kangarooRandomwalk(rangeprior)

# Define initial proposal for theta (to avoid sampling from vague prior)
# Model 1
rinitial_theta1 = function(Ntheta) {
  sigma = runif(Ntheta,0,2)
  tau = runif(Ntheta,0,0.15)
  r = runif(Ntheta,0,8)
  b = runif(Ntheta,0,0.02)
  return (cbind(sigma,tau,r,b))
}
dinitial_theta1 = function(theta, log = TRUE){
  sigma = theta[1]
  tau = theta[2]
  r = theta[3]
  b = theta[4]
  if (log==TRUE){
    return (dunif(sigma,0,2,log) + dunif(tau,0,0.15,log) + dunif(r,0,8,log) + dunif(b,0,0.02,log))
  }
  else{
    return (dunif(sigma,0,2,log) * dunif(tau,0,0.15,log) * dunif(r,0,8,log) * dunif(b,0,0.02,log))
  }
}
# Model 2
rinitial_theta2 = function(Ntheta){
  sigma = runif(Ntheta,0,2)
  tau = runif(Ntheta,0,0.15)
  r = runif(Ntheta,0,8)
  return (cbind(sigma,tau,r))
}
dinitial_theta2 = function(theta, log = TRUE){
  sigma = theta[1]
  tau = theta[2]
  r = theta[3]
  if (log==TRUE){
    return (dunif(sigma,0,2,log) + dunif(tau,0,0.15,log) + dunif(r,0,8,log))
  }
  else{
    return (dunif(sigma,0,2,log) * dunif(tau,0,0.15,log) * dunif(r,0,8,log))
  }
}
# Model 3
rinitial_theta3 = function(Ntheta){
  sigma = runif(Ntheta,0,2)
  tau = runif(Ntheta,0,0.15)
  return (cbind(sigma,tau))
}
dinitial_theta3 = function(theta, log = TRUE){
  sigma = theta[1]
  tau = theta[2]
  if (log==TRUE){
    return (dunif(sigma,0,2,log) + dunif(tau,0,0.15,log))
  }
  else{
    return (dunif(sigma,0,2,log) * dunif(tau,0,0.15,log))
  }
}
# Define algorithmic parameters for each model
common = list(Ntheta = 2^10, Nx = 2^10,
              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
              progress = TRUE)
algorithmic_parameters1 <- common
algorithmic_parameters1$rinitial_theta = rinitial_theta1
algorithmic_parameters1$dinitial_theta = dinitial_theta1
algorithmic_parameters2 <- common
algorithmic_parameters2$rinitial_theta = rinitial_theta2
algorithmic_parameters2$dinitial_theta = dinitial_theta2
algorithmic_parameters3 <- common
algorithmic_parameters3$rinitial_theta = rinitial_theta3
algorithmic_parameters3$dinitial_theta = dinitial_theta3


<<<<<<< HEAD
r1 = hscore_discrete(observations, model1, algorithmic_parameters1)
#takes time !!! (about 10h with 2048 particles in both theta and X)
=======
r1 = hscore_discrete(observations, model1, algorithmic_parameters1) #takes time !!!
>>>>>>> a4d4748b10e61093c13e10230fc4ec164a6d4cdb
# r1 = list(hscore=0,logevidence=0) #run that instead to allow plot while skipping the logistic model
r2 = hscore_discrete(observations, model2, algorithmic_parameters2)
r3 = hscore_discrete(observations, model3, algorithmic_parameters3)



#sanity check ESS
par(mfrow=c(3,1))
plot(r1$ESS,type='l')
plot(r2$ESS,type='l')
plot(r3$ESS,type='l')
par(mfrow=c(1,1))

#sanity check posterior of theta
p11 = qplot(x = r1$thetas[,3], weight = r1$thetanormw, geom = "histogram",xlim = c(-2,10),xlab = "r",binwidth=0.7)
p12 = qplot(x = r1$thetas[,4], weight = r1$thetanormw, geom = "histogram",xlim = c(0,0.02),xlab = "b",binwidth=0.0015)
p21 = qplot(x = r1$thetas[,1], weight = r1$thetanormw, geom = "histogram",xlim = c(0,2),xlab = "sigma",binwidth=0.17)
p22 = qplot(x = r1$thetas[,2], weight = r1$thetanormw, geom = "histogram",xlim = c(0,0.15),xlab = "tau",binwidth=0.01)
grid.arrange(p11,p12,p21,p22, ncol = 2, nrow = 2)

p11 = qplot(x = r1$thetas[,3], weight = r1$thetanormw, geom = "histogram",xlim = c(-2,10),xlab = "r")
p12 = qplot(x = r1$thetas[,4], weight = r1$thetanormw, geom = "histogram",xlim = c(0,0.02),xlab = "b")
p21 = qplot(x = r1$thetas[,1], weight = r1$thetanormw, geom = "histogram",xlim = c(0,2),xlab = "sigma")
p22 = qplot(x = r1$thetas[,2], weight = r1$thetanormw, geom = "histogram",xlim = c(0,0.15),xlab = "tau")
grid.arrange(p11,p12,p21,p22, ncol = 2, nrow = 2)

#sanity check posterior of theta
p11 = qplot(x = r1$thetas[,3], weight = r1$thetanormw, geom = "density",adjust = 3, xlim = c(-2,10),xlab = "r")
p12 = qplot(x = r1$thetas[,4], weight = r1$thetanormw, geom = "density",adjust = 3, xlim = c(0,0.02),xlab = "b")
p21 = qplot(x = r1$thetas[,1], weight = r1$thetanormw, geom = "density",adjust = 3, xlim = c(0,2),xlab = "sigma")
p22 = qplot(x = r1$thetas[,2], weight = r1$thetanormw, geom = "density",adjust = 3, xlim = c(0,0.15),xlab = "tau")
grid.arrange(p11,p12,p21,p22, ncol = 2, nrow = 2)

p11 = qplot(x = r1$thetas[,3], weight = r1$thetanormw, geom = "density", xlim = c(-2,10),xlab = "r")
p12 = qplot(x = r1$thetas[,4], weight = r1$thetanormw, geom = "density", xlim = c(0,0.02),xlab = "b")
p21 = qplot(x = r1$thetas[,1], weight = r1$thetanormw, geom = "density", xlim = c(0,2),xlab = "sigma")
p22 = qplot(x = r1$thetas[,2], weight = r1$thetanormw, geom = "density", xlim = c(0,0.15),xlab = "tau")
grid.arrange(p11,p12,p21,p22, ncol = 2, nrow = 2)

#sanity check hscore
m = min(r1$hscore,r2$hscore,r3$hscore)
M = max(r1$hscore,r2$hscore,r3$hscore)
plot(r3$hscore,type='l',ylim = c(m,M))
points(r3$hscore,pch=17,col=1,lwd=1)
lines(r2$hscore)
points(r2$hscore,pch=19,col=1,lwd=1)
lines(r1$hscore)
points(r1$hscore,pch=15,col=1,lwd=1)
legend("bottomleft",legend=c("M3 - Random Walk","M2 - Exponential","M1 - Logistic"),pch=c(17,19,15))

#sanity check logevidence
m = min(r1$logevidence,r2$logevidence,r3$logevidence)
M = max(r1$logevidence,r2$logevidence,r3$logevidence)
plot(r3$logevidence,type='l',ylim = c(m,M))
points(r3$logevidence,pch=17,col=1,lwd=1)
lines(r2$logevidence)
points(r2$logevidence,pch=19,col=1,lwd=1)
lines(r1$logevidence)
points(r1$logevidence,pch=15,col=1,lwd=1)
legend("bottomleft",legend=c("M3 - Random Walk","M2 - Exponential","M1 - Logistic"),pch=c(17,19,15))

print(r1$logevidence)
print(r2$logevidence)
print(r3$logevidence)

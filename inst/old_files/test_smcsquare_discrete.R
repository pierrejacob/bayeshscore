rm(list = ls())
library(HyvarinenSSM)
library(gridExtra)
set.seed(17)

# Define model and data
observations <- data_kangaroo[c("y1","y2"),]
rangeprior = 10
model1 <- get_model_kangarooLogistic(rangeprior)
model2 <- get_model_kangarooExponential()
model3 <- get_model_kangarooRandomwalk()

# Define initial proposal for theta (to avoid sampling from vague prior)
# Define algorithmic parameters for each model
algorithmic_parameters = list(Ntheta = 2^10, Nx = 2^10,
              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
              progress = TRUE)


r1 = hscore_discrete(observations, model1, algorithmic_parameters) #takes time !!!
# r1 = list(hscore=0,logevidence=0) #run that instead to allow plot while skipping the logistic model
r2 = hscore_discrete(observations, model2, algorithmic_parameters)
r3 = hscore_discrete(observations, model3, algorithmic_parameters)



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

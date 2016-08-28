rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

observations <- data_kangaroo[c("y1","y2"),]
model1 <- get_model_kangarooLogistic()
model2 <- get_model_kangarooExponential()
model3 <- get_model_kangarooRandomwalk()
algorithmic_parameters <- list(Ntheta = 2^10, Nx = 2^10,
                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               progress = TRUE)
# result1 = hscore_discrete(observations, model1, algorithmic_parameters)
result2 = hscore_discrete(observations, model2, algorithmic_parameters)
result3 = hscore_discrete(observations, model3, algorithmic_parameters)


#sanity check posterior of theta (first component)#
setmytheme()
# qplot(x = result1$thetas[,1], weight = result1$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result2$thetas[,1], weight = result2$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result3$thetas[,1], weight = result3$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))


#sanity check
# par(mfrow=c(3,1))
# plot(result3$ESS,type='l')
# plot(result3$logevidence,type='l')
plot(result3$hscore,type='l')
points(result3$hscore,pch=17,col=1,lwd=1)
lines(result2$hscore)
points(result2$hscore,pch=19,col=1,lwd=1)
# lines(result1$hscore)
# points(result1$hscore,pch=15,col=1,lwd=1)
legend("topleft",legend=c("M3 - Random Walk","M2 - Exponential","M1 - Logistic"),pch=c(17,19,15))
par(mfrow=c(1,1))

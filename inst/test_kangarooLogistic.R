rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

nobservations <- 20
model <- get_model_kangarooLogistic()
#simulate data
sim = simulateData(model, theta = c(0.1,0.1,0.1,1), nobservations)
X = sim$X
Y = sim$Y

plot(1:nobservations,X,ylim = c(min(X,Y),max(X,Y)),type = 'l')
points(1:nobservations,Y[1,])
points(1:nobservations,Y[2,])


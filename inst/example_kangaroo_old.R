#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file applies the Hyvarinen score to kangaroos data (Knape et al. 2011 and 2012)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
setwd("C:/Users/jianfushao/Documents/Harvard/MyHarvard/Research/Model Selection for SSM (Pierre and Jie)")

# load the different functions
source("functions_model.R")
source("functions_hscore_discrete.R")

# load the different models (NOTE: the priors are specified inside the model file)
source("model_kangaroos_logistic.R")
source("model_kangaroos_randomwalk.R")
source("model_kangaroos_exponentialgrowth.R")

# Load Data from Caughley et al. 1987.
dataset = read.table("data_kangaroos.txt")
Y = dataset[1:2,]
timestep = unlist(dataset[3,])
T = length(timestep)

# # Simulate Data
# sim = simulateDataKangaroo(model.kangarooLogistic,c(0.1,1,1,0.1),timestep)
# X = sim$X
# Y = sim$Y
# plot(timestep,X,type='l',ylim = c(min(X,Y),max(X,Y)))
# points(timestep,Y[1,])
# points(timestep,Y[2,])


#---------- Compute the prequential Hyvarinen score for each model ----------#
Nx = 100
Ntheta = 100
hscoreM1 = hyvarinenDiscrete(Y,T,model.kangarooLogistic,Ntheta,Nx,timestep)
hscoreM2 = hyvarinenDiscrete(Y,T,model.kangarooExponential,Ntheta,Nx,timestep)
hscoreM3 = hyvarinenDiscrete(Y,T,model.kangarooRandomwalk,Ntheta,Nx,timestep)

#--------- Plot the prequential Hyvarinen scores for each model -------------#
ymin = min(hscoreM1$Hscore,hscoreM2$Hscore,hscoreM3$Hscore)
ymax = max(hscoreM1$Hscore,hscoreM2$Hscore,hscoreM3$Hscore)
plot(hscoreM1$Hscore,type = 'l',xlab = "Time horizon T",ylab = "SH(T)",ylim = c(ymin,ymax))
points(hscoreM1$Hscore,pch=15,col=1,lwd=1)
lines(hscoreM2$Hscore)
points(hscoreM2$Hscore,pch=19,col=1,lwd=1)
lines(hscoreM3$Hscore)
points(hscoreM3$Hscore,pch=17,col=1,lwd=1)
legend("bottomleft",legend=c("M1 - Logistic","M2 - Exponential","M3 - Random Walk"),pch=c(15,19,17))


rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(9)

#--------------------------------------------------------------------------------------------
# create data
nobservations <- 15
# model <- get_model_simplerlineargaussian()
# theta_star <- c(0.8,1,model$psi,model$sigmaV2)
model <- get_model_lineargaussian()
theta_star <- c(0.8,1,1,1)

sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^6
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.3
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL # this forces the use of SMC2
model_withoutlikelihood$dpredictive = NULL # this forces the use of SMC2

#--------------------------------------------------------------------------------------------

repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel

# Run SMC2
print(paste("Started at:",Sys.time()))
time_start = proc.time()
results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
  TreeClass <<- module_tree$Tree
  hscore(observations, model_withoutlikelihood, algorithmic_parameters)
}
time_end = proc.time()-time_start
cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(algorithmic_parameters$Ntheta),
          ", Nx = ",toString(algorithmic_parameters$Nx),"\n",sep = ""))
print(time_end)

#--------------------------------------------------------------------------------------------
# Combine and rearrange results
results.df = data.frame()
posterior.df = data.frame()
for (r in 1:repl){
  results.df = rbind(results.df, data.frame(time = 1:nobservations,
                                            logevidence = results[[r]]$logevidence,
                                            hscore = results[[r]]$hscore,
                                            rep = r,
                                            from = factor("smc2")))
  posterior.df = rbind(posterior.df, data.frame(theta1 = results[[r]]$thetas_history[[nobservations+1]][1,],
                                                theta2 = results[[r]]$thetas_history[[nobservations+1]][2,],
                                                theta3 = results[[r]]$thetas_history[[nobservations+1]][3,],
                                                theta4 = results[[r]]$thetas_history[[nobservations+1]][4,],
                                                rep = r,
                                                from = factor("smc2")))
}

#--------------------------------------------------------------------------------------------
# Plot results
ggplot(results.df) + geom_line(aes(time, logevidence/time, color = from, group = rep))
ggplot(results.df) + geom_line(aes(time, hscore/time, color = from, group = rep))
ggplot(posterior.df) + geom_density(aes(theta1, fill = from, group = rep), alpha = 0.6)









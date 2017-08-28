library(doParallel)
library(HyvarinenSSM)
library(gridExtra)

# Define model and data
nobservations <- 50
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)# observations in a matrix of dimensions dimy x nobservations

# Plot data
observations.df = data.frame(time = 1:nobservations, X = t(X),Y = t(Y))
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=Y),size=2) +
  geom_line(aes(y=X),linetype=2) +
  xlab("\n Time")
plot(g)


# Define algorithmic parameters for each model
Ntheta = 2^12
Nx = 2^8
algorithmic_parameters = list(Ntheta = Ntheta, Nx = Nx,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE)

repl = 6

registerDoParallel(cores=3)
print(paste("Started at:",Sys.time()))
time_start = proc.time()
results = foreach(i=1:repl,.packages='HyvarinenSSM',.verbose = TRUE) %dopar% {
  hscore_continuous(observations, model, algorithmic_parameters)
}
time_end = proc.time()-time_start
cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
          ", Nx = ",toString(Nx),"\n",sep = ""))
print(time_end)


results.df = data.frame()
posterior.df = data.frame()
for (r in 1:repl){
  results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                            logevidence = results[[r]]$logevidence,
                                            hscore = results[[r]]$hscore,
                                            ess = results[[r]]$ESS,
                                            rep = r))
  posterior.df = rbind(posterior.df, data.frame(phi = results[[r]]$thetas[,1],
                                                psi = results[[r]]$thetas[,2],
                                                sigmaV2 = results[[r]]$thetas[,3],
                                                sigmaW2 = results[[r]]$thetas[,4],
                                                w = results[[r]]$thetanormw,
                                                rep=r))
}

#plot ESS
g <- ggplot(results.df, aes(x = time, y = ess, group = rep)) + geom_line(colour = "blue")
plot(g)


#plot Bayes factor
g <- ggplot(results.df, aes(x = time, y = logevidence, group = rep)) + geom_line(colour = "blue")
plot(g)

#plot Hscore
g <- ggplot(results.df, aes(x = time, y = hscore, group = rep)) + geom_line(colour = "blue")
plot(g)

#plot Posterior
g11 <- ggplot(posterior.df, aes(x = phi, weight = w, group = rep)) + geom_density()
g12 <- ggplot(posterior.df, aes(x = psi, weight = w, group = rep)) + geom_density()
g21 <- ggplot(posterior.df, aes(x = sigmaV2, weight = w, group = rep)) + geom_density()
g22 <- ggplot(posterior.df, aes(x = sigmaW2, weight = w, group = rep)) + geom_density()
grid.arrange(g11,g12,g21,g22, ncol = 2, nrow = 2)

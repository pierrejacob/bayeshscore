# library(doParallel)
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

print(paste("Started at:",Sys.time()))
time_start = proc.time()
Rprof(tmp <- tempfile())
results <-  hscore_continuous(observations, model, algorithmic_parameters)
Rprof()
summaryRprof(tmp)
unlink(tmp)

time_end = proc.time()-time_start
cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
          ", Nx = ",toString(Nx),"\n",sep = ""))
aprint(time_end)
#
results.df = data.frame(time = 1:ncol(observations),
                                          logevidence = results$logevidence,
                                          hscore = results$hscore,
                                          ess = results$ESS)
posterior.df = data.frame(phi = results$thetas[,1],
                                              psi = results$thetas[,2],
                                              sigmaV2 = results$thetas[,3],
                                              sigmaW2 = results$thetas[,4],
                                              w = results$thetanormw)
g <- ggplot(results.df, aes(x = time, y = ess)) + geom_line(colour = "blue") + ylim(0, algorithmic_parameters$Ntheta)
plot(g)


#plot Bayes factor
g <- ggplot(results.df, aes(x = time, y = logevidence)) + geom_line(colour = "blue")
plot(g)

#plot Hscore
g <- ggplot(results.df, aes(x = time, y = hscore)) + geom_line(colour = "blue")
plot(g)

#plot Posterior
g11 <- ggplot(posterior.df, aes(x = phi, weight = w)) + geom_density()
g12 <- ggplot(posterior.df, aes(x = psi, weight = w)) + geom_density()
g21 <- ggplot(posterior.df, aes(x = sigmaV2, weight = w)) + geom_density()
g22 <- ggplot(posterior.df, aes(x = sigmaW2, weight = w)) + geom_density()
grid.arrange(g11,g12,g21,g22, ncol = 2, nrow = 2)


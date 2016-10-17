library(doParallel)
library(HyvarinenSSM)
library(gridExtra)

# Define model and data
nobservations <- 100
model <- get_model_lineargaussian_iid()
sim = simulateData(model, theta = c(1), nobservations)
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
Ntheta = 2^7
Nx = 5
algorithmic_parameters = list(Ntheta = Ntheta, Nx = Nx,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE)


# results = list(hscore_continuous(observations, model, algorithmic_parameters))
# repl = 1

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
                                            rep = r,
                                            sim = 1))
  posterior.df = rbind(posterior.df, data.frame(sigmaV2 = results[[r]]$thetas[,1],
                                                w = results[[r]]$thetanormw,
                                                rep = r,
                                                sim = 1))
}

#plot ESS
g <- ggplot(results.df, aes(x = time, y = ess, group = rep)) + geom_line(colour = "blue")
plot(g)


#plot Posterior
nu = model$df
nu_post = nu + nobservations
s2_post = (nu+sum(Y^2))/nu_post
g11 <- ggplot(posterior.df, aes(x = sigmaV2, weight = w, group = rep)) +
  stat_function(fun = function(y)dinvchisq(y,nu_post,s2_post,FALSE), colour = "red",size=1.5,linetype=1) +
  geom_density() + ylab("")
grid.arrange(g11,ncol = 1, nrow = 1)

#Compute exact h-score
hscore_exact = rep(NA,nobservations)
for (t in 1:nobservations){
  s = sum(Y[,1:t]^2)
  hscore_exact[t] = ((nu+t)/((nu+s)^2))*((nu+t+4)*Y[,t]^2-2*(nu+s))
}
hscore_exact = cumsum(hscore_exact)

#compute exact log-evidence
logevidence_exact = rep(NA,nobservations)
for (t in 1:nobservations) {
  nu_t = nu + (t-1)
  st2 = (nu + sum(Y[,1:(t-1)]^2))/nu_t
  logevidence_exact[t] = dtscaled(Y[,t],nu_t,st2,TRUE)
}
logevidence_exact = cumsum(logevidence_exact)


#Append exact results to dataframe
results.df = rbind(results.df, data.frame(time = 1:nobservations,
                                          logevidence = logevidence_exact,
                                          hscore = hscore_exact,
                                          ess = rep(1,nobservations),
                                          rep = 0,
                                          sim = 0))
#CHECK LOG-EVIDENCE
g <- ggplot(subset(results.df,sim>0.5), aes(x = time, y = logevidence, group = rep)) + geom_line() +
  geom_line(aes(x = time, y=logevidence),data = subset(results.df,sim<0.5),colour = "red",size=1.5,linetype=2)
plot(g)

#CHECK HSCORE
g <- ggplot(subset(results.df,sim>0.5), aes(x = time, y = hscore, group = rep)) + geom_line() +
  geom_line(aes(x = time, y=hscore),data = subset(results.df,sim<0.5),colour = "red",size=1.5,linetype=2) +
  ylab("prequential hscore")
plot(g)

#ABSOLUTE DIFFERENCE ESTIMATED VS EXACT LOG-EVIDENCE
g <- ggplot(subset(results.df,sim>0.5), aes(x = time, y = abs(logevidence-logevidence_exact), group = rep)) +
  geom_line(colour="blue") + geom_hline(aes(yintercept=0),linetype=2)
plot(g)

#ABSOLUTE DIFFERENCE ESTIMATED VS EXACT HSCORE
g <- ggplot(subset(results.df,sim>0.5), aes(x = time, y = abs(hscore-hscore_exact), group = rep)) +
  geom_line(colour="blue") + geom_hline(aes(yintercept=0),linetype=2)
plot(g)

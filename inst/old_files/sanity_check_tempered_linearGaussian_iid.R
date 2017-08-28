rm(list = ls())
library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree

# Define model and data
nobservations <- 50
model <- get_model_lineargaussian_iid()
true_sigmav2 = 1
sim = simulateData(model, theta = c(true_sigmav2), nobservations)
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
algorithmic_parameters = list(Ntheta = Ntheta,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE, nmoves = 1)
algorithmic_parameters$Nx = NULL # this uses adaptive Nx, starting with Nx = 128
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$ess_threshold = 0.5

repl = 5
results.df = data.frame()
posterior.df = data.frame()
results_temp.df = data.frame()
posterior_temp.df = data.frame()
for (r in 1:repl){
  # Tempered SMC^2
  results_temp = hscore_continuous(observations, model, algorithmic_parameters)
  # Regular SMC^2
  results = hscore_continuous_no_tempering(observations, model, algorithmic_parameters)
  results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                            logevidence = results$logevidence,
                                            hscore = results$hscore,
                                            ess = results$ESS,
                                            rep = r,
                                            sim = 1))
  results_temp.df = rbind(results_temp.df, data.frame(time = 1:ncol(observations),
                                            logevidence = results_temp$logevidence,
                                            hscore = results_temp$Hscore,
                                            rep = r,
                                            sim = 1))
  posterior.df = rbind(posterior.df, data.frame(sigmaV2 = results$thetas[,1],
                                                w = results$thetanormw,
                                                rep = r,
                                                sim = 1))
  posterior_temp.df = rbind(posterior_temp.df, data.frame(sigmaV2 = results_temp$thetas_history[[nobservations+1]][,1],
                                                w = results_temp$normw_history[[nobservations+1]],
                                                rep = r,
                                                sim = 1))
}


#plot Posterior
nu = model$df
nu_post = nu + nobservations
s2_post = (nu+sum(Y^2))/nu_post
g11 <- ggplot(data = posterior.df, aes(x = sigmaV2, weight = w, group = rep)) +
  stat_function(fun = function(y)dinvchisq(y,nu_post,s2_post,FALSE), colour = "red",size=1.5,linetype=1) +
  geom_density(color = "blue", linetype = "dashed") + ylab("") +
  geom_density(data = posterior_temp.df, aes(x = sigmaV2, weight = w, group = rep),color = "blue")
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

#CHECK LOG-EVIDENCE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = logevidence/time,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=3) +
  geom_line(aes(x = time, y = logevidence/time),data = subset(results.df,sim==1), colour="blue",size = 2,linetype="dashed") +
  geom_line(aes(x = time, y = logevidence/time, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)

#CHECK HSCORE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = hscore,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=3) +
  geom_line(aes(x = time, y = hscore),data = subset(results.df,sim==1), colour="blue",size = 2,linetype="dashed") +
  geom_line(aes(x = time, y = hscore, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)


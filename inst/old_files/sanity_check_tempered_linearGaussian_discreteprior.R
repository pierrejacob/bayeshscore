rm(list = ls())
library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
library(numDeriv)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree


# Define model and data
nobservations <- 50
model <- get_model_lineargaussian_discreteprior()
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
algorithmic_parameters = list(Ntheta = Ntheta,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE, nmoves = 1)
algorithmic_parameters$Nx = NULL # this uses adaptive Nx, starting with Nx = 128
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$ess_threshold = 0.5

repl = 1
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



#----------------------------------------------------------------------------------------------
kalmanfilter = function(Y,theta,initialmean,initialvar) {
  nobserv = ncol(Y)
  phi = theta[1]
  psi = theta[2]
  sigmaV2 = theta[3]
  sigmaW2 = theta[4]
  #initialize containers
  muX_t_t_1 = vector("numeric",nobserv) #contains the means of X_t given Y_1,...,Y_t-1
  muX_t_t = vector("numeric",nobserv) #contains the means of X_t given Y_1,...,Y_t
  PX_t_t_1 = vector("numeric",nobserv) #contains the variances of X_t given Y_1,...,Y_t-1
  PX_t_t = vector("numeric",nobserv) #contains the variances of X_t given Y_1,...,Y_t
  muY_t_t_1 = vector("numeric",nobserv) #contains the means of Y_t given Y_1,...,Y_t-1
  PY_t_t_1 = vector("numeric",nobserv) #contains the variances of Y_t given Y_1,...,Y_t-1
  #initialize recursion
  muX_t_t_1[1] = initialmean
  PX_t_t_1[1] = initialvar
  Kt = PX_t_t_1[1]*psi/(psi*PX_t_t_1[1]*psi + sigmaV2)
  muX_t_t[1] = muX_t_t_1[1] + Kt*(Y[,1] - psi*muX_t_t_1[1])
  PX_t_t[1] = (1-Kt*psi)*PX_t_t_1[1]
  muY_t_t_1[1] = psi*muX_t_t_1[1]
  PY_t_t_1[1] = (psi^2)*PX_t_t_1[1] + sigmaV2
  #iterate
  for (t in 2:nobserv) {
    muX_t_t_1[t] = phi*muX_t_t[t-1]
    PX_t_t_1[t] = (phi^2)*PX_t_t[t-1] + sigmaW2
    Kt = PX_t_t_1[t]*psi/(psi*PX_t_t_1[t]*psi + sigmaV2)
    muX_t_t[t] = muX_t_t_1[t] + Kt*(Y[,t] - psi*muX_t_t_1[t])
    PX_t_t[t] = (1-Kt*psi)*PX_t_t_1[t]
    muY_t_t_1[t] = psi*muX_t_t_1[t]
    PY_t_t_1[t] = (psi^2)*PX_t_t_1[t] + sigmaV2
  }
  return (list(muX_t_t_1 = muX_t_t_1, muX_t_t = muX_t_t, muY_t_t_1 = muY_t_t_1,
               PX_t_t_1 = PX_t_t_1, PX_t_t = PX_t_t, PY_t_t_1 = PY_t_t_1))
}

phi = model$phi
psi = model$psi
sigmaW2 = model$sigmaW2
X_t_t_1 = matrix(NA,nrow = length(model$supportprior),ncol = nobservations)
P_t_t_1 = matrix(NA,nrow = length(model$supportprior),ncol = nobservations)
for (i in 1:length(model$supportprior)){
  kalman = kalmanfilter(Y,c(phi,psi,model$supportprior[i],sigmaW2),model$initialmean,model$initialvar)
  X_t_t_1[i,] = kalman$muX_t_t_1
  P_t_t_1[i,] = kalman$PX_t_t_1
}


#Compute exact posterior
posterior_exact = matrix(NA,nrow = length(model$supportprior), ncol = nobservations)
for (i in 1:length(model$supportprior)){
  sigV2 =  model$supportprior[i]
  posterior_exact[i,1] = dnorm(Y[,1],psi*X_t_t_1[i,1],sqrt(sigV2 + (psi^2)*P_t_t_1[i,1]))
}
posterior_exact[,1] = (posterior_exact[,1])/sum(posterior_exact[,1])
for (t in 2:nobservations){
  for (i in 1:length(model$supportprior)){
    sigV2 =  model$supportprior[i]
    posterior_exact[i,t] = posterior_exact[i,t-1]*dnorm(Y[,t],psi*X_t_t_1[i,t],sqrt(sigV2 + (psi^2)*P_t_t_1[i,t]))
  }
  posterior_exact[,t] = (posterior_exact[,t])/sum(posterior_exact[,t])
}

#Compute exact log-predictive and evidence
log_py_t_t_1 = rep(NA,nobservations)
log_py_t_t_1_func = list()
log_py_t_t_1_func[[1]] = function(y){
  temp = 0
  for (i in 1:length(model$supportprior)){
    sigV2 =  model$supportprior[i]
    temp = temp + (1/length(model$supportprior))*dnorm(y,phi*X_t_t_1[i,1],sqrt(sigV2 + (psi^2)*P_t_t_1[i,1]))
  }
  return (log(temp))
}
log_py_t_t_1[1] = log_py_t_t_1_func[[1]](Y[,1])
for (t in 2:nobservations){
  log_py_t_t_1_func[[t]] = function(y){
    temp = 0
    for (i in 1:length(model$supportprior)){
      sigV2 =  model$supportprior[i]
      temp = temp + posterior_exact[i,t-1]*dnorm(y,phi*X_t_t_1[i,t],sqrt(sigV2 + (psi^2)*P_t_t_1[i,t]))
    }
    return (log(temp))
  }
  log_py_t_t_1[t] = log_py_t_t_1_func[[t]](Y[,t])
}
exact_logevidence = cumsum(log_py_t_t_1)


#Compute exact (numerically via NumDeriv) hscore
exact_hscore = rep(NA,nobservations)
for (t in 1:nobservations){
  exact_hscore[t] = 2*hessian(log_py_t_t_1_func[[t]],Y[,t]) + (grad(log_py_t_t_1_func[[t]],Y[,t]))^2
}
exact_preq_hscore = cumsum(exact_hscore)


#Append exact log-evidence
results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                          logevidence = exact_logevidence,
                                          hscore = exact_preq_hscore,
                                          ess = rep(1,nobservations),
                                          rep = 0,
                                          sim = 0))

#Check Posterior
dx = c(model$supportprior[1],diff(model$supportprior)) #used to renormalize discrete exact posterior (by approximating area by riemann sum)
g11 <- ggplot(subset(posterior.df,sim>0.5)) +
  geom_line(aes(x = sigmaV2, weight = w, group = rep, ..density..),stat = "density", color = "blue", linetype = "dashed") +
  geom_line(aes(x=x,y=y/dx),data = data.frame(x=model$supportprior,y=posterior_exact[,nobservations]),colour="red",size=1.5,linetype=1) +
  geom_density(data = posterior_temp.df, aes(x = sigmaV2, weight = w, group = rep),color = "blue")
grid.arrange(g11,ncol = 1, nrow = 1)

#CHECK LOG-EVIDENCE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = logevidence,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=3) +
  geom_line(aes(x = time, y = logevidence),data = subset(results.df,sim==1), colour="blue",size = 3,linetype="dashed") +
  geom_line(aes(x = time, y = logevidence, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)

#CHECK HSCORE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = hscore,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=3) +
  geom_line(aes(x = time, y = hscore),data = subset(results.df,sim==1), colour="blue",size = 3,linetype="dashed") +
  geom_line(aes(x = time, y = hscore, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)


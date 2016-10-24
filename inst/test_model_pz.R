rm(list = ls())
library(HyvarinenSSM)
library(doMC)

set.seed(17)
get_model_pz_other <- function(){

  model.pz = list()

  # dimension of parameter, observations, and states
  model.pz$dimtheta = 5
  model.pz$dimY = 1
  model.pz$dimX = 2

  # sampler from the prior distribution on parameters
  model.pz$rprior = function(Ntheta){
    theta1 <- runif(Ntheta)
    theta2 <- runif(Ntheta)
    theta3 <- runif(Ntheta)
    theta4 <- runif(Ntheta)
    theta5 <- runif(Ntheta)
    return (cbind(theta1,theta2,theta3,theta4,theta5))
  }

  # density the prior distribution on parameters
  model.pz$dprior = function(theta, log = TRUE){
    density <- dunif(theta[1],0,1, log = TRUE) + dunif(theta[2],0,1, log = TRUE) +
      dunif(theta[4],0,1, log = TRUE) + dunif(theta[4],0,1, log = TRUE) + dunif(theta[5],0,1, log = TRUE)
    if (log){
      return(density)
    } else {
      return(exp(density))
    }
  }

  # sampler from the initial distribution of the states
  model.pz$rinitial = function(theta, N){
    return (matrix(rlnorm(N*2,meanlog = log(2), sdlog = 1), ncol = 2))
  }

  model.pz$rtransition = function(Xt, t, theta){
    N <- nrow(Xt)
    alphas <- theta[1] + theta[2] * rnorm(N, 0, 1)
    xparticles <- pz_transition(Xt, alphas, t-1, c(theta[3:5], 0.1)) # the difference is here: m_q is set to 0
    return(xparticles)
  }

  # density of the observations
  model.pz$dobs = function(Yt,Xt,t,theta, log = TRUE){
    return(dnorm(Yt[1], mean = log(Xt[,1]), sd = 0.2, log = log))
  }

  model.pz$derivativelogdobs = function(Yt,Xt,t,theta,k){
    N = nrow(Xt)
    d1 = (log(Xt[,1])-matrix(Yt, N))/(0.2^2)
    d2 = matrix(-1/(0.2^2),nrow = N)
    return (list(d1log = d1, d2log = d2))
  }
  model.pz$robs = function(Xt,t,theta){
    N = nrow(Xt)
    return (log(Xt[,1]) + rnorm(N, mean = 0, sd = 0.2))
  }
  return(model.pz)
}
model.pz <- get_model_pz()
model.pz$theta = c(0.7,0.5,0.25,0.3)
nobservations <- 100
sim = simulateData(model.pz, theta = c(0.7,0.5,0.25,0.3), nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model.pz$dimY)# observations in a matrix of dimensions dimy x nobservations

model.pz_other <- get_model_pz_other()

# Plot data
observations.df = data.frame(time = 1:nobservations, X = t(log(X[1,,drop=FALSE])),Y = t(Y))
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=Y),size=2) +
  geom_line(aes(y=X),linetype=2) +
  xlab("\n Time")
plot(g)

algorithmc_parameters <- list(Nx = 2^10, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))
theta <- model.pz$theta

mus <- seq(from = 0, to = 1, length.out = 100)
registerDoMC(cores = 6)
loglikelihoods <- foreach(mu = mus, .combine = c) %dorng% {
  bpfresults <- bootstrap_particle_filter(observations, model.pz_other, c(mu, theta[2:4], 0.1), algorithmc_parameters)
  bpfresults$log_p_y_hat
}

plot(x = mus, y = loglikelihoods)


algorithmic_parameters <- list(Ntheta = 2^8, Nx = 2^8, nmoves = 3,
                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               progress = TRUE)


result = hscore_continuous(observations, model.pz_other, algorithmic_parameters)
result$rejuvenation_times
result$rejuvenation_accepts

#sanity check posterior of theta (first component)
length(unique(result$thetas[,4]))
setmytheme()
qplot(x = result$thetas[,1], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result$thetas[,2], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result$thetas[,3], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result$thetas[,4], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))
qplot(x = result$thetas[,5], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))


# registerDoParallel(cores=3)
nrep <- 3
results = foreach(i=1:nrep,.packages='HyvarinenSSM',.verbose = TRUE) %dorng% {
  hscore_continuous(observations, model.pz_other, algorithmic_parameters)
}

names(results[[1]])
hscores <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results[[irep]]$hscore
}
evidences <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results[[irep]]$logevidence
}

ESSs <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results[[irep]]$ESS
}

matplot(hscores, type = "l")
matplot(evidences, type = "l")
matplot(ESSs, type = "l")

length(unique(results[[1]]$thetas))
#sanity check
# par(mfrow=c(3,1))
# plot(result$ESS,type='l')
# plot(result$logevidence,type='l')
# plot(result$hscore,type='l')
# par(mfrow=c(1,1))

# print(result$logevidence)
# print(result$hscore)

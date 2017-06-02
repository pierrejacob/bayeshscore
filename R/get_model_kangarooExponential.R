##################################################################################################
# Kangaroo Model 2 - Exponential (Knape and Valpine, 2012)
##################################################################################################
#'@rdname get_model_kangarooExponential
#'@title get_model_kangarooExponential
#'@description Population model for counts of kangaroos: model 2 in Knape and Valpine (2012).
#'This model has 3 parameters: \code{sigma}, \code{tau}, \code{r}.
#'The priors are independent uniforms
#'on (0,\code{range_sigma}) for \code{sigma},
#'on (0,\code{range_tau}) for \code{tau},
#'on (-\code{range_r},\code{range_r}) for \code{r}.
#'@export
get_model_kangarooExponential = function(timesteps = data_kangaroo["time",],
                                      range_sigma = 10,
                                      range_tau = 10,
                                      range_r = 10){
  model = list()
  # Type of observations (string): "continuous" or "discrete"
  model$observation_type = "discrete"
  # Dimension of parameter, observations, and possibly latent states (int)
  model$dimtheta = 3
  model$dimY = 2
  model$dimX = 1
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # If the observations are discrete, lower and upper bounds for each component of the observations
  # are required to compute the Hyvarinen score.
  # model$lower and model$upper are vectors of length model$dimY, such that model$lower[j] and
  # model$upper[j] are respectively the lower and upper bound defining the support of the j-th
  # component of the observations
  # Note: if no bounds are provided, the lower/upper bounds are set to -Inf/Inf by default
  model$lower = c(0,0)
  model$upper = c(Inf,Inf)
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # Sampler from the prior distribution on parameters
  # inputs: Ntheta (int)
  # outputs: matrix (dimtheta by Ntheta) of prior draws
  model$rprior = function(Ntheta){
    sigma = runif(Ntheta,0,range_sigma)
    tau = runif(Ntheta,0,range_tau)
    r = runif(Ntheta,-range_r,range_r)
    return (rbind(sigma,tau,r))
  }
  # prior density on parameters
  # inputs: theta (single vector), log (TRUE by default)
  # outputs: prior (log)-density theta (double)
  model$dprior = function(theta, log = TRUE){
    llsigma = dunif(theta[1],0,range_sigma,TRUE)
    lltau = dunif(theta[2],0,range_tau,TRUE)
    llr = dunif(theta[3],-range_r,range_r,TRUE)
    lp = llsigma + lltau + llr
    if (log==TRUE) {return (lp)}
    else {return (exp(lp))}
  }
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # Note: if no likelihood nor predictive is provided, the method will be SMC2, which requires
  # specifying the transition kernel and the observation density
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # Sampler from the initial distribution of the latent states
  # inputs: theta (single vector), Nx (int)
  # outputs: matrix (dimX by Nx) of latent states
  model$rinitial = function(theta,Nx){
    return (matrix(rlnorm(Nx,meanlog = 5,sdlog = sqrt(10)),ncol = Nx))
  }
  # Sampler from the transition distribution of the latent states
  # inputs: current states Xs at time (t-1) (dimX by Nx matrix), time t (int), theta (single vector)
  # outputs: updated states (dimX by Nx)
  model$rtransition = function(Xs,t,theta){
    Nx = ncol(Xs)
    sigma = theta[1]
    r = theta[3]
    dt = timesteps[t] - timesteps[t-1]
    return (matrix(Xs*exp(r*dt+sigma*rnorm(Nx,mean = 0,sd = sqrt(dt))),ncol = Nx))
  }
  # observation density
  # inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector), log (TRUE by default)
  # outputs: observation (log)-densities ("vectorized" with respect to the states Xts)
  model$dobs = function(Yt,Xts,t,theta,log = TRUE){
    tau = theta[2]
    n = 1/tau
    ld = dnbinom(Yt[1],size = n,mu = Xts,log = TRUE) + dnbinom(Yt[2],size = n,mu = Xts,log = TRUE)
    if (log==TRUE){return (ld)}
    else{return (exp(ld))}
  }
  # sampler from the observation disctribution
  # inputs: states Xts (dimX by Nx), time t, theta (single vector), log (TRUE by default)
  # outputs: observations (dimY by Nx matrix)
  model$robs = function(Xt,t,theta){
    tau = theta[2]
    return (matrix(rnbinom(2, size = 1/tau, mu = Xt), nrow = model$dimY))
  }

  return(model)
}

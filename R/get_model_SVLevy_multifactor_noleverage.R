##################################################################################################
# Levy-driven Stochastic volatility model: Multi-factor without leverage
# (section 4.1 in Chopin, Jacob, Papaspiliopoulos, 2013)
##################################################################################################
#'@rdname get_model_SVLevy_multifactor_noleverage
#'@title get_model_SVLevy_multifactor_noleverage
#'@description Levy-driven Stochastic volatility model with multi-factor, as described by equations
#'(14) and (15) in Chopin, Jacob, Papaspiliopoulos (2013), and discussed in more detail
#'in Barndorff-Nielsen and Shephard (2002).
#'The states are Xt = (v1t, z1t, v2t, z2t).
#'@export
get_model_SVLevy_multifactor_noleverage <- function(timesteps,
                                                    mu0mu = 0, sigma02mu = 10,
                                                    mu0beta = 0, sigma02beta = 10,
                                                    r0xi = 1/5, r0w2 = 1/5,
                                                    r0lambda1 = 1, r0lambda2_1 = 1/2,
                                                    alpha0w = 1, beta0w = 1){
  model = list()
  # Type of observations (string): "continuous" or "discrete"
  model$observation_type = "continuous"
  # Dimension of parameter, observations, and possibly latent states (int)
  model$dimtheta = 7
  model$dimY = 1
  model$dimX = 4
  # Sampler from the prior distribution on parameters
  # inputs: Ntheta (int)
  # outputs: matrix (dimtheta by Ntheta) of prior draws
  model$rprior = function(Ntheta){
    mu = rnorm(Ntheta, mu0mu, sqrt(sigma02mu))
    beta = rnorm(Ntheta, mu0beta, sqrt(sigma02beta))
    xi = rexp(Ntheta, r0xi)
    w2 = rexp(Ntheta, r0w2)
    lambda1 = rexp(Ntheta, r0lambda1)
    lambda2 = lambda1 + rexp(Ntheta, r0lambda2_1)
    w = rbeta(Ntheta, alpha0w, beta0w)
    return (rbind(mu, beta, xi, w2, lambda1, lambda2, w))
  }
  # prior density on parameters
  # inputs: theta (single vector), log (TRUE by default)
  # outputs: prior (log)-density theta (double)
  model$dprior = function(theta, log = TRUE){
    lmu = dnorm(theta[1], mu0mu, sqrt(sigma02mu), log = TRUE)
    lbeta = dnorm(theta[2], mu0beta, sqrt(sigma02beta), log = TRUE)
    lxi = dexp(theta[3], r0xi, log = TRUE)
    lw2 = dexp(theta[4], r0w2, log = TRUE)
    llambda = dexp(theta[5], r0lambda1, log = TRUE) + dexp(theta[6]-theta[5], r0lambda2_1, log = TRUE)
    lw = dbeta(theta[7], alpha0w, beta0w, log = TRUE)
    lp = lmu + lbeta + lxi + lw2 + llambda + lw
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
    # xi = theta[3]
    # w2 = theta[4]
    # lambda1 = theta[5]
    # lambda2 = theta[6]
    # w = theta[7]
    return (rinitial_SVLevy_multifactor_cpp(Nx,timesteps[1],theta[3],theta[4],theta[5],theta[6],theta[7]))
  }
  # Sampler from the transition distribution of the latent states
  # inputs: current states Xs at time (t-1) (dimX by Nx matrix), time t (int), theta (single vector)
  # outputs: updated states (dimX by Nx)
  model$rtransition = function(Xs,t,theta){
    # xi = theta[3]
    # w2 = theta[4]
    # lambda1 = theta[5]
    # lambda2 = theta[6]
    # w = theta[7]
    return (rtransition_SVLevy_multifactor_cpp(Xs,timesteps[t-1],timesteps[t],theta[3],theta[4],theta[5],theta[6],theta[7]))
  }
  # observation density
  # inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector), log (TRUE by default)
  # outputs: observation (log)-densities ("vectorized" with respect to the states Xt)
  model$dobs = function(Yt,Xts,t,theta,log = TRUE){
    # mu = theta[1]
    # beta = theta[2]
    return (dobs_SVLevy_multifactor_cpp(Yt, Xts, theta[1], theta[2], log))
  }
  # OPTIONAL: first and second partial derivatives of the observation log-density
  # The function is vectorized with respect to the states Xts (dimX by Nx)
  # inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector)
  # outputs: list with the following fields
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  # NB: if missing, this field is automatically filled with numerical derivatives
  # via set_default_model in util_default.R)
  model$derivativelogdobs = function(Yt,Xts,t,theta){
    # mu = theta[1]
    # beta = theta[2]
    return (list(jacobian = d1logdobs_SVLevy_multifactor_cpp(Yt,Xts,theta[1],theta[2]),
                 hessiandiag = d2logdobs_SVLevy_multifactor_cpp(Yt,Xts,theta[1],theta[2])))
  }
  # sampler from the observation disctribution
  # inputs: single state Xt (dimX by 1), time t, theta (single vector)
  # outputs: single observation (dimY by 1 matrix)
  model$robs = function(Xt,t,theta){
    mu = theta[1]
    beta = theta[2]
    vt = Xt[1,] + Xt[3,]
    return (rnorm(1, mu + beta*vt, sqrt(vt)))
  }

  return(model)
}

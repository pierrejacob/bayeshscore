rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# calculation
# likelihood f(y|theta) = exp(-y/theta) * (1 / theta)
# log likelihood log f(y|theta) = -y/theta - log(theta)
# derivatives wrt y
# * 1st: y maps to -1/theta
# * 2nd: y maps to 0
# Jeffrey's prior: p(theta) proportional to 1/theta (equal sqrt of determinant of variance of score)

# xx <- rinvgamma(100000, 2.5, 2)
# hist(xx, prob = TRUE, nclass = 1000)
# curve(sapply(x, function(v) dinvgamma(v, 2.5, 2, log=F)), add = TRUE, col = "red", n = 10000)

##################################################################################################
# Generic template to define a model
##################################################################################################
get_exponential <- function(){
  model = list()
  # Type of observations (string): "continuous" or "discrete"
  model$observation_type = "continuous"
  # Dimension of parameter, observations, and possibly latent states (int)
  model$dimtheta = 1
  model$dimY = 1
  # model$dimX = 0
  # Sampler from the prior distribution on parameters
  # inputs: Ntheta (int)
  # outputs: matrix (dimtheta by Ntheta) of prior draws
  model$rprior = function(Ntheta){
    return (matrix(rinvgamma(Ntheta, 1, 1), nrow = 1))
  }
  # prior density on parameters
  # inputs: theta (single vector), log (TRUE by default)
  # outputs: prior (log)-density theta (double)
  model$dprior = function(theta, log = TRUE){
    return(dinvgamma(theta, 1, 1, log = log))
  }
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # Note: to use SMC, one may specify either the likelihood OR the one-step ahead predictive
  # (one is automatically filled given the other, via set_default_model in util_default.R)
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # OPTIONAL: one-step predicitve density of the observation at t given the past from 1 to (t-1) and theta
  # inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector),
  #         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter),
  #         log (TRUE by default)
  # outputs: log-likelihood of the observations from time 1 to t given theta (double)
  # WARNING: must be an explicit function of the observation at time t to allow the
  # computation of the derivative of the log-predictive density
  model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
    y <- observations[1,t]
    lp = -y/theta - log(theta)
    if (log) {return(lp)}
    else {return(exp(lp))}
  }
  # OPTIONAL: derivatives of the predicitve density with respect to the observation at time t
  # inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector),
  #         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter)
  # outputs: list with the following fields
  # jacobian >> the transpose of the gradient (1 by dimY)
  # hessiandiag >> the Hessian diagonal coefficients (1 by dimY)
  # NB: if missing, this field is automatically filled with numerical derivatives
  # via set_default_model in util_default.R)
  model$derivativelogdpredictive = function(observations,t,theta,byproduct) {
    y <- observations[1,t]
    deriv1 <- -1/theta
    deriv2 <- 0
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }
  return(model)
}


#--------------------------------------------------------------------------------------------
#--------------------------------- TREE MODULE (for SMC2) -----------------------------------
#--------------------------------------------------------------------------------------------
# if SMC2 will be used, need tree module (NB: automatically called upon loading HyvarinenSSM package)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
#---------------------------------        DATA            -----------------------------------
#--------------------------------------------------------------------------------------------
nobservations = 20
Y = rexp(nobservations, 1)
# observations in a matrix of dimensions dimy x nobservations
observations = matrix(Y, nrow = 1)

model <- get_exponential()


#--------------------------------------------------------------------------------------------
#---------------------               ALGORITHMIC PARAMETERS         -------------------------
#--------------------------------------------------------------------------------------------
# Unspecified fields are set to their default values (in parentheses)
# via set_default_algorithmic_parameters in util_default.R
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10 # <<<< Number of particles theta (2^7)
# algorithmic_parameters$ess_threshold = ... # <<<< ESS threshold for rejuvenation (0.5)
# algorithmic_parameters$nmoves = ... # <<<< Number of moves per rejuvenation step (1)
# algorithmic_parameters$Nx = ... # <<<< (ONLY for SMC2) Number of particles X (big O of nobservations)
# algorithmic_parameters$adaptNx = ... # <<<< (ONLY for SMC2) Adaptive number of particles X (TRUE)
# algorithmic_parameters$min_acceptance_rate = ... # <<<< (ONLY for SMC2 if adaptNx) Acceptance rae threshold to increase Nx (0.2)
# algorithmic_parameters$Nx_max = ... # <<<< (ONLY for SMC2 if adaptNx) Maximum value for Nx
# algorithmic_parameters$hscore = ... # <<<< if TRUE, computes the prequential hscore (TRUE)
algorithmic_parameters$verbose = TRUE # <<<< if TRUE, displays tempering steps and acceptance rates (TRUE)
# algorithmic_parameters$progress = ... # <<<< if TRUE, displays progress bar (FALSE)
# algorithmic_parameters$store_theta = ... # <<<< if TRUE, stores the history of particles thetas (TRUE for smc and smc2, FALSE for hscore)
# algorithmic_parameters$store_X = ... # <<<< (ONLY for SMC2) if TRUE, stores the history of particles X (FALSE)
# algorithmic_parameters$store_byproducts = ... # <<<< (ONLY for SMC) if TRUE, stores the history of byproducts (FALSE)
# algorithmic_parameters$save = ... # <<<< if TRUE, saves intermediary results in savefilename (FALSE)
# algorithmic_parameters$savefilename = ... # <<<< (ONLY if save) filename to save intermediary results (RDS file with time stamp in working directory)
# algorithmic_parameters$time_budget = ... # <<<< time budget for computation in number of seconds (NULL)
# algorithmic_parameters$resampling = ... # <<<< resampling function (systematic resampling)
# algorithmic_parameters$proposalmove = ... # <<<< proposal for rejuvenation moves (get_proposal_mixture())
#--------------------------------------------------------------------------------------------
#---------------------                HSCORE + LOGEVIDENCE          -------------------------
#--------------------------------------------------------------------------------------------
# Run SMC or SMC2
results = hscore(observations, model, algorithmic_parameters)
# exctract particles thetas and weights
if (results$algorithmic_parameters$store_theta) {
  thetas = results$thetas_history[[nobservations+1]]
  normw = results$normw_history[[nobservations+1]]
}
# exctract logevidence
logevidence = results$logevidence
# exctract hscore
if (results$algorithmic_parameters$hscore) {
  hscore = results$hscore
}
hscore


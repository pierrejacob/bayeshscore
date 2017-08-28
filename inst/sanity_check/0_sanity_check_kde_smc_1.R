##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in an iid Normal case (with conjugate prior so that everything can be computed analytically)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# Define model and data
nu0 = 1
s02 = 1
nobservations = 5
observations = matrix(rnorm(nobservations), nrow = 1)# observations in a matrix of dimensions dimy x nobservations
model = get_model_iid_gaussian_unknown_variance(nu0,s02)
model = set_default_model(model)
# CDF of a scaled t distribution
ptscaled = function(x,nu_t,st2) {
  return (pt(x/sqrt(st2),nu_t))
}
# #--------------------------------------------------------------------------------------------
# # set algorithmic parameters
# algorithmic_parameters = list()
# algorithmic_parameters$Ntheta = 2^10
# algorithmic_parameters$verbose = TRUE
# algorithmic_parameters$store_thetas_history = TRUE
#
# algorithmic_parameters$use_kde = TRUE
# algorithmic_parameters$parameters_kde = list(Ny = 2^10)
# algorithmic_parameters = set_default_algorithmic_parameters(observations,model,algorithmic_parameters)
# #--------------------------------------------------------------------------------------------
# ### Run SMC
# results = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#----------------------------------------------------------------------------------------
Ny = 10^3
tol = 1
hscore_kde = rep(NA,nobservations)
hscore_exact = rep(NA,nobservations)
hscore_numderiv = rep(NA,nobservations)

for (t in 1:nobservations){
  upper = observations[,t]+tol
  lower = observations[,t]-tol
  # Draw one Yt for each particle theta
  Yt_sim = matrix(NA, nrow = model$dimY, ncol = Ny)
  nu_t = nu0 + (t-1)
  if (t == 1){st2 = nu0/nu_t}
  if (t >= 2){st2 = (nu0+sum(observations[,1:(t-1)]^2))/nu_t}
  for (i in 1:Ny){
    Ydraw = rtscaled(1,nu_t,st2)
    while ((Ydraw < lower)||(Ydraw > upper)){
      Ydraw = rtscaled(1,nu_t,st2)
    }
    Yt_sim[,i] = Ydraw
  }
  # grid of values where to compute the density and its derivatives
  xgrid = seq(lower,upper,length.out = 100)
  # exact log-derivatives
  d1ll_exact = -(nu_t+1)*xgrid/(nu_t*st2 + xgrid^2)
  d2ll_exact = -(nu_t+1)*(nu_t*st2 - xgrid^2)/((nu_t*st2 + xgrid^2)^2)
  # numDeriv log-derivatives (as a sanity check of the previous analytical expressions)
  ll_exact_fn = function(x) {dtscaled(x,nu_t,st2,TRUE)-log(ptscaled(upper,nu_t,st2)-ptscaled(lower,nu_t,st2))}
  l_exact = exp(ll_exact_fn(xgrid))
  d1ll_numd = sapply(xgrid, function(x)grad(ll_exact_fn,x))
  d2ll_numd = sapply(xgrid, function(x)hessian(ll_exact_fn,x))
  # density and derivatives via ks package
  l_kde = kdde(t(Yt_sim),deriv.order = 0,eval.points = xgrid)$estimate
  d1_kde = kdde(t(Yt_sim),deriv.order = 1,eval.points = xgrid)$estimate
  d2_kde = kdde(t(Yt_sim),deriv.order = 2,eval.points = xgrid)$estimate
  # log-density and log-derivatives via ks package
  d1ll_kde = d1_kde/l_kde
  d2ll_kde = d2_kde/l_kde - (d1ll_kde)^2
  #################################### plots #########################################
  g0 = ggplot() + geom_histogram(aes(x=c(Yt_sim), y=..density..),bins = 100, alpha = 0.6) +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_line(aes(xgrid, l_kde), col = "blue") +
    geom_line(aes(xgrid, l_exact),col = "red",size=2,linetype="dashed")
  g0bis = ggplot() +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_line(aes(xgrid, abs(l_kde-l_exact)),col="blue") +
    scale_y_log10()
  g1 = ggplot() +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_line(aes(xgrid, d1ll_kde),col="blue") +
    geom_line(aes(xgrid, d1ll_numd), col = "red", linetype = "dashed", size = 2) +
    geom_line(aes(xgrid, d1ll_exact), col = "red")
  g1bis = ggplot() +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_line(aes(xgrid, abs(d1ll_kde-d1ll_exact)),col="blue") +
    scale_y_log10()
  g2 = ggplot() +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_line(aes(xgrid, d2ll_kde),col="blue") +
    geom_line(aes(xgrid, d2ll_numd), col = "red", linetype = "dashed", size = 2) +
    geom_line(aes(xgrid, d2ll_exact), col="red")
  g2bis = ggplot() +
    geom_vline(xintercept = observations[,t,drop=FALSE]) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_line(aes(xgrid, abs(d2ll_kde-d2ll_exact)),col="blue") +
    scale_y_log10()

  grid.arrange(g0,g0bis,g1,g1bis,g2,g2bis,nrow=3)

  # compute derivatives at the observation
  # h = hlscv(t(Yt_sim),tol/2,deriv.order = 2)
  # h = hpi(t(Yt_sim),deriv.order = 2)
  # h = hscv(t(Yt_sim))
  h = tol/4

  l_kde_t = kdde(t(Yt_sim),deriv.order = 0,eval.points = observations[,t], h = h)$estimate
  d1_kde_t = kdde(t(Yt_sim),deriv.order = 1,eval.points = observations[,t], h = h)$estimate
  d2_kde_t = kdde(t(Yt_sim),deriv.order = 2,eval.points = observations[,t], h = h)$estimate
  # compute log-derivative at the observation
  ll_kde_t = log(l_kde_t)
  d1ll_kde_t = d1_kde_t/l_kde_t
  d2ll_kde_t = d2_kde_t/l_kde_t - (d1ll_kde_t)^2
  # compute h score via KDE
  hscore_kde[t] = (2*d2ll_kde_t + (d1ll_kde_t)^2)
  # compute h score exact
  s = sum(observations[,1:t]^2)
  hscore_exact[t] = ((nu0+t)/((nu0+s)^2))*((nu0+t+4)*observations[,t]^2-2*(nu0+s))
  # compute h score exact via numDeriv
  d1ll_t = grad(ll_exact_fn,observations[,t])
  d2ll_t = hessian(ll_exact_fn,observations[,t])
  hscore_numderiv[t] = (2*d2ll_t + (d1ll_t)^2)
}

ggplot() +
  geom_line(aes(1:nobservations, hscore_kde), size = 1, color = "blue") +
  geom_line(aes(1:nobservations, hscore_numderiv),color="red") +
  geom_line(aes(1:nobservations, hscore_exact),color="red",size=2,linetype=2)


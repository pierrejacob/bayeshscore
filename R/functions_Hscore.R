#----------------------------------------------------------------------------------------#
#---------------------------- CONTINUOUS OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
# smc2 version when predictive density is intractable
hincrement_continuous_smc2_ = function(t,model,observationt,thetas,Wtheta,PFs) {
  Ntheta = ncol(thetas)
  hincrement = 0
  d1log = list()
  d2log = list()
  for (k in 1:model$dimY) {
    Ek_theta = vector("numeric",Ntheta)
    Fk_theta = vector("numeric",Ntheta)
    for (m in 1:Ntheta){
      Ek_theta[m] = 0
      Fk_theta[m] = 0
      if (Wtheta[m] == 0) {
        next
      }
      PF = PFs[[m]]
      X = PF$X
      WX = PF$xnormW
      if (is.null(dim(X))){
        X = matrix(X,nrow = model$dimX)
      }
      # compute the derivatives for all k and all X at once during the first loop
      if (k==1){
        derivatives = model$derivativelogdobs(observationt,X,t,thetas[,m])
        d1log[[m]] = derivatives$jacobian
        d2log[[m]] = derivatives$hessiandiag
      }
      # we automatically discard particles with weight 0 to avoid artificial apparition of NaN
      # (which could happen when performing 0*Inf)
      nonzeroWX = (WX>0)
      Ek_theta[m] = sum(WX[nonzeroWX]*(d2log[[m]][nonzeroWX,k] + (d1log[[m]][nonzeroWX,k])^2))
      Fk_theta[m] = sum(WX[nonzeroWX]*d1log[[m]][nonzeroWX,k])
    }
    # we automatically discard particles with weight 0 to avoid artificial apparition of NaN
    # (which could happen when performing 0*Inf)
    nonzeroWtheta = (Wtheta>0)
    Ek = sum(Wtheta[nonzeroWtheta]*Ek_theta[nonzeroWtheta])
    Fk = (sum(Wtheta[nonzeroWtheta]*Fk_theta[nonzeroWtheta]))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}
#-------------------------------------------------------------------------------------------
# smc version when predictive density is available
hincrement_continuous_smc_ = function(t,model,observations,thetas,Wtheta,byproducts) {
  Ntheta = ncol(thetas)
  hincrement = 0
  d1log = list()
  d2log = list()
  for (k in 1:model$dimY) {
    Ek_theta = array(0,dim = Ntheta)
    Fk_theta = array(0,dim = Ntheta)
    for (m in 1:Ntheta){
      if (Wtheta[m] == 0) {
        next
      }
      # compute the derivatives for all k at once during the first loop
      if (k==1){
        if (!is.null(byproducts)){
          derivatives = model$derivativelogdpredictive(observations,t,thetas[,m],byproducts[[m]])
        } else {
          derivatives = model$derivativelogdpredictive(observations,t,thetas[,m])
        }
        d1log[[m]] = derivatives$jacobian
        d2log[[m]] = derivatives$hessiandiag
      }
      Ek_theta[m] = d2log[[m]][k] + (d1log[[m]][k])^2
      Fk_theta[m] = d1log[[m]][k]
    }
    # we automatically discard particles with weight 0 to avoid artificial apparition of NaN
    # (which could happen when performing 0*Inf)
    nonzeroWtheta = (Wtheta>0)
    Ek = sum(Wtheta[nonzeroWtheta]*Ek_theta[nonzeroWtheta])
    Fk = (sum(Wtheta[nonzeroWtheta]*Fk_theta[nonzeroWtheta]))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}
#----------------------------------------------------------------------------------------#
#------------------------------ DISCRETE OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
# This function computes the approximation pt_hat(y) for a given set of particles
phat_smc2 = function(t,model,y,thetas,thetanormw,Ntheta,Xpred,xprednormw) {
  py = 0
  for (m in 1:Ntheta){
    if (thetanormw[m]==0){
      next
    }
    py = py + thetanormw[m]*sum(xprednormw[,m]*model$dobs(y,Xpred[,,m],t,thetas[,m],log = FALSE))
  }
  return (py)
}
# This function computes the partial score term Hdk (forward difference version)
# a,b are vectors of componentwise lower and upper bounds of the observations
Hdk_forward = function(k,a,b,d,y,py_minusek,py,py_plusek) {
  if (y[k]==b[k]) {
    return (-2*(py-py_minusek)/py_minusek)
  }
  else {
    if (y[k]==a[k]) {
      return (2*(py_plusek-py)/py + ((py_plusek-py)/py)^2)
    }
    else {
      return (2*((py_plusek-py)/py-(py-py_minusek)/py_minusek) + ((py_plusek-py)/py)^2)
    }
  }
}
# This function computes the partial score term Hdk (central difference version)
# a,b are vectors of componentwise lower and upper bounds of the observations
Hdk_central = function(k,a,b,d,y,py_minus2ek,py_minusek,py,py_plusek,py_plus2ek) {
  if (y[k]==b[k]) {
    return (-(py-py_minus2ek)/(2*py_minusek))
  }
  if (y[k]==b[k]-1) {
    return (-(py-py_minus2ek)/(2*py_minusek) + ((py_plusek-py_minusek)/(2*py))^2)
  }
  if (y[k]==a[k]) {
    return ((py_plus2ek-py)/(2*py_plusek))
  }
  if (y[k]==a[k]+1) {
    return ((py_plus2ek-py)/(2*py_plusek) + ((py_plusek-py_minusek)/(2*py))^2)
  }
  if (y[k]!=b[k] && y[k]!=b[k]-1 && y[k]!=a[k] && y[k]!=a[k]+1){
    return ((py_plus2ek-py)/(2*py_plusek) - (py-py_minus2ek)/(2*py_minusek) + ((py_plusek-py_minusek)/(2*py))^2)
  }
}
# This function computes the partial score term Hdk and is a wrapper of Hdk_forward and Hdk_central
# a,b are vectors of componentwise lower and upper bounds of the observations
Hdk = function(k,a,b,d,y,py_minusek,py,py_plusek,diff_type = "central",py_minus2ek = NULL,py_plus2ek = NULL) {
  if (diff_type == "central"){
    return (Hdk_central(k,a,b,d,y,py_minus2ek,py_minusek,py,py_plusek,py_plus2ek))
  }
  if (diff_type == "forward"){
    return (Hdk_forward(k,a,b,d,y,py_minusek,py,py_plusek))
  }
}
# This function computes the partial score term Hd
# a,b are vectors of componentwise lower and upper bounds of the observations
# d is the dimension of y
Hd_smc2 = function(t,model,yt,thetas,thetanormw,Xpred,xprednormw, diff_type = "central") {
  Ntheta = ncol(thetas)
  a = model$lower
  b = model$upper
  d = model$dimY
  result = 0
  for (k in 1:d) {
    ek = rep(0,d)
    ek[k] = 1
    py = phat_smc2(t,model,yt,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    py_minusek = phat_smc2(t,model,yt-ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    py_plusek = phat_smc2(t,model,yt+ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    if (diff_type == "central") {
      py_minus2ek = phat_smc2(t,model,yt-2*ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
      py_plus2ek = phat_smc2(t,model,yt+2*ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    } else {
      py_minus2ek = NULL
      py_plus2ek = NULL
    }
    result = result + Hdk(k,a,b,d,yt,py_minusek,py,py_plusek, diff_type, py_minus2ek, py_plus2ek)
  }
  return (result)
}
#-------------------------------------------------------------------------------------------
# smc version when predictive density is available
# This function computes the approximation qt_hat(y) for a given set of particles
phat_smc = function(t,model,observations,thetas,thetanormw,Ntheta,byproducts) {
  py = 0
  if (!is.null(byproducts)){
    for (m in 1:Ntheta){
      if (thetanormw[m]==0){
        next
      }
      py = py + thetanormw[m]*model$dpredictive(observations,t,thetas[,m],byproducts[[m]],log = FALSE)
    }
  } else {
    for (m in 1:Ntheta){
      if (thetanormw[m]==0){
        next
      }
      py = py + thetanormw[m]*model$dpredictive(observations,t,thetas[,m],log = FALSE)
    }
  }
  return (py)
}
# This function computes the partial score term Hd
# a,b are vectors of componentwise lower and upper bounds of the observations
# d is the dimension of y
Hd_smc = function(t,model,observations,thetas,thetanormw,byproducts, diff_type = "central") {
  Ntheta = ncol(thetas)
  a = model$lower
  b = model$upper
  d = model$dimY
  result = 0
  for (k in 1:d) {
    ek = rep(0,d)
    ek[k] = 1
    obsminusek = observations
    obsplusek = observations
    obsminusek[,t] = obsminusek[,t]-ek
    obsplusek[,t] = obsplusek[,t]+ek
    if (diff_type == "central") {
      obsminus2ek = observations
      obsplus2ek = observations
      obsminus2ek[,t] = obsminus2ek[,t]-2*ek
      obsplus2ek[,t] = obsplus2ek[,t]+2*ek
    } else {
      obsminus2ek = NULL
      obsplus2ek = NULL
    }
    py = phat_smc(t,model,observations,thetas,thetanormw,Ntheta,byproducts)
    py_minusek = phat_smc(t,model,obsminusek,thetas,thetanormw,Ntheta,byproducts)
    py_plusek = phat_smc(t,model,obsplusek,thetas,thetanormw,Ntheta,byproducts)
    if (diff_type == "central") {
      py_minus2ek = phat_smc(t,model,obsminus2ek,thetas,thetanormw,Ntheta,byproducts)
      py_plus2ek = phat_smc(t,model,obsplus2ek,thetas,thetanormw,Ntheta,byproducts)
    } else {
      py_minus2ek = NULL
      py_plus2ek = NULL
    }
    result = result + Hdk(k,a,b,d,observations[,t],py_minusek,py,py_plusek, diff_type, py_minus2ek, py_plus2ek)
  }
  return (result)
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Wrappers: the following function performs all the preprocessing needed to compute
# the Hyvarinen score (e.g. generate x-particles from the one-step predictive distribution,
# generate more theta-particles for variance reduction if needed, etc ...)
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Discrete case, SMC2: wrapper of Hd_smc2
hincrement_discrete_smc2 = function(thetas, normw, PFs, t, observations, model,
                                    logtargetdensities, algorithmic_parameters) {
  discrete_diff_type = algorithmic_parameters$discrete_diff_type
  # If reduce_variance: temporarily increase the number of particles before computing
  # the Hyvarinen score.
  if (algorithmic_parameters$reduce_variance) {
    # generate additional theta-particles (and their associated particle filters)
    larger_pool =  get_additional_particles_smc2(algorithmic_parameters$Nc, thetas, normw, PFs, t,
                                                 observations, model, logtargetdensities,
                                                 algorithmic_parameters, algorithmic_parameters$Ncx)
    thetas_pool = larger_pool$thetas
    PFs_pool = larger_pool$PFs
    normw_pool = rep(1/ncol(thetas_pool), ncol(thetas_pool))
    # Get the new number of particles
    Ntheta_pool = ncol(thetas_pool)
    Nx_pool = PFs_pool[[1]]$Nx
    # generate x-particles targeting the one-step predictive distribution which are already
    # contained in the PFs after initialization
    if (t==1) {
      Xpred = array(NA,dim = c(model$dimX,Nx_pool,Ntheta_pool)) # (need to reconstruct since size Nx might change)
      XnormW_previous = matrix(1/Nx_pool, nrow = Nx_pool, ncol = Ntheta_pool)
      for (itheta in 1:Ntheta_pool){
        Xpred[,,itheta] = PFs_pool[[itheta]]$X
      }
    }
    # Construct particles targeting the one-step-ahead predictive
    if (t > 1){
      Xpred = array(NA,dim = c(model$dimX,Nx_pool,Ntheta_pool))
      XnormW_previous = matrix(NA, nrow = Nx_pool, ncol = Ntheta_pool)
      # At this point in the code, PFs[[itheta]]$X contains the x-particles from the filtering
      # distribution at time (t-1)
      for (itheta in 1:Ntheta_pool){
        if (is.null(dim(PFs_pool[[itheta]]$X))){
          XnormW_previous[,itheta] = PFs_pool[[itheta]]$xnormW
          Xpred[,,itheta] = model$rtransition(matrix(PFs_pool[[itheta]]$X,ncol = Nx_pool), t, thetas_pool[,itheta])
        }
        else{
          XnormW_previous[,itheta] = PFs_pool[[itheta]]$xnormW
          Xpred[,,itheta] = model$rtransition(PFs_pool[[itheta]]$X, t, thetas_pool[,itheta])
        }
      }
    }
    # compute incremental H score (with theta from time t-1, see formula in the paper)
    return (Hd_smc2(t,model,observations[,t],thetas_pool,normw_pool,Xpred,XnormW_previous, discrete_diff_type))
  } else {
    Ntheta = ncol(thetas)
    Nx = PFs[[1]]$Nx
    # Construct particles targeting the one-step-ahead predictive (need to reconstruct since size Nx might change)
    if (t==1) {
      Xpred = array(NA,dim = c(model$dimX, Nx, Ntheta))
      XnormW_previous = matrix(1/Nx, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for X at previous step
      for (itheta in 1:Ntheta){
        Xpred[,,itheta] = PFs[[itheta]]$X
      }
    }
    if (t > 1){
      Nx = PFs[[1]]$Nx
      Xpred = array(NA,dim = c(model$dimX,Nx,Ntheta)) # (need to reconstruct since size Nx might change)
      XnormW_previous = matrix(NA, nrow = Nx, ncol = Ntheta)
      # At this point in the code, PFs[[itheta]]$X contains the x-particles from the filtering
      # distribution at time (t-1)
      for (itheta in 1:Ntheta){
        if (is.null(dim(PFs[[itheta]]$X))){
          XnormW_previous[,itheta] = PFs[[itheta]]$xnormW
          Xpred[,,itheta] = model$rtransition(matrix(PFs[[itheta]]$X,ncol = Nx), t, thetas[,itheta])
        }
        else{
          XnormW_previous[,itheta] = PFs[[itheta]]$xnormW
          Xpred[,,itheta] = model$rtransition(PFs[[itheta]]$X, t, thetas[,itheta])
        }
      }
    }
    # compute incremental H score (with theta from time t-1, see formula in the paper)
    return (Hd_smc2(t,model,observations[,t],thetas,normw,Xpred,XnormW_previous, discrete_diff_type))
  }
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Discrete case, SMC: wrapper of Hd_smc
hincrement_discrete_smc = function(thetas, normw, byproducts, t, observations, model,
                                   logtargetdensities, algorithmic_parameters) {
  discrete_diff_type = algorithmic_parameters$discrete_diff_type
  # compute incremental H score (with theta from time t-1, see formula in the paper)
  if (algorithmic_parameters$reduce_variance) {
    # generate additional particles
    larger_pool =  get_additional_particles_smc(algorithmic_parameters$Nc, thetas, normw,
                                                byproducts, t, observations, model,
                                                logtargetdensities, algorithmic_parameters)
    thetas_pool = larger_pool$thetas
    byproducts_pool = larger_pool$byproducts
    normw_pool = rep(1/ncol(thetas_pool), ncol(thetas_pool))
    # compute Hyvarinen score with more particles
    return (Hd_smc(t,model,observations,thetas_pool,normw_pool,byproducts_pool,diff_type = discrete_diff_type))
  } else {
    return (Hd_smc(t,model,observations,thetas,normw,byproducts,diff_type = discrete_diff_type))
  }
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Continuous case, SMC2: wrapper of hincrement_continuous_smc2_
hincrement_continuous_smc2 = function(thetas, normw, PFs, t, observations, model,
                                      logtargetdensities, algorithmic_parameters) {
  # If reduce_variance: temporarily increase the number of particles before computing
  # the Hyvarinen score.
  if (algorithmic_parameters$reduce_variance) {
    # generate additional theta-particles (and their associated particle filters)
    larger_pool =  get_additional_particles_smc2(algorithmic_parameters$Nc, thetas, normw, PFs, t,
                                                 observations, model, logtargetdensities,
                                                 algorithmic_parameters, algorithmic_parameters$Ncx)
    thetas_pool = larger_pool$thetas
    PFs_pool = larger_pool$PFs
    normw_pool = rep(1/ncol(thetas_pool), ncol(thetas_pool))
    return (hincrement_continuous_smc2_(t,model,observations[,t,drop=FALSE],thetas_pool,normw_pool,PFs_pool))
  } else {
    return (hincrement_continuous_smc2_(t,model,observations[,t,drop=FALSE],thetas,normw,PFs))
  }
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Continuous case, SMC: wrapper of hincrement_continuous_smc_
hincrement_continuous_smc = function(thetas, normw, byproducts, t, observations, model,
                                     logtargetdensities, algorithmic_parameters) {
  # If reduce_variance: temporarily increase the number of particles before computing
  # the Hyvarinen score.
  if (algorithmic_parameters$reduce_variance) {
    # generate additional particles
    larger_pool =  get_additional_particles_smc(algorithmic_parameters$Nc, thetas, normw,
                                                byproducts, t, observations, model,
                                                logtargetdensities, algorithmic_parameters)
    thetas_pool = larger_pool$thetas
    byproducts_pool = larger_pool$byproducts
    normw_pool = rep(1/ncol(thetas_pool), ncol(thetas_pool))
    # compute Hyvarinen score with more particles
    return (hincrement_continuous_smc_(t, model, observations,thetas_pool,normw_pool,byproducts_pool))
  } else {
    return (hincrement_continuous_smc_(t, model, observations,thetas,normw,byproducts))
  }
}

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Additional implementations using Kernel density estimators
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# smc version when predictive density is available, using density estimators
# via local regression. WARNING: only implemented for univariate observations.
hincrementContinuous_smc_dde = function(t,model,observations,thetas,normw,byproducts,
                                        logtargetdensities, algorithmic_parameters){
  if (t > algorithmic_parameters$dde_options$nb_steps){
    return (NA)
  }
  if (is.null(model$rpredictive)) {
    print("Can't perform density estimation: no sampler 'rpredictive' provided")
    return (NA)
  }
  if (model$dimY!=1) {
    print("WARNING: density estimation has only been implemented for univariate observations")
    return (NA)
  } else {
    Ny = algorithmic_parameters$dde_options$Ny
    sigma20 = algorithmic_parameters$dde_options$sigma2_order0
    sigma21 = algorithmic_parameters$dde_options$sigma2_order1
    sigma22 = algorithmic_parameters$dde_options$sigma2_order2
    if (t == 1){
      thetas_sim = model$rprior(Ny)
      Yt_sim = apply(thetas_sim,MARGIN = 2, function(theta)model$rpredictive(1,1,theta,NULL))
    }
    if (t >=2){
      # generate additional particles until we have at least Ny particles
      larger_pool =  get_additional_particles_smc(Ny, thetas, normw, byproducts, t, observations, model,
                                                  logtargetdensities, algorithmic_parameters)
      # Note that the particles thetas returned by get_additional_particles_smc are equally weighted
      thetas_pool = larger_pool$thetas
      # Draw one Yt for each particle theta
      Yt_sim = apply(thetas_pool,MARGIN = 2, function(theta)model$rpredictive(1,t,theta,observations[,1:(t-1),drop=FALSE]))
    }
    d0_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma20,order = 0)
    d1_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma21,order = 1)
    d2_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma22,order = 2)
    # compute h score via local estimation
    return ((2*d2_est/d0_est)-((d1_est/d0_est)^2))
  }
}
# smc version when predictive density is available, using density estimators
# via local regression. WARNING: only implemented for univariate observations.
hincrementContinuous_smc2_dde = function(t,model,observations,thetas,normw,PFs,
                                         logtargetdensities, algorithmic_parameters){
  if (t > algorithmic_parameters$dde_options$nb_steps){
    return (NA)
  }
  if (is.null(model$robs)) {
    print("Can't perform density estimation: no sampler 'robs' provided")
    return (NA)
  }
  if (model$dimY!=1) {
    print("WARNING: density estimation has only been implemented for univariate observations")
    return (NA)
  } else {
    Ny = algorithmic_parameters$dde_options$Ny
    sigma20 = algorithmic_parameters$dde_options$sigma2_order0
    sigma21 = algorithmic_parameters$dde_options$sigma2_order1
    sigma22 = algorithmic_parameters$dde_options$sigma2_order2
    if (t == 1){
      thetas_pool = model$rprior(Ny)
      Xt_sim = apply(thetas_pool,MARGIN = 2, function(theta)model$rinitial(theta,1))
      Ntheta_larger = ncol(thetas_pool)
      Yt_sim = rep(NA,Ntheta_larger)
      for (i in 1:Ntheta_larger){
        Yt_sim[i] = model$robs(Xt_sim[,i,drop=FALSE],1,thetas_pool[,i,drop=FALSE])
      }
    }
    if (t >=2){
      # generate additional particles until we have at least Ny particles
      larger_pool =  get_additional_particles_smc2(Ny, thetas, normw, PFs, t, observations, model,
                                                   logtargetdensities, algorithmic_parameters, NULL)
      # Note that the particles thetas returned by get_additional_particles_smc are equally weighted
      thetas_pool = larger_pool$thetas
      PFs_pool = larger_pool$PFs
      Ntheta_larger = ncol(thetas_pool)
      Yt_sim = rep(NA,Ntheta_larger)
      for (i in 1:Ntheta_larger){
        Xprev = PFs_pool[[i]]$X[,sample(1:PFs_pool[[i]]$Nx,size = 1,prob = PFs_pool[[i]]$xnormW),drop=FALSE]
        Xt_sim = model$rtransition(Xprev,t,thetas_pool[,i])
        Yt_sim[i] = model$robs(Xt_sim,t,thetas_pool[,i])
      }
    }
    d0_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma20,order = 0)
    d1_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma21,order = 1)
    d2_est = get_derivative_RBFlocal(Yt_sim,observations[,t],sigma22,order = 2)
    # compute h score via local estimation
    return ((2*d2_est/d0_est)-((d1_est/d0_est)^2))
  }
}





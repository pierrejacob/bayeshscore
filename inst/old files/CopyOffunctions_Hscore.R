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
# This function computes the partial score term Hdk
# a,b are vectors of componentwise lower and upper bounds of the observations
Hdk = function(k,a,b,d,y,py_minusek,py,py_plusek) {
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
# This function computes the partial score term Hd
# a,b are vectors of componentwise lower and upper bounds of the observations
# d is the dimension of y
Hd_smc2 = function(t,model,yt,thetas,thetanormw,Xpred,xprednormw) {
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
    result = result + Hdk(k,a,b,d,yt,py_minusek,py,py_plusek)
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
Hd_smc = function(t,model,observations,thetas,thetanormw,byproducts) {
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
    py = phat_smc(t,model,observations,thetas,thetanormw,Ntheta,byproducts)
    py_minusek = phat_smc(t,model,obsminusek,thetas,thetanormw,Ntheta,byproducts)
    py_plusek = phat_smc(t,model,obsplusek,thetas,thetanormw,Ntheta,byproducts)
    result = result + Hdk(k,a,b,d,observations[,t],py_minusek,py,py_plusek)
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
  # If reduce_variance: temporarily increase the number of particles before computing
  # the Hyvarinen score.
  if (algorithmic_parameters$reduce_variance) {
    # generate additional theta-particles (and their associated particle filters)
    larger_pool =  get_additional_particles_smc2(thetas, normw, PFs, t, observations, model,
                                                 logtargetdensities, algorithmic_parameters)
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
    return (Hd_smc2(t,model,observations[,t],thetas_pool,normw_pool,Xpred,XnormW_previous))
  } else {
    Ntheta = ncol(thetas)
    Nx = PFs_pool[[1]]$Nx
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
    return (Hd_smc2(t,model,observations[,t],thetas,normw,Xpred,XnormW_previous))
  }
}
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Discrete case, SMC: wrapper of Hd_smc
hincrement_discrete_smc = function(thetas, normw, byproducts, t, observations, model,
                                   logtargetdensities, algorithmic_parameters) {
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
    return (Hd_smc(t,model,observations,thetas_pool,normw_pool,byproducts_pool))
  } else {
    return (Hd_smc(t,model,observations,thetas,normw,byproducts))
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
    larger_pool =  get_additional_particles_smc2(thetas, normw, PFs, t, observations, model,
                                                 logtargetdensities, algorithmic_parameters)
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
# smc version when predictive density is available, using Kernel density estimators
# from the ks package
hincrementContinuous_smc_kde = function(t,model,observations,thetas,normw,byproducts,
                                        logtargetdensities, algorithmic_parameters){
  Ny = algorithmic_parameters$parameters_kde$Ny
  # generate additional particles until we have at least Ny particles
  larger_pool =  get_additional_particles_smc(Ny, thetas, normw, byproducts, t, observations, model,
                                              logtargetdensities, algorithmic_parameters)
  # Note that the particles thetas returned by get_additional_particles_smc are equally weighted
  thetas_pool = larger_pool$thetas
  Ntheta = ncol(thetas_pool)
  # Draw one Yt for each particle theta
  Yt_sim = matrix(NA, nrow = model$dimY, ncol = Ntheta)
  if (t == 1){
    for (i in 1:Ntheta){
      Yt_sim[,i] = model$rpredictive(1,t,thetas_pool[,i],NULL)
    }
  } else if (t>=2){
    for (i in 1:Ntheta){
      Yt_sim[,i] = model$rpredictive(1,t,thetas_pool[,i],observations[,1:(t-1),drop=FALSE])
    }
  }


  xgrid = seq(-5,5,0.1)
  # print(Yt_sim)
  g = ggplot() +  geom_histogram(aes(x=c(Yt_sim), y=..density..), alpha = 0.6) +
    geom_line(aes(xgrid, kdde(t(Yt_sim), eval.points = xgrid)$estimate)) +
    geom_vline(xintercept = observations[,t,drop=FALSE])
  plot(g)



  pred_kde = kdde(t(Yt_sim),deriv.order = 0,eval.points = observations[,t,drop=FALSE])$estimate
  d1pred_kde = kdde(t(Yt_sim),deriv.order = 1,eval.points = observations[,t,drop=FALSE])$estimate
  d2pred_kde = kdde(t(Yt_sim),deriv.order = 2,eval.points = observations[,t,drop=FALSE])$estimate
  # return the h score increment
  return (2*d2pred_kde/pred_kde - (d1pred_kde/pred_kde)^2)
  # return (2*sum(diag(d2pred_kde))/pred_kde - (d1pred_kde%*%t(d1pred_kde)/pred_kde))
}





#----------------------------------------------------------------------------------------#
#---------------------------- CONTINUOUS OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
# smc2 version when predictive density is intractable
hincrementContinuous_smc2 = function(t,model,observationt,thetas,Wtheta,PFs,Ntheta) {
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
        derivatives = model$derivativelogdobs(observationt,X,t,thetas[,m],model$dimY)
        d1log[[m]] = derivatives$jacobian
        d2log[[m]] = derivatives$hessiandiag
      }
      Ek_theta[m] = sum(WX*(d2log[[m]][,k] + (d1log[[m]][,k])^2))
      Fk_theta[m] = sum(WX*d1log[[m]][,k])
    }
    Ek = sum(Wtheta*Ek_theta)
    Fk = (sum(Wtheta*Fk_theta))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}
#-------------------------------------------------------------------------------------------
# smc version when predictive density is available
hincrementContinuous_smc = function(t,model,observations,thetas,Wtheta,byproducts,Ntheta) {
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
          derivatives = model$derivativelogdpredictive(observations,t,thetas[,m],byproducts[[m]],model$dimY)
        } else {
          derivatives = model$derivativelogdpredictive(observations,t,thetas[,m],model$dimY)
        }
        d1log[[m]] = derivatives$jacobian
        d2log[[m]] = derivatives$hessiandiag
      }
      Ek_theta[m] = d2log[[m]][k] + (d1log[[m]][k])^2
      Fk_theta[m] = d1log[[m]][k]
    }
    Ek = sum(Wtheta*Ek_theta)
    Fk = (sum(Wtheta*Fk_theta))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}

#----------------------------------------------------------------------------------------#
#------------------------------ DISCRETE OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
# This function computes the approximation pt_hat(y) for a given set of particles
phat = function(t,model,y,thetas,thetanormw,Ntheta,Xpred,xprednormw) {
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
Hd = function(t,model,yt,thetas,thetanormw,Ntheta,Xpred,xprednormw) {
  a = model$lower
  b = model$upper
  d = model$dimY
  result = 0
  for (k in 1:d) {
    ek = rep(0,d)
    ek[k] = 1
    py = phat(t,model,yt,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    py_minusek = phat(t,model,yt-ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
    py_plusek = phat(t,model,yt+ek,thetas,thetanormw,Ntheta,Xpred,xprednormw)
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
Hd_smc = function(t,model,observations,thetas,thetanormw,Ntheta,byproducts) {
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

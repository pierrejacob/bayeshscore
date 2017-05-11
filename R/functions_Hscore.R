#----------------------------------------------------------------------------------------#
#---------------------------- CONTINUOUS OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
# smc2 version when predictive density is intractable
hincrementContinuous_smc2 = function(t,model,observationt,thetas,Wtheta,PFs,Ntheta) {
  d = model$dimY
  hincrement = 0
  for (k in 1:d) {
    Ek_theta = vector("numeric",Ntheta)
    Fk_theta = vector("numeric",Ntheta)
    for (m in 1:Ntheta){
      PF = PFs[[m]]
      X = PF$X
      if (is.null(dim(X))){
        X = matrix(X,nrow = model$dimX)
      }
      WX = PF$xnormW
      Ek_theta[m] = 0
      Fk_theta[m] = 0
      if (Wtheta[m] == 0) {
        next
      }
      dlogobs_k = model$derivativelogdobs(observationt,X,t,thetas[,m],k)
      Ek_theta[m] = sum(WX*(dlogobs_k$d2log + (dlogobs_k$d1log)^2))
      Fk_theta[m] = sum(WX*dlogobs_k$d1log)
    }
    Ek = sum(Wtheta*Ek_theta)
    Fk = (sum(Wtheta*Fk_theta))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}
# smc version when predictive density is available
hincrementContinuous_smc = function(t,model,observations,thetas,Wtheta,Ntheta,byproducts) {
  d = model$dimY
  hincrement = 0
  gradients = list()
  hessians = list()
  for (k in 1:d) {
    Ek_theta = array(0,dim = Ntheta)
    Fk_theta = array(0,dim = Ntheta)
    for (m in 1:Ntheta){
      if (Wtheta[m] == 0) {
        next
      }
      # compute the derivatives for all k at once during the first loop
      if (k==1){
        derivatives = model$derivativelogdpredictive(observations,t,thetas[,m],byproducts[[m]])
        gradients[[m]] = derivatives$grad
        hessians[[m]] = derivatives$hessian
      }
      Ek_theta[m] = hessians[[m]][k,k] + (gradients[[m]][k])^2
      Fk_theta[m] = gradients[[m]][k]
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
# This function computes the approximation qt_hat(y) for a given set of particles
phat = function(t,model,y,thetas,thetanormw,Xpred,xprednormw,Ntheta) {
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
Hd = function(t,model,yt,thetas,thetanormw,Xpred,xprednormw,Ntheta) {
  a = model$lower
  b = model$upper
  d = model$dimY
  result = 0
  for (k in 1:d) {
    ek = rep(0,d)
    ek[k] = 1
    py = phat(t,model,yt,thetas,thetanormw,Xpred,xprednormw,Ntheta)
    py_minusek = phat(t,model,yt-ek,thetas,thetanormw,Xpred,xprednormw,Ntheta)
    py_plusek = phat(t,model,yt+ek,thetas,thetanormw,Xpred,xprednormw,Ntheta)
    result = result + Hdk(k,a,b,d,yt,py_minusek,py,py_plusek)
  }
  return (result)
}

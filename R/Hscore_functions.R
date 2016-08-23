hincrementContinuous = function(t,model,observationt,thetas,Wtheta,X,WX,Ntheta,Nx) {
  d = model$dimY
  hincrement = 0
  for (k in 1:d) {
    Ek_theta = vector("numeric",Ntheta)
    Fk_theta = vector("numeric",Ntheta)
    for (m in 1:Ntheta){
      Ek_theta[m] = 0
      Fk_theta[m] = 0
      if (is.null(dim(X[,,m]))){
        dlogobs_k = model$derivativelogdobs(observationt,matrix(X[,,m],ncol = model$dimX),t,thetas[m,],k)
      }
      else{
        dlogobs_k = model$derivativelogdobs(observationt,X[,,m],t,thetas[m,],k)
      }
      Ek_theta[m] = sum(WX[,m]*(dlogobs_k$d2log + (dlogobs_k$d1log)^2))
      Fk_theta[m] = sum(WX[,m]*dlogobs_k$d1log)
    }
    Ek = sum(Wtheta*Ek_theta)
    Fk = (sum(Wtheta*Fk_theta))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}

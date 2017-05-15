#----------------------------------------------------------------------------------------#
#---------------------------- CONTINUOUS OBSERVATIONS -----------------------------------#
#----------------------------------------------------------------------------------------#
hincrementContinuous_no_tempering = function(t,model,observationt,thetas,Wtheta,X,WX,Ntheta,Nx) {
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
      if (k==1){
        derivatives = model$derivativelogdobs(observationt,matrix(X[,,m], nrow = model$dimX),t,thetas[,m],model$dimY)
        d1log[[m]] = derivatives$jacobian
        d2log[[m]] = derivatives$hessiandiag
      }
      Ek_theta[m] = sum(WX[,m]*(d2log[[m]][k] + (d1log[[m]][,k])^2))
      Fk_theta[m] = sum(WX[,m]*d1log[[m]][,k])
    }
    Ek = sum(Wtheta*Ek_theta)
    Fk = (sum(Wtheta*Fk_theta))^2
    hincrement = hincrement + 2*Ek - Fk
  }
  return (hincrement)
}
# #----------------------------------------------------------------------------------------#
# #------------------------------ DISCRETE OBSERVATIONS -----------------------------------#
# #----------------------------------------------------------------------------------------#
# # This function computes the approximation qt_hat(y) for a given set of particles
# qhat = function(t,model,y,thetas,thetanormw,Xpred,xprednormw,Ntheta,Nx) {
#   qy = 0
#   for (m in 1:Ntheta){
#     if (thetanormw[m]==0){
#       next
#     }
#     qy = qy + thetanormw[m]*sum(xprednormw[,m]*model$dobs(y,Xpred[,,m],t,thetas[,m],log = FALSE))
#   }
#   return (qy)
# }
# # This function computes the partial score term SBk
# # a,b are vectors of componentwise lower and upper bounds of the observations
# SBk = function(k,a,b,d,y,qy_minusek,qy,qy_plusek) {
#   if (y[k]==b[k]) {
#     return (-2*(qy-qy_minusek)/qy_minusek)
#   }
#   else {
#     if (y[k]==a[k]) {
#       return (2*(qy_plusek-qy)/qy + ((qy_plusek-qy)/qy)^2)
#     }
#     else {
#       return (2*((qy_plusek-qy)/qy-(qy-qy_minusek)/qy_minusek) + ((qy_plusek-qy)/qy)^2)
#     }
#   }
# }
# # This function computes the partial score term SHd
# # a,b are vectors of componentwise lower and upper bounds of the observations
# # d is the dimension of y
# SHd = function(t,model,yt,thetas,thetanormw,Xpred,xprednormw,Ntheta,Nx) {
#   a = model$lower
#   b = model$upper
#   d = model$dimY
#   result = 0
#   for (k in 1:d) {
#     ek = rep(0,d)
#     ek[k] = 1
#     qy = qhat(t,model,yt,thetas,thetanormw,Xpred,xprednormw,Ntheta,Nx)
#     qy_minusek = qhat(t,model,yt-ek,thetas,thetanormw,Xpred,xprednormw,Ntheta,Nx)
#     qy_plusek = qhat(t,model,yt+ek,thetas,thetanormw,Xpred,xprednormw,Ntheta,Nx)
#     result = result + SBk(k,a,b,d,yt,qy_minusek,qy,qy_plusek)
#   }
#   return (result)
# }

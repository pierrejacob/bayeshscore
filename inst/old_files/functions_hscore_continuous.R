#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains functions to compute the prequential Hyvarinen score in CONTINUOUS case
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# WARNING: the current code only works for univariate states !!!
library(Matrix) #contains nearPD (nearest positive definite matrix)
library(mvtnorm) #contains Multivariate Normal
source("functions_resampling.R")
source("functions_proposal.R")
source("functions_smc.R")
source("functions_finite_differences.R")

# This function computes the ESS given a set of normalized weights
getESS = function(normalized_weights) {
  return (1/sum(normalized_weights^2))
}

# This function computes the prequential Hyvarinen score increment
hincrementContinuous = function(t,model,yt,theta,Wtheta,X,WX,Ntheta,Nx) {
  d = length(yt)
  hincrement = 0
  for (k in 1:d) {
    term1theta = vector("numeric",Ntheta)
    term2theta = vector("numeric",Ntheta)
    for (m in 1:Ntheta){
      term1theta[m] = 0
      term2theta[m] = 0
      for (n in 1:Nx) {
        g = function(y) model$dobs(y,X[n,m],t,theta[m,])
        d1k = d1log(g,yt,k)
        d2k = d2log(g,yt,k)
        term1theta[m] = term1theta[m] + WX[n,m]*(d2k + (d1k^2))
        term2theta[m] = WX[n,m]*d1k
      }
    }
    hincrement = hincrement + 2*sum(Wtheta*term1theta) - (sum(Wtheta*term2theta))^2
  }
  return (hincrement)
}

# This function compute the total prequential Hyvarinen Score
hyvarinenContinuous = function(y, T, model,  Ntheta, Nx, smc = SMC, resamplingSMC = multinomial, 
                               proposalSMC = transition, progress = TRUE) {
  if (progress) {
    progbar = txtProgressBar(min = 0,max = Ntheta*T,style=3)
    count = 0
    time_start = proc.time()
  }
  theta = array(NA,dim = c(Ntheta,model$dimtheta,T)) #Ntheta particles of dimension dimtheta for T steps
  w = array(NA,dim = c(Ntheta,T)) #unnormalized weights for particles theta
  Wtheta = array(NA,dim = c(Ntheta,T)) #normalized weights for particles theta
  X = array(NA,dim = c(Nx,Ntheta)) #Nx particles of dimension 1 for each theta (most recent)
  z = array(NA,dim = c(Ntheta)) #matrix of likelihood estimates
  Wx = array(NA,dim = c(Nx,Ntheta)) #matrix of normalized weights for X (most recent)
  ESS = array(NA,dim = c(T)) #ESS at successive times t
  Hscore = array(NA,dim = c(T)) #prequential Hyvarinen score at successive times t
  #initialization t = 1
  theta[,,1] = model$rprior(Ntheta)
  for (m in 1:Ntheta) {
    X[,m] = model$rinitial(theta[m,,1],Nx)
    for (n in 1:Nx) {
      Wx[n,m] = model$dobs(y[1],X[n,m],1,theta[m,,1])
    }
    sum_Wx = sum(Wx[,m])
    meanWx = sum_Wx/Nx
    z[m] = meanWx
    w[m,1] = meanWx
    #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
    if (sum_Wx>0) {
      Wx[,m] = Wx[,m]/sum_Wx
    }
    else {
      Wx[,m] = 1/Nx #WARNING: NOT SURE ABOUT THIS CASE !!
    }
    if (progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
    }
  }
  #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
  sum_w = sum(w[,1])
  if (sum_w>0) {
    Wtheta[,1] = w[,1]/sum_w
  }
  else {
    Wtheta[,1] = 1/Ntheta #WARNING: NOT SURE ABOUT THIS CASE !!
  }
  Hscore[1] = hincrementContinuous(1,model,y[1],theta[,,1],Wtheta[,1],X,Wx,Ntheta,Nx)
  ESS[1] = getESS(Wtheta[,1])
  for (t in 2:T) {
    theta[,,t] = theta[,,t-1]
    for (m in 1:Ntheta) {
      A = resamplingSMC(Nx,Wx[,m])
      X[,m] = model$rtransition(X[A,m],t,theta[m,,t])
      for (n in 1:Nx) {
        Wx[n,m] = model$dobs(y[t],X[n,m],t,theta[m,,t])
      }
      sum_Wx = sum(Wx[,m])
      meanWx = sum_Wx/Nx
      z[m] = z[m]*meanWx
      w[m,t] = w[m,t-1]*meanWx
      #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
      if (sum_Wx>0) {
        Wx[,m] = Wx[,m]/sum_Wx
      }
      else {
        Wx[,m] = 1/Nx #WARNING: NOT SURE ABOUT THIS CASE !!
      }
      if (progress) {
        count = count + 1
        setTxtProgressBar(progbar, count)
      }
    }
    #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
    sum_w = sum(w[,t])
    if (sum_w>0) {
      Wtheta[,t] = w[,t]/sum_w
    }
    else {
      Wtheta[,t] = 1/Ntheta #WARNING: NOT SURE ABOUT THIS CASE !!
    }
    Hscore[t] = Hscore[t-1] + hincrementContinuous(t,model,y[t],theta[,,t],Wtheta[,t],X,Wx,Ntheta,Nx)
    ESS[t] = getESS(Wtheta[,t])
    if ((ESS[t]/Ntheta) < 0.5) {
      #compute parameters for proposal move step
      covariance = cov.wt(theta[,,t],wt = Wtheta[,t])
      mean_t = covariance$center
      cov_t = matrix(nearPD(covariance$cov)$mat,nrow = model$dimtheta)
      #resample
      resampled_index = sample(1:Ntheta,Ntheta,replace = TRUE,prob = Wtheta[,t])
      for (i in 1:Ntheta) {
        theta_old = theta[resampled_index[i],,t]
        theta_new = rmvnorm(n = 1,mean_t,cov_t)
        z_old = z[resampled_index[i]]
        prior_theta_new = model$dprior(theta_new)
        if (prior_theta_new == 0) {
          theta[i,,t] = theta_old
          X[,i] = X[,resampled_index[i]]
          Wx[,i] = Wx[,resampled_index[i]]
          z[i] = z_old
          next 
        }
        else {
          PF = smc(Nx,y[1:t],model,theta_new,resamplingSMC, proposalSMC, smoothing = FALSE)
          z_new = PF$py_hat
          num = prior_theta_new*(z_new)*(dmvnorm(theta_old,mean_t,cov_t))
          denom = (model$dprior(theta_old))*(z_old)*(dmvnorm(theta_new,mean_t,cov_t))
          a_ratio = min(1,num/denom)
          u = runif(1)
          if (u <= a_ratio) {
            theta[i,,t] = theta_new
            X[,i] = PF$Xfilter[,t]
            Wx[,i] = PF$W[,t]
            z[i] = z_new
          }
          else {
            theta[i,,t] = theta_old
            X[,i] = X[,resampled_index[i]]
            Wx[,i] = Wx[,resampled_index[i]]
            z[i] = z_old
          }
        }
      }
      w[,t] = 1
      Wtheta[,t] = 1/Ntheta
    }
  }
  if (progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("Hscore: T = ",toString(T),", Ntheta = ",toString(Ntheta),
              ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return (list(Hscore = Hscore, theta = theta, Wtheta = Wtheta, ESS = ESS))
}

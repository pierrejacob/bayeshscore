#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains functions to compute the prequential Hyvarinen score in DISCRETE case
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# WARNING: the current code only works for univariate states !!!
library(Matrix) #contains nearPD (nearest positive definite matrix)
library(mvtnorm) #contains Multivariate Normal
source("functions_resampling.R")
source("functions_proposal.R")
source("functions_smc.R")


# This function computes the ESS given a set of normalized weights
getESS = function(normalized_weights) {
  return (1/sum(normalized_weights^2))
}

# This function computes the approximation qt_hat(y) for a given set of particles
qthat = function(t,model,y,theta,Wtheta,Xpred,WXpred,Ntheta,Nx) {
  d = model$dimY
  qt = 0
  for (m in 1:Ntheta){
    for (n in 1:Nx) {
      qt = qt + Wtheta[m]*WXpred[n,m]*model$dobs(y,Xpred[n,m],t,theta[m,])
    }
  }
  return (qt)
}

# This function computes the partial score term SBk
# a,b are vectors of componentwise lower and upper bounds
# d is the dimension of y
SBk = function(k,a,b,d,y,qy_minusek,qy,qy_plusek) {
  if (y[k]==b[k]) {
    return (-2*(qy-qy_minusek)/qy_minusek)
  }
  else {
    if (y[k]==a[k]) {
      return (2*(qy_plusek-qy)/qy + ((qy_plusek-qy)/qy)^2)
    }
    else {
      return (2*((qy_plusek-qy)/qy-(qy-qy_minusek)/qy_minusek) + ((qy_plusek-qy)/qy)^2)
    }
  }
}

# This function computes the partial score term SHd
# a,b are vectors of componentwise lower and upper bounds
# d is the dimension of y
SHd = function(t,model,yt,theta,Wtheta,Xpred,WXpred,Ntheta,Nx) {
  d = model$dimY
  result = 0
  for (k in 1:d) {
    ek = rep(0,d)
    ek[k] = 1
    qy = qthat(t,model,yt,theta,Wtheta,Xpred,WXpred,Ntheta,Nx)
    qy_minusek = qthat(t,model,yt-ek,theta,Wtheta,Xpred,WXpred,Ntheta,Nx)
    qy_plusek = qthat(t,model,yt+ek,theta,Wtheta,Xpred,WXpred,Ntheta,Nx)
    result = result + SBk(k,model$a,model$b,d,yt,qy_minusek,qy,qy_plusek)
  }
  return (result)
}

# This function compute the total prequential Hyvarinen Score
hyvarinenDiscrete = function(y, T, model, Ntheta, Nx, timesteps = NA, smc = SMCmvtY, 
                             resamplingSMC = multinomial, proposalSMC = transition, progress = TRUE) {
  if (is.na(timesteps[1])) {
    timesteps = 1:T
  }
  if (progress) {
    progbar = txtProgressBar(min = 0,max = Ntheta*T,style=3)
    count = 0
    time_start = proc.time()
  }
  theta = array(NA,dim = c(Ntheta,model$dimtheta,T)) #particles for p(theta|y_1...t-1)
  w = array(NA,dim = c(Ntheta,T)) #unnormalized weights for particles theta
  wtheta = array(NA,dim = c(Ntheta,T)) #normalized weights for particles theta
  X = array(NA,dim = c(Nx,Ntheta)) #Nx particles of dimension 1 for each theta (most recent)
  Xpred = array(NA,dim = c(Nx,Ntheta)) #particles for p(xt|y_1...t-1) (predictive)
  Wx = array(NA,dim = c(Nx,Ntheta)) #normalized weights for X (most recent)
  WXpred = array(NA,dim = c(Nx,Ntheta)) #normalized weights for Xpred (predictive)
  z = array(NA,dim = c(Ntheta)) #matrix of likelihood estimates
  ESS = array(NA,dim = c(T)) #ESS at successive times t
  Hscore = array(NA,dim = c(T)) #prequential Hyvarinen score at successive times t
  #initialization t = 0
  thetaprior = model$rprior(Ntheta)
  wthetaprior = rep(1/Ntheta,Ntheta)
  theta[,,1] = thetaprior
  for (m in 1:Ntheta) {
    X[,m] = model$rinitial(theta[m,,1],Nx)
    Wx[,m] = model$dobs(y[,1],X[,m],1,theta[m,,1])
    Xpred[,m] = X[,m]
    WXpred[,m] = 1/Nx
    sum_Wx = sum(Wx[,m])
    meanWx = sum_Wx/Nx
    z[m] = meanWx
    w[m,1] = meanWx
    #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
    if (sum_Wx>0) {
      Wx[,m] = Wx[,m]/sum_Wx
    }
    else {
      Wx[,m] = rep(1/Nx,Nx) #WARNING: NOT SURE ABOUT THIS CASE !!
    }
    if (progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
    }
  }
  #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
  sum_w = sum(w[,1])
  if (sum_w>0) {
    wtheta[,1] = w[,1]/sum_w
  }
  else {
    wtheta[,1] = 1/Ntheta #WARNING: NOT SURE ABOUT THIS CASE !!
  }
  Hscore[1] = SHd(1,model,y[,1],thetaprior,wthetaprior,Xpred,WXpred,Ntheta,Nx)
  ESS[1] = getESS(wtheta[,1])
  for (t in 2:T) {
    theta[,,t] = theta[,,t-1]
    dt = timesteps[t]-timesteps[t-1]
    for (m in 1:Ntheta) {
      Xpred[,m] = model$rtransition(X[,m],dt,theta[m,,t-1])
      WXpred[,m] = Wx[,m]
      A = resamplingSMC(Nx,Wx[,m])
      X[,m] = model$rtransition(X[A,m],dt,theta[m,,t])
      Wx[,m] = model$dobs(y[,t],X[,m],t,theta[m,,t])
      sum_Wx = sum(Wx[,m])
      meanWx = sum_Wx/Nx
      z[m] = z[m]*meanWx
      w[m,t] = w[m,t-1]*meanWx
      #Sometimes the weights become all 0 leading to NaN when trying to renormalize...
      if (sum_Wx>0) {
        Wx[,m] = Wx[,m]/sum_Wx
      }
      else {
        Wx[,m] = rep(1/Nx,Nx) #WARNING: NOT SURE ABOUT THIS CASE !!
      }
      if (progress) {
        count = count + 1
        setTxtProgressBar(progbar, count)
      }
    }
    sum_w = sum(w[,t])
    if (sum_w>0) {
      wtheta[,t] = w[,t]/sum_w
    }
    else {
      wtheta[,t] = 1/Ntheta #WARNING: NOT SURE ABOUT THIS CASE !!
    }
    Hscore[t] = Hscore[t-1] + SHd(t,model,y[,t],theta[,,t-1],wtheta[,t-1],Xpred,WXpred,Ntheta,Nx)
    ESS[t] = getESS(wtheta[,t])
    if ((ESS[t]/Ntheta) < 0.5) {
      #compute parameters for proposal move step
      covariance = cov.wt(theta[,,t],wt = wtheta[,t])
      mean_t = covariance$center
      cov_t = matrix(nearPD(covariance$cov)$mat,nrow = model$dimtheta)
      #resample
      resampled_index = sample(1:Ntheta,Ntheta,replace = TRUE,prob = wtheta[,t])
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
          PF = smc(Nx,y[,1:t],model,theta_new,resamplingSMC, proposalSMC, smoothing = FALSE)
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
      wtheta[,t] = 1/Ntheta
    }
  }
  if (progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("Hscore: T = ",toString(T),", Ntheta = ",toString(Ntheta),
              ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return (list(Hscore = Hscore, theta = theta, wtheta = wtheta, ESS = ESS))
}

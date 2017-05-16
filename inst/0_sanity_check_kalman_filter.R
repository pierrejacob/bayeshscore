library(HyvarinenSSM)

#create data
nobservations <- 100
model <- get_model_lineargaussian()
theta = c(0.5,1,1,1)
sim = simulateData(model, theta = theta, nobservations)
X = sim$X
Y = sim$Y

# set parameters
phi = theta[1]
psi = theta[2]
sigmaV2 = theta[3]
sigmaW2 = theta[4]
initial_mean = 0
initial_var = (sigmaW2)/(1-phi^2)

# run Kalman filter
KF = KF_filtering(Y,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)

# run Kalman filter Rcpp
kalman_module <- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
Kalman <- new(kalman_module$Kalman)
Kalman$set_parameters(list(rho = phi, sigma = sqrt(sigmaW2), eta = psi, tau = sqrt(sigmaV2)))
Kalman$set_observations(matrix(Y, ncol = 1))
Kalman$filtering()


# run particle filter
algorithmic_parameters = list()
algorithmic_parameters$Nx = 2^10
algorithmic_parameters$resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1))
BPF = bootstrap_particle_filter(Y,model,theta,algorithmic_parameters)

# run particle filter CPF (equivalent to bootstrap PF when no conditioning path)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
CPF = conditional_particle_filter(Y,model,theta,algorithmic_parameters$Nx)

# Check log-likelihood
ggplot() +
  geom_line(aes(1:nobservations,cumsum(sapply(1:nobservations,function(t)KF_logdpredictive(Y[,t,drop=FALSE],t, KF)))),col='blue',size=2,linetype=2,alpha=0.6) +
  geom_line(aes(1:nobservations,cumsum(Kalman$get_incremental_ll())),col='red',size=1) +
  geom_point(aes(1:nobservations,cumsum(BPF$incremental_ll)),size=2) +
  geom_point(aes(1:nobservations,cumsum(CPF$incremental_ll)),size=3,shape=3)

# Check filtering means
ggplot() +
  geom_line(aes(1:nobservations,sapply(1:nobservations,function(t)KF[[t]]$muX_t_t)),col='blue',size=2,linetype=2,alpha=0.6)+
  geom_line(aes(1:nobservations,sapply(1:nobservations,function(t)Kalman$get_filtering_mean(t))),col='red',size=1)+
  geom_point(aes(1:nobservations,sapply(1:nobservations,function(t)sum(BPF$X_history[[t]]*BPF$weight_history[[t]]))),size=3)

# # Plot some paths
# path.df = data.frame()
# for (i in 1:algorithmic_parameters$Nx) {
#   path.df = rbind(path.df,data.frame(x=c(CPF$tree$get_path(i-1)),time=1:nobservations,index=i))
# }
# ggplot(path.df) + geom_line(aes(time,x,group=index))



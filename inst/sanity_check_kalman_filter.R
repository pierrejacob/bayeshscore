library(HyvarinenSSM)
# set.seed(17)

nobs <- 50
model <- get_model_lineargaussian()
theta = c(0.5,1,1,1)
sim = simulateData(model, theta = theta, nobs)
X = sim$X
Y = sim$Y


library(dlm)
phi <- theta[1]
psi = theta[2]
sigmaV2 = theta[3]
sigmaW2 <- theta[4]
initial_mean = 0
initial_var = (sigmaW2)/(1-phi^2)


# run Kalman filter
KF = KF_filtering(Y,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)

# run particle filter
algorithmic_parameters = list()
algorithmic_parameters$Nx = 2^10
algorithmic_parameters$resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1))
BPF = bootstrap_particle_filter(Y,model,theta,algorithmic_parameters)

ggplot() +
  geom_line(aes(1:nobs,sapply(1:nobs,function(t)KF[[t]]$muX_t_t)),col='blue') +
  geom_point(aes(1:nobs,sapply(1:nobs,function(t)sum(BPF$X_history[[t]]*BPF$weight_history[[t]]))),size=2)





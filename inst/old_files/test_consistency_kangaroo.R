library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
library(wesanderson)

# Define model and data
nobservations <- 20
timestep = seq(from = 0, by = 0.25, length.out = nobservations)
theta_sim1 = c(0.75,0.05,2,0.0025)
theta_sim2 = c(0.15,0.05,0.001)
theta_sim3 = c(0.15,0.05)
model1 <- get_model_kangarooLogistic
model2 <- get_model_kangarooExponential
model3 <- get_model_kangarooRandomwalk
all_models = list(model1,model2,model3)
all_thetas = list(theta_sim1,theta_sim2,theta_sim3)

#True data-generating model and parameter
index = 1
model = all_models[[index]](timestep)
theta = all_thetas[[index]]


#Generate artificial data
sim = simulateData(model, theta = theta, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)# observations in a matrix of dimensions dimy x nobservations

#Plot states
observations.df = data.frame(time = 1:nobservations, X = t(X),Y1 = Y[1,],Y2 = Y[2,])
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=X),size=2) +
  geom_line(aes(y=X)) +
  xlab("\n Time")
plot(g)

# Plot data
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=Y1),size=2) +
  geom_point(aes(y=Y2),size=2) +
  geom_segment(aes(y = Y1, xend = time, yend = Y2),linetype="dashed") +
  xlab("\n Time")
plot(g)

# Define algorithmic parameters for each model
Ntheta = 2^7
Nx = 2^5
algorithmic_parameters = list(Ntheta = Ntheta, Nx = Nx,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE)

repl = 5
results.df = data.frame()
posterior.df = data.frame()

registerDoParallel(cores=5)

for (i in 1:length(all_models)){
  M = all_models[[i]]
  print(paste("Started at:",Sys.time()))
  time_start = proc.time()
  results = foreach(i=1:repl,.packages='HyvarinenSSM',.verbose = TRUE) %dopar% {
    hscore_discrete(observations, M(timestep), algorithmic_parameters)
  }
  time_end = proc.time()-time_start
  cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
            ", Nx = ",toString(Nx),"\n",sep = ""))
  print(time_end)

  for (r in 1:repl){
    results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                              logevidence = results[[r]]$logevidence,
                                              hscore = results[[r]]$hscore,
                                              rep = r,
                                              model = i))
    posterior.df = rbind(posterior.df, data.frame(sigma = results[[r]]$thetas[,1],
                                                  tau = results[[r]]$thetas[,2],
                                                  r = tryCatch({results[[r]]$thetas[,3]},error=function(e){NA}),
                                                  b = tryCatch({results[[r]]$thetas[,4]},error=function(e){NA}),
                                                  w = results[[r]]$thetanormw,
                                                  rep=r,
                                                  model = i))
  }
}
cat(paste("True Model =",toString(index)))

g = ggplot(results.df, aes(x = time, y = logevidence, group = interaction(rep,model),colour = factor(model))) +
  geom_line(size=1.2) + scale_color_manual(name= expression(bold(Model)),values = wes_palette("FantasticFox")[c(3,4,5)]) +
  ylab("Log-evidence") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=24,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold")) + xlab("\n Number of observations")
plot(g)

g = ggplot(results.df, aes(x = time, y = hscore, group = interaction(rep,model),colour = factor(model))) +
  geom_line(size=1.2) + scale_color_manual(name= expression(bold(Model)),values = wes_palette("FantasticFox")[c(3,4,5)]) +
  ylab("Prequential Hyvarinen score") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=24,face="bold"),
        legend.text=element_text(size=24,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold")) + xlab("\n Number of observations")
plot(g)

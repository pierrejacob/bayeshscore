library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
library(wesanderson)


#=======================================================================
#=======================================================================
#=======================================================================
#====== WARNING: current crashes when trying to use parallel foreach ===
#=======================================================================
#========= (not sure how to properly export the Cpp Class TreeClass ) ==
#=======================================================================


# Define model and data
model1 = get_model_kangarooLogistic()
model2 = get_model_kangarooExponential()
model3 = get_model_kangarooRandomwalk()
all_models = list(model1,model2,model3)
dataset = data_kangaroo
observations = dataset[1:2,]
nobservations = ncol(observations)


# Define algorithmic parameters for each model
Ntheta = 2^5
Nx = NULL
min_acceptance_rate = 0.30
algorithmic_parameters = list(Ntheta = Ntheta, Nx = Nx,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE, min_acceptance_rate = min_acceptance_rate)
repl = 1
results.df = data.frame()
posterior.df = data.frame()



module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
# algorithmic_parameters = list(Ntheta = Ntheta, Nx = Nx,
#                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
#                               progress = TRUE, TreeClass = TreeClass)
all_results = list()
for (i in 1:length(all_models)){
  all_results[[i]] = list()
  M = all_models[[i]]
  for (r in 1:repl){
    print(paste("Model",toString(i)))
    results = hscore_discrete(observations, M, algorithmic_parameters)
    all_results[[i]][[r]] = results
    results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                              logevidence = results$logevidence,
                                              hscore = results$hscore,
                                              rep = r,
                                              model = i))
    posterior.df = rbind(posterior.df, data.frame(sigma = results$thetas[,1],
                                                  tau = results$thetas[,2],
                                                  r = tryCatch({results$thetas[,3]},error=function(e){NA}),
                                                  b = tryCatch({results$thetas[,4]},error=function(e){NA}),
                                                  w = results$thetanormw,
                                                  rep=r,
                                                  model = i))
  }
}


#=======================================================================
#=======================================================================
#=======================================================================
#== THE FOLLOWING BLOCK IS THE PARALLEL VERSION BUT ERROR WITH TREES ===
#=======================================================================
#========= (not sure how to properly export the Cpp Class TreeClass ) ==
#=======================================================================
# registerDoParallel(cores=5)
#
# for (i in 1:length(all_models)){
#   M = all_models[[i]]
#   print(paste("Started at:",Sys.time()))
#   time_start = proc.time()
#   results = foreach(i=1:repl,.packages='HyvarinenSSM',.verbose = TRUE) %dopar% {
#
#     ####### THIS DOES NOT WORK !!!!! "error: objet TreeClass introuvable" #######
#     module_tree <- Module("module_tree", PACKAGE = "HyvarinenSSM")
#     TreeClass <- module_tree$Tree
#     ####### THIS DOES NOT WORK !!!!! "error: objet TreeClass introuvable" #######
#     hscore_discrete(observations, M, algorithmic_parameters)
#   }
#   time_end = proc.time()-time_start
#   cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
#             ", Nx = ",toString(Nx),"\n",sep = ""))
#   print(time_end)
#
#   for (r in 1:repl){
#     results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
#                                               logevidence = results[[r]]$logevidence,
#                                               hscore = results[[r]]$hscore,
#                                               rep = r,
#                                               model = i))
#     posterior.df = rbind(posterior.df, data.frame(sigma = results[[r]]$thetas[,1],
#                                                   tau = results[[r]]$thetas[,2],
#                                                   r = tryCatch({results[[r]]$thetas[,3]},error=function(e){NA}),
#                                                   b = tryCatch({results[[r]]$thetas[,4]},error=function(e){NA}),
#                                                   w = results[[r]]$thetanormw,
#                                                   rep=r,
#                                                   model = i))
#   }
# }

#Plot log-evidence across time
g = ggplot(results.df, aes(x = time, y = -logevidence, group = interaction(rep,model),colour = factor(model))) +
  geom_line(size=1.2) +
  scale_color_manual(name= "Model",values = wes_palette("Darjeeling")[c(1,3,5)]) +
  ylab("Negative Log-evidence") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold")) + xlab("\n Number of observations")
plot(g)

#Plot prequential hscore across time
g = ggplot(results.df, aes(x = time, y = hscore, group = interaction(rep,model),colour = factor(model))) +
  geom_line(size=1.2) +
  scale_color_manual(name= "Model",values = wes_palette("Darjeeling")[c(1,3,5)]) +
  ylab("Prequential Hyvarinen score") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold")) + xlab("\n Number of observations")
plot(g)

#Boxplot final log-evidence and final prequential hscore
final_logevidence = subset(results.df,time==41)[,c("logevidence","model")]
final_preq_hscore = subset(results.df,time==41)[,c("hscore","model")]
g = ggplot(final_logevidence,aes(model,-logevidence,group=model,colour=factor(model))) +
  geom_boxplot(outlier.shape = NA,size=1) +
  scale_color_manual(name= "Model",values = wes_palette("Darjeeling")[c(1,3,5)]) +
  ylab("Final negative log-evidence") +
  xlab("") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
plot(g)
g = ggplot(final_preq_hscore,aes(model,hscore,group=model,colour=factor(model))) +
  geom_boxplot(outlier.shape = NA,size=1) +
  scale_color_manual(name= "Model",values = wes_palette("Darjeeling")[c(1,3,5)]) +
  ylab("Final prequential hscore") +
  xlab("") +
  guides(colour = guide_legend(title.hjust = 0.2,keywidth = 3, keyheight = 2)) +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
plot(g)



posterior1 = subset(posterior.df,model==1)
g1 = ggplot(posterior1, aes(x = r, weight = w, group = rep,colour = rep)) +
  stat_density(size=0.8,geom='line',position="identity") +
  ylab("") +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.position='none') +
  xlab(expression(atop("",r)))
g2 = ggplot(posterior1, aes(x = sigma, weight = w, group = rep,colour = rep)) +
  stat_density(size=0.8,geom='line',position="identity") +
  ylab("") +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"),
        legend.position='none') +
  xlab(expression(atop("",bold(sigma))))
g3 = ggplot(posterior1, aes(x = tau, weight = w, group = rep,colour = rep)) +
  stat_density(size=0.8,geom='line',position="identity") +
  ylab("") +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold")) +
  xlab(expression(atop("",bold(tau))))
g4 = ggplot(posterior1, aes(x = b, weight = w, group = rep,colour = rep)) +
  stat_density(size=0.8,geom='line',position="identity") +
  ylab("") +
  theme(legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20,face="bold"),
        axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold")) +
  xlab(expression(atop("",bold(b))))
grid.arrange(g1,g4,g2,g3, ncol = 2, nrow = 2)

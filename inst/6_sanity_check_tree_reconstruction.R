##################################################################################################
# This checks the deconstruction and reconstruction of (RCpp) tree objects
##################################################################################################
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
#--------------------------------------------------------------------------------------------
# create data
nobservations = 5
model = get_model_lineargaussian()
theta_star = c(0.8,1,1,1)
#--------------------------------------------------------------------------------------------
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^5
algorithmic_parameters$Nx_max = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.8 # purposely set high to trigger increase Nx step and test Nx_max
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
### Run SMC_2
model_nolikelihood = model
model_nolikelihood$likelihood = NULL # this forces the use of SMC2
model_nolikelihood$dpredictive = NULL # this forces the use of SMC2
smc2_results = hscore(observations,model_nolikelihood,algorithmic_parameters)
#---------------------------------------------------------------------------------------------
######################################################################################
# Define a function that checks if two trees are identical
######################################################################################
is.sametree = function(tree1, tree2){
  return (tree1$N == tree2$N &&
            tree1$M == tree2$M &&
            tree1$dimx == tree2$dimx &&
            tree1$nsteps == tree2$nsteps &&
            sum(tree1$a_star == tree2$a_star) == tree1$M &&
            sum(tree1$o_star == tree2$o_star) == tree1$M &&
            sum(tree1$x_star == tree2$x_star) == tree1$M &&
            sum(tree1$l_star == tree2$l_star) == tree1$N)
}

######################################################################################
# Take a single arbitrary tree, break it down, and reconstruct it
######################################################################################
example_tree = smc2_results$PF_history[[2]][[1]]$tree
example_tree_attributes = tree_getattributes(example_tree)
example_tree_reconstructed = tree_reconstruct(example_tree_attributes)
# Check that the reconstructed tree is identical to the initial tree
is.sametree(example_tree, example_tree_reconstructed)
#---------------------------------------------------------------------------------------------
######################################################################################
# Take a list of trees, break them down, and reconstruct them
######################################################################################
example_trees = lapply(1:algorithmic_parameters$Ntheta,function(i)smc2_results$PF_history[[2]][[i]]$tree)
example_trees_attributes = trees_getattributes(example_trees)
example_trees_reconstructed = trees_reconstruct(example_trees_attributes)
# Check that the reconstructed tree is identical to the initial tree
check = sum(sapply(1:length(example_trees),function(i)is.sametree(example_trees[[i]],example_trees_reconstructed[[i]])))
cat("Total number of trees : ",length(example_trees),"\n","Number of same trees : ", check, "\n",sep="")
if (check==length(example_trees)) {cat("All the trees match !\n")}

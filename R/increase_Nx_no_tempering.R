#'@rdname increase_Nx_no_tempering
#'@title increase_Nx_no_tempering
#'@description This function runs a CPF in order to increase Nx (the number of particles). The default new number of particles is set to 2*Nx.
#'@export
increase_Nx_no_tempering <- function(observations, t, model, thetas, xnormW, trees, algorithmic_parameters){
  # Extract algorithmic parameters
  Ntheta = algorithmic_parameters$Ntheta
  Nx = algorithmic_parameters$Nx
  Nx_new <- 2*Nx #by default we multiply the number of particles by 2
  algorithmic_parameters$Nx = Nx_new
  # Initialize empty arrays of larger size
  X_new = array(NA,dim = c(Nx_new, model$dimX, Ntheta)) #Nx particles (most recent) for each theta (size = Nx,dimX,Ntheta)
  xnormW_new = matrix(NA, nrow = Nx_new, ncol = Ntheta) #matrix of corresponding normalized X-weights (size = Nx,Ntheta)
  log_z_new = rep(0, Ntheta) #matrix of log-likelihood estimates (size = Ntheta)
  # Construct list of trees to store paths (one tree for each theta)
  trees_new = list()
  for (i in 1:Ntheta){
    trees_new[[i]] = new(TreeClass, Nx_new, 10*Nx_new, model$dimX)
  }
  # Perform a conditional particle filter for each theta
  for (i in 1:Ntheta){
    current_tree = trees[[i]]
    current_path = current_tree$get_path(sample(x = 0:(Nx-1), size = 1, replace = TRUE, prob = xnormW[,i]))
    cpf = conditional_particle_filter(matrix(observations[,1:t],ncol = t), model, thetas[,i], Nx_new, path = current_path)
    X_new[,,i] = cpf$X
    xnormW_new[,i] = cpf$xnormW
    log_z_new[i] = cpf$log_p_y_hat
    trees_new[[i]] = cpf$tree
  }
  return(list(log_z = log_z_new, X = X_new, xnormW = xnormW_new, trees = trees_new, new_Nx = Nx_new))
}

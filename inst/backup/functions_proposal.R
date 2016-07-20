#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains different proposal for particle filter
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# This proposal corresponds to the Markov transition kernel of the latent states
transition = function(Xt,y,t,model,theta) {
  return (model$rtransition(Xt,t,theta))
}
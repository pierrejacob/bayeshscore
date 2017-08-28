#'@rdname get_batch_logevidence
#'@title get_batch_logevidence
#'@description This function computes the logevidence of a batch of observation Y1:b via importance sampling
#'@export
get_batch_logevidence = function(observations, model, proposal, algorithmic_parameters) {
  Ntheta = algorithmic_parameters$Ntheta
  b = algorithmic_parameters$initialbatchsize
  observations_1_to_b = matrix(observations[,1:b], nrow = model$dimY)
  thetas = proposal$r(Ntheta)
  log_z = vector("numeric",Ntheta)
  log_w = vector("numeric",Ntheta)
  for (i in 1:Ntheta) {
    if (model$dprior(thetas[i,],log = TRUE) == -Inf) {
      log_w[i] = -Inf
      next
    }
    else {
      smc = bootstrap_particle_filter(observations_1_to_b, model, thetas[i,], algorithmic_parameters)
      log_z[i] = smc$log_p_y_hat
      log_w[i] = model$dprior(thetas[i,],log = TRUE) + log_z[i] - proposal$logd(thetas[i,])
    }
  }
  log_w_max = max(log_w)
  w = exp(log_w - log_w_max)
  logevidence_batch = log(sum(w)/Ntheta) + log_w_max
  return (logevidence_batch)
}


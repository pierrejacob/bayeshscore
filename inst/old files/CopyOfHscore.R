#'@rdname hscore
#'@title hscore
#'@description This function computes successive prequential Hyvarinen score by running SMC or SMC2.
#'It is a wrapper of the four different cases: continuous or discrete / SMC or SMC2.
#'It also computes the successive log-evidence as a by-product.
#'@export
hscore <- function(observations, model, algorithmic_parameters){
  if (is.null(model$observation_type)){
    cat('ERROR: please specify observation_type in the model (\'dicrete\' or \'continuous\')\n')
    return(NULL)
  } else {
    observation_type = tolower(model$observation_type)
    intractable_likelihood = (is.null(model$likelihood))&(is.null(model$dpredictive))
    algorithmic_parameters = set_default_algorithmic_parameters(observations,model,algorithmic_parameters)
    model = set_default_model(model)
    if (intractable_likelihood) {
      cat(paste('Hscore: type = ',observation_type,', Method = SMC2, ',
                'Ntheta = ',toString(algorithmic_parameters$Ntheta),
                ', Nx = ',toString(algorithmic_parameters$Nx),'\n',sep = ""))
      if(observation_type == 'continuous'){return(hscore_continuous_smc2(observations, model, algorithmic_parameters))}
      else if (observation_type == 'discrete'){return(hscore_discrete_smc2(observations, model, algorithmic_parameters))}
      else {cat('ERROR: observation_type must be \'dicrete\' or \'continuous\'\n'); return(NULL)}
    } else {
      cat(paste('Hscore: type = ', observation_type, ', Method = SMC, ',
                'Ntheta = ',toString(algorithmic_parameters$Ntheta), '\n',sep = ""))
      if(observation_type == 'continuous'){return(hscore_continuous_smc(observations, model, algorithmic_parameters))}
      else if (observation_type == 'discrete'){return(hscore_discrete_smc(observations, model, algorithmic_parameters))}
      else {cat('ERROR: observation_type must be \'dicrete\' or \'continuous\'\n'); return(NULL)}
    }
  }
}

#'@rdname hscore
#'@title hscore
#'@description This function computes successive prequential Hyvarinen score by running SMC or SMC2.
#'It is a wrapper of the four different cases: continuous or discrete / SMC or SMC2.
#'It also computes the successive log-evidence as a by-product.
#'@export
hscore <- function(observations, model, algorithmic_parameters){
  observation_type = model$observation_type
  intractable_likelihood = (is.null(model$likelihood))&(is.null(model$dpredictive))
  if (is.null(observation_type)){
    cat('ERROR: please specify observation_type in the model (\'dicrete\' or \'continuous\')\n')
    return(NULL)
  } else {
    if (tolower(observation_type) == 'continuous'){
      if (intractable_likelihood) {
        cat('Hscore for continuous observations using SMC2\n')
        return(hscore_continuous_smc2(observations, model, algorithmic_parameters))
      } else {
        cat('Hscore for continuous observations using SMC\n')
        return(hscore_continuous_smc(observations, model, algorithmic_parameters))
      }
    }
    else if (tolower(observation_type) == 'discrete'){
      if (intractable_likelihood) {
        cat('Hscore for discrete observations using SMC2\n')
        return(hscore_discrete_smc2(observations, model, algorithmic_parameters))
      } else {
        cat('Hscore for discrete observations using SMC\n')
        return(hscore_discrete_smc(observations, model, algorithmic_parameters))
      }
    }
    else {
      cat('ERROR: observation_type must be \'dicrete\' or \'continuous\'\n')
      return(NULL)
    }
  }
}

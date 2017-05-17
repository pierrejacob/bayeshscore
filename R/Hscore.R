#'@rdname hscore
#'@title hscore
#'@description This function computes successive prequential Hyvarinen score by running SMC or SMC2.
#'It is a wrapper of the four different cases: continuous or discrete / SMC or SMC2.
#'It also computes the successive log-evidence as a by-product.
#'@export
hscore <- function(observations, model, algorithmic_parameters){
  if (is.null(model$observation_type)||!(tolower(model$observation_type)%in%c("discrete","continuous"))){
    cat("ERROR: please specify observation_type in the model (\"dicrete\" or \"continuous\")\n")
    return(NULL)
  } else {
    intractable_likelihood = (is.null(model$likelihood))&&(is.null(model$dpredictive))
    # set missing algorithmic parameters to default values
    algorithmic_parameters = set_default_algorithmic_parameters(observations,model,algorithmic_parameters)
    algorithmic_parameters$hscore = TRUE
    # set missing fields to automatic values (e.g. numerical derivatives)
    model = set_default_model(model)
    cat(paste("Hscore: type = ",tolower(model$observation_type),", Method = ",sep=""))
    # Compute hscore using SMC2 (likelihood unavailable) or SMC (likelihood available)
    if (intractable_likelihood) {
      cat(paste("SMC2", "Ntheta = ",toString(algorithmic_parameters$Ntheta),
                ", Nx (initial) = ",toString(algorithmic_parameters$Nx),"\n",sep = ""))
      return (smc2(observations, model, algorithmic_parameters))
    } else {
      cat(paste("SMC", "Ntheta = ",toString(algorithmic_parameters$Ntheta),"\n",sep = ""))
      return (smc(observations, model, algorithmic_parameters))
    }
  }
}

#'@rdname hscore_resume
#'@title hscore_resume
#'@description This function resumes a previously interrupted SMC / SMC2 run to compute the hscore.
#'It is a wrapper of the four different cases: continuous or discrete / SMC or SMC2.
#'It also computes the successive log-evidence as a by-product.
#'@export
hscore_resume <- function(RDSsave=NULL, savefilename=NULL, next_observations=NULL, new_algorithmic_parameters=NULL) {
  ########################################################################
  ################ Load partial results from save file ###################
  ########################################################################
  # load RDS file
  if (is.null(RDSsave)) {
    if (is.null(savefilename)) {
      cat("partial results (loaded from RDS file) or path to RDS file (savefilename) must be provided\n")
      return (NULL)
    } else {
      RDSsave = readRDS(savefilename)
    }
  }
  method = RDSsave$method
  cat("Hscore resume: type = ",tolower(RDSsave$model$observation_type),", Method = ",method,"\n", sep="")
  if (method == "SMC"){
    return (smc_resume(RDSsave=RDSsave, next_observations=next_observations, new_algorithmic_parameters=new_algorithmic_parameters))
  } else if (method == "SMC2") {
    return (smc2_resume(RDSsave=RDSsave, next_observations=next_observations, new_algorithmic_parameters=new_algorithmic_parameters))
  } else {
    cat("Invalid RDS save file, no method specified\n")
    return (NULL)
  }
}



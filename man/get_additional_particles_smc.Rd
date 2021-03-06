% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_variance_reduction.R
\name{get_additional_particles_smc}
\alias{get_additional_particles_smc}
\title{get_additional_particles_smc}
\usage{
get_additional_particles_smc(Nc, thetas, normw, byproducts, t, observations,
  model, logtargetdensities, algorithmic_parameters)
}
\description{
(SMC version) Generates additional particles given current ones, in order to reduce the variance
of the associated particle estimators. Nc is the desired total number of particles
to be used in the computation / estimation (so about Nc - Ntheta new particles are generated,
used in the estimation, and then discarded). The generated particles are equally weigthed.
}

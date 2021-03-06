% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_model_lineargaussiansimpler_discreteprior.R
\name{get_model_lineargaussiansimpler_discreteprior}
\alias{get_model_lineargaussiansimpler_discreteprior}
\title{get_model_lineargaussiansimpler_discreteprior}
\usage{
get_model_lineargaussiansimpler_discreteprior()
}
\description{
Univariate linear Gaussian model with one parameter (\code{sigmaV2})and discrete prior
Latent states: X[t] = \code{phi}*X[t-1] + N(0,\code{sigmaW2}).
Observations: Y[t] = \code{psi}*X[t] + N(0,\code{sigmaV2}).
}

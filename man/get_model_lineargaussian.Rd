% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_model_lineargaussian.R
\name{get_model_lineargaussian}
\alias{get_model_lineargaussian}
\title{get_model_lineargaussian}
\usage{
get_model_lineargaussian()
}
\description{
Univariate linear Gaussian model with 4 unknown parameters
(\code{phi}, \code{psi}, \code{sigmaW2}, \code{sigmaV2}).
Latent states: X[t] = \code{phi}*X[t-1] + N(0,\code{sigmaW2}).
Observations: Y[t] = \code{psi}*X[t] + N(0,\code{sigmaV2}).
}

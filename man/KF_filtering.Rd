% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_kalmanfilter.R
\name{KF_filtering}
\alias{KF_filtering}
\title{KF_filtering}
\usage{
KF_filtering(Y, phi, psi, sigmaV2, sigmaW2, initial_mean, initial_var)
}
\description{
This function runs a basic kalman filter and returns the predictive means and variances
of both states and observations, and the filtering means and variances
Prototype: state parameters = phi, sigmaV2 // observation parameter: psi, sigmaW2
}

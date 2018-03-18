#ifndef SSP_H_
#define SSP_H_
#include <RcppEigen.h>
using namespace Rcpp;

IntegerVector SSP_resampling_n_(const NumericVector & weights, int ndraws, const NumericVector & u);

#endif


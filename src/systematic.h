#ifndef SYSTEMATIC_H_
#define SYSTEMATIC_H_
#include <RcppEigen.h>
using namespace Rcpp;

IntegerVector systematic_resampling_n_(const NumericVector & weights, int ndraws, double u);
  
#endif


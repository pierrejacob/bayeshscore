#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
// Functions that compute the matrices G and h needed for estimating the derivatives of the density
// by following Sasaki, Noh, SUgiyama (2015)
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// RBF kernel
// [[Rcpp::export]]
double phi_cpp(double y, double yi, double sigma2) {
  return (exp(-pow(y-yi,2)/(2*sigma2)));
}
// first derivative RBF kernel
// [[Rcpp::export]]
double dphi_cpp(double y, double yi, double sigma2) {
  return (-phi_cpp(y,yi,sigma2)*(y-yi)/sigma2);
}
// second derivative RBF kernel
// [[Rcpp::export]]
double d2phi_cpp(double y, double yi, double sigma2) {
  return ((-1/sigma2 + pow((y-yi)/sigma2,2))*phi_cpp(y,yi,sigma2));
}
// Matrix G
// [[Rcpp::export]]
NumericMatrix get_G_RBF_cpp(NumericVector x, double sigma2) {
  int N = x.length();
  double c = sqrt(M_PI*sigma2);
  NumericMatrix result(N,N);
  for (int i = 0; i < N; i ++){
    result(i,i) = c;
    for (int j = 0; j < i; j ++){
      result(i,j) = c*phi_cpp(x(j),x(i),2*sigma2);
    }
  }
  for (int i = 0; i < N; i ++){
    for (int j = (i+1); j < N; j ++){
      result(i,j) = result(j,i);
    }
  }
  return result;
}
// Vector h for p'(y) (first derivative)
// [[Rcpp::export]]
NumericVector get_h1_RBF_cpp(NumericVector x, double sigma2) {
  int N = x.length();
  NumericMatrix h_matrix(N,N);
  NumericVector h(N);
  for (int i = 0; i < N; i ++){
    for (int j = 0; j < i; j ++){
      h_matrix(i,j) = dphi_cpp(x(j),x(i),sigma2);
    }
  }
  for (int i = 0; i < N; i ++){
    for (int j = (i+1); j < N; j ++){
      // WARNING: h1 is not symmetric because of the (yi-yj) term, hence the minus sign
      h_matrix(i,j) = -h_matrix(j,i);
    }
  }
  for (int j = 0; j < N; j ++){
    h = h + h_matrix(_,j);
  }
  return (h/N);
}
// Vector h for p''(y) (second derivative)
// [[Rcpp::export]]
NumericVector get_h2_RBF_cpp(NumericVector x, double sigma2) {
  int N = x.length();
  double diag = 1/sigma2;
  NumericMatrix h_matrix(N,N);
  NumericVector h(N);
  for (int i = 0; i < N; i ++){
    h_matrix(i,i) = diag;
    for (int j = 0; j < i; j ++){
      // Warning: there is a minus sign
      h_matrix(i,j) = -d2phi_cpp(x(j), x(i), sigma2);
    }
  }
  for (int i = 0; i < N; i ++){
    for (int j = (i+1); j < N; j ++){
      h_matrix(i,j) = h_matrix(j,i);
    }
  }
  for (int j = 0; j < N; j ++){
    h = h + h_matrix(_,j);
  }
  return (h/N);
}
// Compute the loss for the first derivative given a fitted theta, a bandwidth sigma2,
// a training sample and a test sample
// [[Rcpp::export]]
double loss_d1(NumericVector theta, double sigma2, NumericMatrix Gtraining,
               NumericVector xtraining, NumericVector xtest) {
  int Ntraining = xtraining.length();
  int Ntest = xtest.length();
  double result = 0;
  for (int j = 0; j < Ntest; j ++){
    for (int i = 0; i < Ntraining; i ++){
      result = result + theta(i)*dphi_cpp(xtest(j), xtraining(i), sigma2);
    }
  }
  result = 2*result/Ntest;
  for (int i = 0; i < Ntraining; i ++){
    for (int j = 0; j < Ntraining; j ++){
      result = result + Gtraining(i,j)*theta(i)*theta(j);
    }
  }
  return result;
}
// Compute the loss for the second derivative given a fitted theta, a bandwidth sigma2,
// a training sample and a test sample
// [[Rcpp::export]]
double loss_d2(NumericVector theta, double sigma2, NumericMatrix Gtraining,
               NumericVector xtraining, NumericVector xtest) {
  int Ntraining = xtraining.length();
  int Ntest = xtest.length();
  double result = 0;
  for (int j = 0; j < Ntest; j ++){
    for (int i = 0; i < Ntraining; i ++){
      result = result + theta(i)*d2phi_cpp(xtest(j), xtraining(i), sigma2);
    }
  }
  result = -2*result/Ntest;
  for (int i = 0; i < Ntraining; i ++){
    for (int j = 0; j < Ntraining; j ++){
      result = result + Gtraining(i,j)*theta(i)*theta(j);
    }
  }
  return result;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variation with locally weigthed loss to focus in the neighborhood of a specific ystar.
//
// Gaussian kernel
// [[Rcpp::export]]
double K_cpp(double y, double ystar, double sigma2star) {
  return (exp(-pow(y-ystar,2)/(2*sigma2star))/sqrt(2*M_PI*sigma2star));
}
// first derivative Gaussian kernel
// [[Rcpp::export]]
double dK_cpp(double y, double ystar, double sigma2star) {
  return (-K_cpp(y,ystar,sigma2star)*(y-ystar)/sigma2star);
}
// second derivative Gaussian kernel
// [[Rcpp::export]]
double d2K_cpp(double y, double ystar, double sigma2star) {
  return ((-1/sigma2star + pow((y-ystar)/sigma2star,2))*K_cpp(y,ystar,sigma2star));
}
// Estimate derivatives via a constant
// [[Rcpp::export]]
double get_derivative_cpp(NumericVector ys, double ystar, double sigma2star, int order) {
  int N = ys.length();
  double result = 0;
  if (order==0){
    for (int j = 0; j < N; j ++){
      result = result + K_cpp(ys(j),ystar,sigma2star);
    }
  }
  if (order==1){
    for (int j = 0; j < N; j ++){
      result = result - dK_cpp(ys(j),ystar,sigma2star);
    }
  }
  if (order==2){
    for (int j = 0; j < N; j ++){
      result = result + d2K_cpp(ys(j),ystar,sigma2star);
    }
  }
  return (result/N);
}

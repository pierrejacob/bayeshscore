#ifndef LinearGaussian_H_
#define LinearGaussian_H_
#include <RcppEigen.h>
using namespace std;
using namespace Eigen;
using namespace Rcpp;

/**
 * Linear Gaussian models
 * Notation taken from Petris, Petrone and Campagnoli: Dynamical Linear Models with R.
 * (page 41)
 * Y_t = F_t Theta_t + v_t
 * Theta_t = G_t Theta_t-1 + w_t
 * v_t ~ N(0, V_t); w_t ~ N(0, W_t)
 * Theta_0 ~ N(m_0, C_0)
 *
 * Kalman recursion: Proposition 2.2, page 53 for filtering
 * Kalman smoothing: Proposition 2.4, page 61 for smoothing
 */


class LinearGaussian {
public:
  LinearGaussian();
  virtual ~LinearGaussian();
  MatrixXd F, G, V, W, C_0, invV, sdW;
  VectorXd m_0;
  void set_F(const NumericMatrix & F);
  void set_G(const NumericMatrix & F);
  void set_V(const NumericMatrix & V);
  void set_W(const NumericMatrix & W);
  void setLinearGaussianMatrices();
  void set_parameters(const List & parameters);
  void set_multivariate_parameters(double alpha, int d);
  
  void set_observations(const NumericMatrix & observations);
  //NumericMatrix observations;
  //NumericMatrix hidden_states;
  MatrixXd observations;
  int datalength;
  int XDimension;
  int YDimension;

};

#endif /* LinearGaussian_H_ */


#include "lineargaussian.h"

/**
 * Linear Gaussian models
 * Notation taken from Petris, Petrone and Campagnoli: Dynamical Linear Models with R.
 * (page 41)
 * for all t >= 1, Y_t = F_t Theta_t + v_t
 * for all t >= 1, Theta_t = G_t Theta_t-1 + w_t
 * where v_t ~ N(0, V_t); w_t ~ N(0, W_t)
 * initial condition
 * Theta_0 ~ N(m_0, C_0)
 * (there are no observations associated with the initial state)
 * Kalman recursion: Proposition 2.2, page 53 for filtering
 * Kalman smoothing: Proposition 2.4, page 61 for smoothing
 */

LinearGaussian::LinearGaussian() {
  XDimension = 1;
  YDimension = 1;
}

LinearGaussian::~LinearGaussian() {
}

void LinearGaussian::set_observations(const NumericMatrix & observations){
  Map<MatrixXd> map_obs = as<Map<MatrixXd> >(observations);
  this->observations = map_obs;
}

void LinearGaussian::setLinearGaussianMatrices(){
  this->F = MatrixXd::Identity(this->YDimension, this->XDimension);
  this->G = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->V = MatrixXd::Identity(this->YDimension, this->YDimension);
  this->invV = MatrixXd::Identity(this->YDimension, this->YDimension);
  this->W = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->sdW = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->m_0 = VectorXd::Zero(this->XDimension);
  this->C_0 = MatrixXd::Identity(this->XDimension, this->XDimension);
}

void LinearGaussian::set_F(const NumericMatrix & F){
  Map<MatrixXd> map_F = as<Map<MatrixXd> >(F);
  this->F = map_F;
  this->YDimension = this->F.rows();
  this->XDimension = this->F.cols();
}

void LinearGaussian::set_G(const NumericMatrix & G){
  Map<MatrixXd> map_G = as<Map<MatrixXd> >(G);
  this->G = map_G;
  this->XDimension = this->G.rows();

}

void LinearGaussian::set_V(const NumericMatrix & V){
  Map<MatrixXd> map_V = as<Map<MatrixXd> >(V);
  this->V = map_V;
}

void LinearGaussian::set_W(const NumericMatrix & W){
  Map<MatrixXd> map_W = as<Map<MatrixXd> >(W);
  this->W = map_W;
}

void LinearGaussian::set_parameters(const List & parameters){
  this->setLinearGaussianMatrices();
  this->G(0,0) = as<double>(parameters["rho"]);
  this->W(0,0) = as<double>(parameters["sigma"])  * as<double>(parameters["sigma"]);
  this->F(0,0) = as<double>(parameters["eta"]);
  this->V(0,0) = as<double>(parameters["tau"]) * as<double>(parameters["tau"]);
  // this->m_0(0) = 0.;
  // if (abs(this->G(0,0)) > pow(10.0, -10)){
  this->C_0(0, 0) = (1 - this->W(0,0)) / (this->G(0,0) * this->G(0,0));
  // } else {
    // this->C_0(0, 0) = 1.0;
  // }
}

void LinearGaussian::set_multivariate_parameters(double alpha, int d){
  this->XDimension = d;
  this->YDimension = d;
  this->setLinearGaussianMatrices();
  for (int i = 0; i < d; i++){
    for (int j = 0; j < d; j++){
      this->G(i,j) = std::pow(alpha, 1 + std::abs(i - j));
    }
  }
}


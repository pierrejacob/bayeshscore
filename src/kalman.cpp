#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <RcppEigen.h>
#include "kalman.h"
#include "lineargaussian.h"

using namespace std;
using namespace Eigen;

Kalman::Kalman() {
  XDimension = 1;
  YDimension = 1;
  this->F = MatrixXd::Identity(this->YDimension, this->XDimension);
  this->G = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->V = MatrixXd::Identity(this->YDimension, this->YDimension);
  this->invV = MatrixXd::Identity(this->YDimension, this->YDimension);
  this->W = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->sdW = MatrixXd::Identity(this->XDimension, this->XDimension);
  this->m_0 = VectorXd::Zero(this->XDimension);
  this->C_0 = MatrixXd::Identity(this->XDimension, this->XDimension);
}

Kalman::~Kalman() {
}

void Kalman::set_observations(const NumericMatrix & observations){
  Map<MatrixXd> map_obs = as<Map<MatrixXd> >(observations);
  this->observations = map_obs;
  this->nobs = this->observations.rows();
  incrLL = ArrayXd::Zero(this->nobs);
}

void Kalman::set_parameters(const List & parameters){
  this->G(0,0) = as<double>(parameters["rho"]);
  this->W(0,0) = as<double>(parameters["sigma"])  * as<double>(parameters["sigma"]);
  this->F(0,0) = as<double>(parameters["eta"]);
  this->V(0,0) = as<double>(parameters["tau"]) * as<double>(parameters["tau"]);
  this->C_0(0, 0) = this->W(0,0) / (1 - this->G(0,0) * this->G(0,0));
  this->m_0(0) = 0;
}

void Kalman::filtering(){
  /**
   * * KF computes the Kalman filter for the given observations and the given parameter
   * * Notation taken from Petris, Petrone and Campagnoli: Dynamical Linear Models with R.
   * * (page 41)
   * * Y_t = F_t Theta_t + v_t
   * * Theta_t = G_t Theta_t-1 + w_t
   * * v_t ~ N(0, V_t); w_t ~ N(0, W_t)
   * * Theta_0 ~ N(m_0, C_0)
   * *
   * * Kalman recursion: Proposition 2.2, page 53 for filtering
   * * Kalman smoothing: Proposition 2.4, page 61 for smoothing
   * */
  this->xFilterMeans.clear(); this->xFilterVariances.clear();
  histm.clear(); histC.clear(); hista.clear(); histR.clear();
  histf.clear(); histQ.clear();
  // VectorXd m(model->XDimension);
  // MatrixXd C(model->XDimension, model->XDimension);
  // VectorXd a(model->XDimension);
  // MatrixXd R(model->XDimension, model->XDimension);
  // VectorXd f(model->YDimension);
  // MatrixXd Q(model->YDimension, model->YDimension);

  // Initial sep
  m = this->m_0;
  C = this->C_0;
  this->xFilterMeans.push_back(this->m_0.array());
  this->xFilterVariances.push_back(C.array());
  histm.push_back(this->m_0); histC.push_back(C);


  for (unsigned int t = 0; t < this->nobs; t ++){
    this->filtering_step(t);
    hista.push_back(a); histR.push_back(R);
    histf.push_back(f); histQ.push_back(Q);
    histm.push_back(m); histC.push_back(C);
    this->xFilterMeans.push_back(m.array());
    this->xFilterVariances.push_back(C.array());

    // a = model->G * m;
    // R = model->G * C * model->G.transpose() + model->W;
    // hista.push_back(a); histR.push_back(R);
    // f = model->F * a;
    // Q = model->F * R * model->F.transpose() + model->V;
    // histf.push_back(f); histQ.push_back(Q);
    // e = this->model->observations.row(t);
    // e = e - f;
    // m = a + R * (model->F).transpose() * Q.inverse() * e;
    // C = R - R * (model->F).transpose() * Q.inverse() * (model->F) * R;
    // histm.push_back(m); histC.push_back(C);
    // this->xFilterMeans.push_back(m.array());
    // this->xFilterVariances.push_back(C.array());
    // incrLL(t) = - model->YDimension / 2. * log(2 * M_PI) - 0.5 * log(Q.determinant())
    //   - 0.5 * e.transpose() * Q.inverse() * e;
  }
}

void Kalman::first_step(){
  m = this->m_0;
  C = this->C_0;
}

void Kalman::filtering_step(int t){
  a = this->G * m;
  R = this->G * C * this->G.transpose() + this->W;
  f = this->F * a;
  Q = this->F * R * this->F.transpose() + this->V;
  VectorXd e = this->observations.row(t);
  e = e - f;
  m = a + R * (this->F).transpose() * Q.inverse() * e;
  C = R - R * (this->F).transpose() * Q.inverse() * (this->F) * R;
  incrLL(t) = - this->YDimension / 2. * log(2 * M_PI) - 0.5 * log(Q.determinant())
    - 0.5 * e.transpose() * Q.inverse() * e;
}

// void Kalman::smoothing(){
//   // this function assumes that the filtering function has been called already
//   hists.clear(); histS.clear();
//   // quantities introduced in Proposition 2.4 for smoothing
//   VectorXd s(model->XDimension);
//   MatrixXd S(model->XDimension, model->XDimension);
//   s = this->histm[this->nobs];
//   S = this->histC[this->nobs];
//   hists.push_back(s.array());
//   histS.push_back(S.array());
//   MatrixXd tmpinvR;
//   for (int t = (this->nobs - 1); t >= 0; t --){
//     tmpinvR = histR[t].inverse();
//     s = histm[t] + histC[t] * (model->G).transpose() * tmpinvR * (s - hista[t]);
//     S = histC[t] - histC[t] * (model->G).transpose() * tmpinvR * (histR[t] - S) * tmpinvR * model->G * histC[t];
//     hists.push_back(s.array());
//     histS.push_back(S.array());
//   }
//   vector<ArrayXd> reversehists = this->hists;
//   reverse(reversehists.begin(), reversehists.end());
//   xSmoothMeans = reversehists;
//   vector<ArrayXXd> reversehistS = this->histS;
//   reverse(reversehistS.begin(), reversehistS.end());
//   xSmoothVariances = reversehistS;
// }


ArrayXd Kalman::getIncrementalLL(void){
  return this->incrLL;
}

double Kalman::getLL(void){
  return this->getIncrementalLL().sum();
}

NumericVector Kalman::get_filtering_mean(int time_step){
  return wrap(this->xFilterMeans[time_step]);
}

NumericMatrix Kalman::get_filtering_variance(int time_step){
  return wrap(this->xFilterVariances[time_step]);
}

// NumericVector Kalman::get_smoothing_mean(int time_step){
//   return wrap(this->xSmoothMeans[time_step]);
// }
//
// NumericMatrix Kalman::get_smoothing_variance(int time_step){
//   return wrap(this->xSmoothVariances[time_step]);
// }

// NumericMatrix get_F(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->F);}
// NumericMatrix get_G(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->G);}
// NumericMatrix get_V(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->V);}
// NumericMatrix get_W(LinearGaussian* LinearGaussian){return wrap(LinearGaussian->W);}
// int get_nobs(LinearGaussian* LinearGaussian){return LinearGaussian->observations.rows();}
NumericVector get_incremental_ll(Kalman* kalman){return wrap(kalman->getIncrementalLL());}


// RCPP_EXPOSED_CLASS(LinearGaussian)
RCPP_EXPOSED_CLASS(Kalman)

  RCPP_MODULE(kalman_mod) {
    // class_<LinearGaussian>( "LinearGaussian" )
    // .constructor()
    // .method( "set_F", &LinearGaussian::set_F)
    // .method( "set_G", &LinearGaussian::set_G)
    // .method( "set_V", &LinearGaussian::set_V)
    // .method( "set_W", &LinearGaussian::set_W)
    // .method( "get_F", &get_F)
    // .method( "get_G", &get_G)
    // .method( "get_V", &get_V)
    // .method( "get_W", &get_W)
    // .method( "setLinearGaussianMatrices", &LinearGaussian::setLinearGaussianMatrices)
    // .method( "set_parameters", &LinearGaussian::set_parameters)
    // .method( "set_multivariate_parameters", &LinearGaussian::set_multivariate_parameters)
    // .method( "set_observations", &LinearGaussian::set_observations)
    // .method( "get_nobs", &get_nobs)
    // ;

    class_<Kalman>( "Kalman" )
      .constructor()
      // .method( "setLinearGaussian", &Kalman::setLinearGaussian)
      .method( "filtering", &Kalman::filtering)
      .method( "first_step", &Kalman::first_step)
      .method( "filtering_step", &Kalman::filtering_step)
      .method( "set_parameters", &Kalman::set_parameters)
      .method( "set_observations", &Kalman::set_observations)
      // .method( "smoothing", &Kalman::smoothing)
      .method( "get_filtering_mean", &Kalman::get_filtering_mean)
      .method( "get_filtering_variance", &Kalman::get_filtering_variance)
      // .method( "get_smoothing_mean", &Kalman::get_smoothing_mean)
      // .method( "get_smoothing_variance", &Kalman::get_smoothing_variance)
      .method( "getLL", &Kalman::getLL)
      .method( "get_incremental_ll", &get_incremental_ll)
      .field( "nobs", &Kalman::nobs)
    ;

  }

#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// One step of Kalman filter
// [[Rcpp::export]]
List KF_assimilate_one_cpp(NumericVector Yt, int t,
                           NumericVector initial_mean, NumericMatrix initial_var,
                           NumericMatrix phi, NumericMatrix psi,
                           NumericMatrix sigmaV2, NumericMatrix sigmaW2,
                           List KF_current){
  // Convert all the R objects in Eigen class objects for easier matrix manipulation
  const Eigen::Map<Eigen::VectorXd> Yt_(as<Eigen::Map<Eigen::VectorXd> >(Yt));
  const Eigen::Map<Eigen::VectorXd> initial_mean_(as<Eigen::Map<Eigen::VectorXd> >(initial_mean));
  const Eigen::Map<Eigen::MatrixXd> initial_var_(as<Eigen::Map<Eigen::MatrixXd> >(initial_var));
  const Eigen::Map<Eigen::MatrixXd> phi_(as<Eigen::Map<Eigen::MatrixXd> >(phi));
  const Eigen::Map<Eigen::MatrixXd> psi_(as<Eigen::Map<Eigen::MatrixXd> >(psi));
  const Eigen::Map<Eigen::MatrixXd> sigmaV2_(as<Eigen::Map<Eigen::MatrixXd> >(sigmaV2));
  const Eigen::Map<Eigen::MatrixXd> sigmaW2_(as<Eigen::Map<Eigen::MatrixXd> >(sigmaW2));
  List KF_updated = KF_current;
  // Create containers for the results
  int dimX = initial_var.nrow();
  int dimY = sigmaV2.nrow();
  Eigen::VectorXd muX_t_t_1(dimX); // contains the means of X_t given Y_1,...,Y_t-1
  Eigen::VectorXd muX_t_t(dimX); // contains the means of X_t given Y_1,...,Y_t
  Eigen::VectorXd muY_t_t_1(dimY); // contains the means of Y_t given Y_1,...,Y_t-1
  Eigen::MatrixXd PX_t_t_1(dimX,dimX); // contains the variances of X_t given Y_1,...,Y_t-1
  Eigen::MatrixXd PX_t_t(dimX,dimX); // contains the variances of X_t given Y_1,...,Y_t
  Eigen::MatrixXd PY_t_t_1(dimY,dimY); // contains the variances of Y_t given Y_1,...,Y_t-1
  Eigen::MatrixXd Kt; // Kalman gain matrix
  if (t==0) {
    muX_t_t_1 = initial_mean_;
    PX_t_t_1 = initial_var_;
    Kt = PX_t_t_1*(psi_.transpose())*((psi_*PX_t_t_1*(psi_.transpose())+sigmaV2_).inverse());
    muX_t_t = muX_t_t_1 + Kt * (Yt_ - psi_ * muX_t_t_1);
    PX_t_t = (Eigen::MatrixXd::Identity(dimX,dimX) - Kt*psi_) * PX_t_t_1;
    muY_t_t_1 = psi_ * muX_t_t_1;
    PY_t_t_1 = psi_ * PX_t_t_1 * (psi_.transpose()) + sigmaV2_;
  } else {
    // extract and convert previous parameters
    List KF = KF_current[t-1];
    const Eigen::Map<Eigen::VectorXd> muX_t_t_(as<Eigen::Map<Eigen::VectorXd> >(KF["muX_t_t"]));
    const Eigen::Map<Eigen::MatrixXd> PX_t_t_(as<Eigen::Map<Eigen::MatrixXd> >(KF["PX_t_t"]));
    // update parameters
    muX_t_t_1 = phi_ * muX_t_t_;
    PX_t_t_1 = phi_ * PX_t_t_ * (phi_.transpose()) + sigmaW2_;
    Kt = PX_t_t_1*(psi_.transpose())*((psi_*PX_t_t_1*(psi_.transpose())+sigmaV2_).inverse());
    muX_t_t = muX_t_t_1 + Kt * (Yt_ - psi_ * muX_t_t_1);
    PX_t_t = (Eigen::MatrixXd::Identity(dimX,dimX) - Kt*psi_) * PX_t_t_1;
    muY_t_t_1 = psi_ * muX_t_t_1;
    PY_t_t_1 = psi_ * PX_t_t_1 * (psi_.transpose()) + sigmaV2_;
  }
  KF_updated[t] = List::create(Named("muX_t_t_1") = muX_t_t_1,
                               Named("PX_t_t_1") = PX_t_t_1,
                               Named("muX_t_t") = muX_t_t,
                               Named("PX_t_t") = PX_t_t,
                               Named("muY_t_t_1") = muY_t_t_1,
                               Named("PY_t_t_1") = PY_t_t_1);
  return KF_updated;
}
// Kalman filter assimilating all the observations in Y (dimY by nobservations)
// [[Rcpp::export]]
List KF_filtering_cpp(NumericMatrix Y,
                      NumericVector initial_mean, NumericMatrix initial_var,
                      NumericMatrix phi, NumericMatrix psi,
                      NumericMatrix sigmaV2, NumericMatrix sigmaW2){
  int nobservations = Y.ncol();
  List KF(nobservations);
  for (int t = 0; t < nobservations; t ++){
    KF = KF_assimilate_one_cpp(Y(_,t), t, initial_mean, initial_var, phi, psi, sigmaV2, sigmaW2, KF);
  }
  return KF;
}

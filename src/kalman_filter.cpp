#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_;
  P_ = F_ * P_ * F_.transpose()+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd I = MatrixXd::Identity(F_.rows(),F_.cols());

  VectorXd y = z - H_ * x_;

  MatrixXd PHt = P_ * H_.transpose(); 
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  MatrixXd Hj = tools.CalculateJacobian(x_);
  MatrixXd I = MatrixXd::Identity(F_.rows(),F_.cols());

  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  double rho = sqrt(px*px+py*py);
  double phi = atan2(py/px);
  double rho_dot = (px*vx+py*vy)/rho;
  VectorXd z_pred = VectorXd(3);
  z_pred<<rho,phi,rho_dot;

  VectorXd y = z - z_pred;

  while (y(1)>M_PI/2){
    y(1) = y(1) - M_PI;
  }
  while (y(1)<-M_PI){
    y(1) = y(1) + M_PI;
  }

  MatrixXd PHt = P_ * Hj.transpose(); 
  MatrixXd S = Hj * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * Hj) * P_;
}

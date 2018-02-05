#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_<< 1,0,0,0,
            0,1,0,0;
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 0.0225, 0, 0, 0,
        0, 0.0225, 0, 0,
        0, 0, 0.0225, 0,
        0, 0, 0, 0.0225;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  long long t = measurement_pack.timestamp_ - previous_timestamp_;
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double rho_dot = measurement_pack.raw_measurements_(2);
      double px = rho * sin(phi);
      double py = rho * cos(phi);
      double vx = rho_dot/cos(phi) * sin(M_PI/2 - 2*phi);
      double vy = rho_dot/cos(phi) * cos(M_PI/2 - 2*phi);

      ekf_.x_ << px,py,vx,vy;

      


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double px = measurement_pack.raw_measurements_(0);
      double py = measurement_pack.raw_measurements_(1);
      ekf_.x_ << px,py,1,1;

    }
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1,0,t,0,
              0,1,0,t,
              0,0,1,0,
              0,0,0,1;
    ekf_.H_ = H_laser_;
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double noise_ax = 9;
  double noise_ay = 9;
  double c1 = pow(t,4)/4;
  double c2 = pow(t,3)/2;
  double c3 = pow(t,2);
  double x_sq = pow(noise_ax,2);
  double y_sq = pow(noise_ay,2);
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << c1*x_sq, 0, c2*x_sq, 0,
              0, c1*y_sq, 0, c2*y_sq,
              c2*x_sq, 0, c3*x_sq, 0,
              0, c2*y_sq, 0, c3*y_sq;

  ekf_.Predict();

  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    // Hj_ = toos.CalculateJacobian(measurement_pack.raw_measurements_); 
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

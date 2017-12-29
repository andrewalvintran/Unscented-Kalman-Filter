#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  
  // augmented mean state
  x_aug << x_, 0.0, 0.0;

  // augmented covariance
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // creating augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(n_aug_ + lambda_)*L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(n_aug_ + lambda_)*L.col(i);
  }

  double delta_t_sq = delta_t * delta_t;

  // predict sigma points
  for (int i = 0; i < n_aug_*2 + 1; i++) {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);
    
    double px_p, py_p;
    
    if (fabs(yaw_dot) <= 0.001) {
        px_p = p_x + v * delta_t_sq * cos(yaw);
        py_p = p_y + v * delta_t_sq * sin(yaw);
    } else {
        px_p = p_x + v/yaw_dot * (sin(yaw + yaw_dot*delta_t_sq) - sin(yaw));
        py_p = p_y + v/yaw_dot * (-cos(yaw + yaw_dot*delta_t_sq) + cos(yaw));
    }
    
    Xsig_pred_(0, i) = px_p + 0.5*(delta_t_sq)*cos(yaw)*nu_a;
    Xsig_pred_(1, i) = py_p + 0.5*(delta_t_sq)*sin(yaw)*nu_a;
    Xsig_pred_(2, i) = v + delta_t_sq * nu_a;
    Xsig_pred_(3, i) = yaw + yaw_dot*delta_t_sq + 0.5*(delta_t_sq * delta_t_sq) * nu_yawdd;
    Xsig_pred_(4, i) = yaw_dot + delta_t_sq*nu_yawdd;
  }

  // predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0.0, 0.0, 0.0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      double p_x = rho * cos(phi);
      double p_y = rho * sin(phi);

      double v_x = rho_dot * cos(phi);
      double v_y = rho_dot * sin(phi);
      double v = sqrt(v_x*v_x + v_y*v_y); 
      x_ << p_x, p_y, v, 0.0, 0.0;
    }

    //set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2*n_aug_ + 1; i++) {
        weights_(i) = 0.5 / (lambda_ + n_aug_);
    }

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
  } else {
    const double NUM_SEC_IN_MICRO = 1000000.0;
    double delta_t = (meas_package.timestamp_ - time_us_) / NUM_SEC_IN_MICRO;
    time_us_ = meas_package.timestamp_;
    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
  }
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
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    } else {
        px_p = p_x + v/yaw_dot * (sin(yaw + yaw_dot*delta_t) - sin(yaw));
        py_p = p_y + v/yaw_dot * (-cos(yaw + yaw_dot*delta_t) + cos(yaw));
    }

    const double delta_t_sq = delta_t * delta_t;
    Xsig_pred_(0, i) = px_p + 0.5*(delta_t_sq)*cos(yaw)*nu_a;
    Xsig_pred_(1, i) = py_p + 0.5*(delta_t_sq)*sin(yaw)*nu_a;
    Xsig_pred_(2, i) = v + delta_t * nu_a;
    Xsig_pred_(3, i) = yaw + yaw_dot*delta_t + 0.5*(delta_t_sq) * nu_yawdd;
    Xsig_pred_(4, i) = yaw_dot + delta_t*nu_yawdd;
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

void UKF::Update(MeasurementPackage meas_package, MatrixXd Zsig) {
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);

  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while(z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;
    while(z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_, n_z_);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    R << std_radr_*std_radr_, 0.0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
  } else {
    R << std_laspx_*std_laspx_, 0.0,
         0.0, std_laspy_*std_laspy_;
  }
  S += R;

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      while(z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;
      while(z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;

      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;

  VectorXd z_diff = z - z_pred;
  while(z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;
  while(z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  n_z_ = 2;
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z_, n_aug_*2 + 1);
  Update(meas_package, Zsig);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  n_z_ = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_ + 1);

  for (int i = 0; i < n_aug_*2 + 1; i++) {
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);
      
      double v1 = p_x*cos(yaw)*v;
      double v2 = p_y*sin(yaw)*v;
      
      Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig(1, i) = atan2(p_y, p_x);
      Zsig(2, i) = (v1 + v2) / sqrt(p_x*p_x + p_y*p_y);
  }

  Update(meas_package, Zsig);
}

#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define Eps 0.0001 // to check if values are too small

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initialized state information
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // n_x is the number of states
  n_x_ = 5;

  // num_aug is the number of augmented state dimensions (two more dimensions than states)
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu psi] in SI units and rad
  x_aug_ = VectorXd(n_aug_);

  X_sigma_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // initial aug state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Process noise standard variance longitudinal acceleration in m/s^2
  var_a_ = 1;

  // Process noise standard variance yaw acceleration in rad/s^2
  var_yawdd_ = 0.5*0.5;

  // Laser measurement noise variance deviation position1 in m
  var_laspx_ = 0.15*0.15;

  // Laser measurement noise variance deviation position2 in m
  var_laspy_ = 0.15*0.15;

  // Radar measurement noise variance deviation radius in m
  var_radr_ = 0.3*0.3;

  // Radar measurement noise variance deviation angle in rad
  var_radphi_ = 0.03*0.03;

  // Radar measurement noise standard deviation radius change in m/s
  var_radrd_ = 0.3*0.3;

  // lambda is the sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // time_us is the parameter that keeps track of time to generate delta t
  time_us_ = 0.0;

  // initializing matrices
  R_lidar_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);

  R_radar_ << var_radr_, 0, 0,
              0, var_radphi_, 0,
              0, 0,var_radrd_;
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << var_laspx_,0,
              0,var_laspy_;

  NIS_radar_ = 0.0;
  NIS_lidar_ = 0.0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF for first measurement: " << endl;


    // Convert radar from polar to cartesian coordinates and initialize state.
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // polar to cartesian
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v  = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0.0, 0.0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0, 0.0;
    }

    // initialize the augmented state
    x_aug_.head(5) << x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;

    // initial uncertainty, 1 for the position x/y and 1000 for the speed
    P_ <<   1  ,  0  ,  0  ,  0  ,  0  ,
            0  ,  1  ,  0  ,  0  ,  0  ,
            0  ,  0  ,  1  ,  0  ,  0  ,
            0  ,  0  ,  0  ,  1  ,  0  ,
            0  ,  0  ,  0  ,  0  ,  1  ;

    // initial uncertainty for augmented state, fill top left n_x_ with P_ and use
    // longitudinal and yaw accelartion noise for the other two dimensions
    P_aug_.fill(0);
    P_aug_.topLeftCorner(n_x_,n_x_) = P_;
    P_aug_(5,5) = var_a_;
    P_aug_(6,6) = var_yawdd_;

    // set UKF weights
    weights_(0) = lambda_/(lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; i++){
      weights_(i) = 0.5/(n_aug_ + lambda_);
    }

    // Check to see if the initial values of px or py are too small
    if (fabs(x_(0)) < Eps and fabs(x_(1)) < Eps) {
      x_(0) = Eps;
      x_(1) = Eps;
    }

    time_us_ = measurement_pack.timestamp_;

    cout << "Initilization finished" << endl;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // elapsed time, converted to seconds
  float delta_t = (measurement_pack.timestamp_ - time_us_)/1000000.0;
  time_us_ = measurement_pack.timestamp_;

  Prediction(delta_t);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  Tools tools;
  // Here there is no need for augmentation, since the noise does not enter the nonlinear
  // function.

  // TODO: try to generate new sigma points instead of just using the ones we got from
  // prediction

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar
    UpdateRadar(measurement_pack);

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    // Lidar
    UpdateLidar(measurement_pack);
  }

  // print the normalized innovation squared
  cout << "NIS_radar_ = " << NIS_radar_ << endl;
  cout << "NIS_lidar_" << NIS_lidar_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  GenerateSigmaPoints();

  PredictSigmaPoints(delta_t);

  PredictMeanAndCovariance();

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z_ = 2; // The lidar has 2 dimensions, px and py

  // Transform sigma points into measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z_, 2 * n_aug_ + 1); // lidar is linear, no extra steps

  UpdateUKF( meas_package, Zsig, n_z_);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //set measurement dimension, rho, phi, and rho_dot
  int n_z_ = 3; // The radar has 3 dimensions, rho, phi, and rho_dot

  // Transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1); // sigma points

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  UpdateUKF( meas_package, Zsig, n_z_);

}

void UKF::GenerateSigmaPoints() {

  x_aug_.head(n_x_) = x_;

  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = var_a_;
  P_aug_(6,6) = var_yawdd_;

  //create square root matrix
  MatrixXd P_aug_inv_ = P_aug_.llt().matrixL();

  X_sigma_aug_.col(0) = x_aug_;

  double scale = sqrt(lambda_+n_aug_); // Save some computations

  for (int i = 0; i< n_aug_; i++)
  {
    X_sigma_aug_.col(i + 1) = x_aug_ + scale * P_aug_inv_.col(i);
    X_sigma_aug_.col(i + 1 + n_aug_) = x_aug_ - scale * P_aug_inv_.col(i);
  }

}

void UKF::PredictSigmaPoints(double delta_t) {

  for (int i = 0; i < 2* n_aug_ + 1; i++)
  {
    ProcessModel(delta_t, i);
  }

}

void UKF::PredictMeanAndCovariance(){
  Tools tools;
  //predicted state mean
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools.NormalizeAngle( (double) x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::ProcessModel(double delta_t, int i) {
  // This function implements the Constant Turn Rate and Velocity (CTRV) Model

  // states
  double px = X_sigma_aug_(0,i);
  double py = X_sigma_aug_(1,i);
  double v = X_sigma_aug_(2,i);
  double yaw = X_sigma_aug_(3,i);
  double yaw_rate = X_sigma_aug_(4,i);
  // noise
  double nu_a = X_sigma_aug_(5,i);
  double nu_yaw = X_sigma_aug_(6,i);

  // reduce amount of calculations
  double yaw_delta_t = yaw + yaw_rate * delta_t;
  double sin_yaw = sin(yaw);
  double cos_yaw = cos(yaw);
  double delta_t2 = delta_t * delta_t;

  if (fabs(yaw_rate) > Eps) {
    Xsig_pred_(0,i) = px + (v / yaw_rate)*(sin(yaw_delta_t) - sin_yaw) + 0.5*delta_t2*cos_yaw*nu_a;
    Xsig_pred_(1,i) = py + (v / yaw_rate)*(-cos(yaw_delta_t) + cos_yaw) + 0.5*delta_t2*sin_yaw*nu_a;
    Xsig_pred_(2,i) = v + delta_t*nu_a;
    Xsig_pred_(3,i) = yaw_delta_t + 0.5*delta_t2*nu_yaw;
    Xsig_pred_(4,i) = yaw_rate + delta_t*nu_yaw;
  } else{
    Xsig_pred_(0,i) = px + v*cos_yaw*delta_t + 0.5*delta_t2*cos_yaw*nu_a;
    Xsig_pred_(1,i) = py + v*sin_yaw*delta_t + 0.5*delta_t2*sin_yaw*nu_a;
    Xsig_pred_(2,i) = v + delta_t*nu_a;
    Xsig_pred_(3,i) = yaw + 0.5*delta_t2*nu_yaw;
    Xsig_pred_(4,i) = yaw_rate + delta_t*nu_yaw;
  }

}

void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z_){
  // Update the measurement
  Tools tools;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
//    z_pred.fill(0.0);
//    for (int i=0; i < 2*n_aug_+1; i++) {
//        z_pred = z_pred + weights_(i) * Zsig.col(i);
//    }
  z_pred  = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    z_diff(1) = tools.NormalizeAngle( (double) z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_, n_z_);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    R = R_lidar_;
  }

  S = S + R;


  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools.NormalizeAngle( (double) z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools.NormalizeAngle( (double) x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_meas = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z_meas - z_pred;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    //angle normalization
    z_diff(1) = tools.NormalizeAngle( (double) z_diff(1));
  }


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    NIS_radar_ = z_meas.transpose() * S.inverse() * z_meas;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    NIS_lidar_ = z_meas.transpose() * S.inverse() * z_meas;
  }
}



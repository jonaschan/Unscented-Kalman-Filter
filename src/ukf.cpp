#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() 
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  //***********************************************************************
  // Set initial example covariance matrix (Values obtained from lectures)
  P_ = MatrixXd(5, 5);
  
  // Here we initialise the covariance matrix using the identity matrix to
  // allow these initialised values to be used across other datasets.
  P_ << 1,		0,		0,		0,		 0,
		0,		1,		0,		0,		 0,
	    0,		0,		1,		0,		 0,
		0,		0,		0,		1,		 0,
		0,		0,		0,		0,		 1;
  //***********************************************************************

  //***********************************************************************
  //                  Initialise variables with values	                 //
  //***********************************************************************	  
  
  //This value is tuned to give the expected RMSE values
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  //This value is tuned to give the expected RMSE values
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Initiate matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
 
  // Weights of sigma points
  weights_ = VectorXd(n_sig_);
  
  //Set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i = 1; i < n_sig_; i++)
  {
	  weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

	// In this function we, perform the following checks:
	//	1.	Check if the states have been initialised, if yes, predict and update
	//		the using the new values, if no, initialise the values.
	//	2.	Check if the measurement type; LIDAR or RADAR

	// Initialise the state vector if the is_initialised flag is false
	if (!is_initialized_)
	{
		time_us_ = meas_package.timestamp_;

		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			float px = meas_package.raw_measurements_(0);
			float py = meas_package.raw_measurements_(1);

			x_ << px, py, 0, 0, 0;
		}

		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			**/
			float rho = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float rho_dot = meas_package.raw_measurements_(2);

			/**
			Here we perform the calculations to convert polar coordinates to
			cartesian coordinates using the Pythogoras theorem
			**/
			float px = rho * cos(phi);
			float py = rho * sin(phi);
			
			x_ << px, py, 0, 0, 0;
		}

		is_initialized_ = true;
	}

	// Here we predict and update the values if the is_initialise flag is true
	// (state vector has been initialised)
	else
	{
		double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;
		
		Prediction(delta_t);

		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			UpdateLidar(meas_package);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
		{
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
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  **/

	// Create augmented mean vector
	VectorXd x_aug_(7);
	x_aug_.fill(0);
	x_aug_.head(n_x_) = x_;
	
	// Create augmented state covariance
	MatrixXd P_aug_(n_aug_, n_aug_);
	P_aug_.fill(0);
	P_aug_.topLeftCorner(n_x_, n_x_) = P_;
	P_aug_(5, 5) = pow(std_a_, 2);
	P_aug_(6, 6) = pow(std_yawdd_, 2);

	// Create square root matrix
	MatrixXd L_= P_aug_.llt().matrixL();

	// Create augmented sigma points
	MatrixXd Xsig_aug_(n_aug_, n_sig_);
	
	Xsig_aug_.col(0) = x_aug_;
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug_.col(i + 1)		  = x_aug_ + sqrt(lambda_ + n_aug_) * L_.col(i);
		Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L_.col(i);
	}

	//Predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x		= Xsig_aug_(0, i);
		double p_y		= Xsig_aug_(1, i);
		double v		= Xsig_aug_(2, i);
		double yaw		= Xsig_aug_(3, i);
		double yawd		= Xsig_aug_(4, i);
		double nu_a		= Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001)
		{
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		else
		{
			px_p = p_x + v * delta_t * cos(yaw);
			py_p = p_y + v * delta_t * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		// Add the noise component
		px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		// Write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	// Predict state mean
	x_ = Xsig_pred_ * weights_;

	// Predicted state covariance matrix
	P_.fill(0.0);
		
	// Iterate over sigma points
	for (int i = 0; i < n_sig_; i++)
	{  	// state difference
		VectorXd x_difff = Xsig_pred_.col(i) - x_;
			
		x_difff(3) = x_difff(3) - ceil((x_difff(3) - M_PI) / (2.*M_PI))*2.*M_PI;
		P_ = P_ + weights_(i) * x_difff * x_difff.transpose();
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

	MatrixXd Zsig_ = MatrixXd(2, 15);
	VectorXd z_pred_ = VectorXd(2);

	Zsig_.row(0) = Xsig_pred_.row(0);
	Zsig_.row(1) = Xsig_pred_.row(1);

	// Get measurement states
	z_pred_ = Zsig_ * weights_;

	// Create measurement covariance matrix, S
	// Covariance matrix is a square matrix of which its dimensions
	// is equal to the number of the states. Since LIDAR has two states,
	// px and py, hence the dimension of S is 2 X 2.
	MatrixXd S_ = MatrixXd(2, 2);
	S_.fill(0);

	// Create cross correlation matrix, Tc
	MatrixXd Tc_ = MatrixXd(5, 2);
	Tc_.fill(0);

	// Calculate the cross correlation matrix
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;

		// Calculate the measurement covariance matrix
		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//	Calculate the cross correlation matrix
		Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
	}

	// Calculate the LIDAR measurement noise covariance matrix
	MatrixXd noise_laser(2, 2);

	noise_laser <<	std_laspx_ * std_laspx_, 0,
					0, std_laspy_ * std_laspy_;

	// Add measurement noise covariance matrix
	S_ = S_ + noise_laser;

	// Kalman gain, K
	MatrixXd K = Tc_ * (S_.inverse());

	// Residual
	VectorXd z_laser_residual_ = meas_package.raw_measurements_ - z_pred_;

	// Update state mean
	x_ = x_ + K * z_laser_residual_;

	// Update NIS
	NIS_laser_ = z_laser_residual_.transpose() * S_.inverse() * z_laser_residual_;

	// Update covariance matrix
	P_ = P_ - K * S_ * K.transpose();
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	MatrixXd Zsig_Radar_ = MatrixXd(3, n_sig_);
	VectorXd z_pred_ = VectorXd(3);

	for (int i = 0; i < n_sig_; i++) 
	{
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
		double yawd = Xsig_pred_(4, i);

		// Calculate the directional velocities in the x and y 
		// vector using Pythogoras Theorem
		double vx = cos(yaw) * v;
		double vy = sin(yaw) * v;

		// Calculate the predicted states
		Zsig_Radar_(0, i) = sqrt(px * px + py * py);
		Zsig_Radar_(1, i) = atan2(py, px);
		
		// We are still calculating the predicted state, 
		// but at this point we are checking for division
		// by zero
		if (sqrt(px * px + py * py) > 0.001)
		{
			Zsig_Radar_(2, i) = (px * vx + py * vy) / sqrt(px * px + py * py);
		}
		else 
		{
			Zsig_Radar_(2, i) = 0.0;
		}
	}

	// Get radar measurement states
	z_pred_ = Zsig_Radar_ * weights_;

	// Create measurement covariance matrix, S
	// Covariance matrix is a square matrix of which its dimensions
	// is equal to the number of the states. Since RADAR has three states,
	// rho, phi and rho dot, hence the dimension of S is 3 X 3.
	MatrixXd S = MatrixXd(3, 3);
	S.fill(0);

	// Create cross relation matrix
	MatrixXd Tc = MatrixXd(5, 3);
	Tc.fill(0);

	// Calculate the correlation matrix
	for (int i = 0; i<n_sig_; i++) 
	{
		// residual
		VectorXd z_diff = Zsig_Radar_.col(i) - z_pred_;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// Calculate the measurement correlation matrix
		S = S + weights_(i) * z_diff * z_diff.transpose();

		// Calculate the cross correlation matrix
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	// Calculate the RADAR measurement noise covariance matrix
	MatrixXd noise_radar(3, 3);
	
	noise_radar <<	std_radr_ * std_radr_,		0,	0,
					0, std_radphi_ * std_radphi_,	0,
					0,	0,	  std_radrd_ * std_radrd_;

	// Add measurement noise covariance matrix
	S = S + noise_radar;

	// Kalman gain, K
	MatrixXd K = Tc * S.inverse();

	// Residual
	VectorXd z_radar_residual = meas_package.raw_measurements_ - z_pred_;

	// Update state mean
	x_ = x_ + K * z_radar_residual;

	// Update NIS
	NIS_radar_ = z_radar_residual.transpose() * S.inverse() * z_radar_residual;

	// Update covariance matrix
	P_ = P_ - K * S * K.transpose();
}

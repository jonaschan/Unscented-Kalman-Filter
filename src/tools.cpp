#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  /**
  TODO:
    * Calculate the RMSE here.
  */

  //Initialise an empty vector for the calculation of RMSE values
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  //Check the validity of the inputs:
  // 1. Check the estimations vector size should not be zero
  // 2. Check the estimations vector size should be equal to the size of the ground_truth vector size

  if (estimations.size() == 0 || estimations.size() != ground_truth.size())
  {
	  cout << "Estimations vector is zero or estimations vector size not equal to the size of the ground truth vector" << endl;

	  return rmse;
  }

  //Accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); i++)
  {
	  VectorXd error = estimations[i] - ground_truth[i];

	  //Get the square of the error
	  error = error.array() * error.array();

	  //Sum the squared error
	  rmse += error;
  }

  //Calculate the mean squared error
  rmse = rmse / estimations.size();

  //Calculate the square root mean squared error
  rmse = rmse.array().sqrt();

  return rmse;
}
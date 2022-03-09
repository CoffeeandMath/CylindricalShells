/*
 * referenceconfiguration.h
 *
 *  Created on: Jul 14, 2021
 *      Author: kevin
 */

#ifndef SRC_REFERENCE_CONFIGURATION_H_
#define SRC_REFERENCE_CONFIGURATION_H_




#include <boost/math/quadrature/trapezoidal.hpp>
#include <iostream>
using boost::math::quadrature::trapezoidal;
#include <math.h>
#include <Eigen/Dense>

class Reference_Configuration {
public:
	Reference_Configuration();
	virtual ~Reference_Configuration();
	void set_point(double);
	void set_deformation_param(double);
	Eigen::Matrix2d get_Ac();
	Eigen::Matrix2d get_Bc();
	Eigen::Matrix2d get_AC();

	void set_R0(double);
	double get_R();
	double get_phi();
	double get_dphi();
	Eigen::Matrix2d get_P();



private:
	void calc_covariants();

	double phifun(double);
	double dphifun(double);
	double defmag = 0.0;

	double S_Point = 0.0;
	double Rval = 0.0;
	double dRdSval = 0.0;
	double ddRddSval = 0.0;
	double phival = 0.0;
	double dphidSval = 0.0;

	double R0 = 2.0;



	Eigen::Matrix2d Cov;
	Eigen::Matrix2d Form2;

	Eigen::Matrix2d AC;


	Eigen::Matrix2d Proj;
};

#endif /* SRC_REFERENCE_CONFIGURATION_H_ */

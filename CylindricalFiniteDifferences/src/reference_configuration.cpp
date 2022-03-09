/*
 * referenceconfiguration.cpp
 *
 *  Created on: Jul 14, 2021
 *      Author: kevin
 */

#include "reference_configuration.h"

Reference_Configuration::Reference_Configuration() {
	// TODO Auto-generated constructor stub
	Cov.setZero();
	Form2.setZero();
	Proj.setZero();
}

Reference_Configuration::~Reference_Configuration() {
	// TODO Auto-generated destructor stub
}




void Reference_Configuration::set_point(double S){
	S_Point = S;

	calc_covariants();


}

void Reference_Configuration::calc_covariants(){

	auto sinphifun = [this](double s){
		return sin(phifun(s));
	};




	double dR = 0.0;
	if (S_Point>0){
		double tol = 1e-6;
		int max_refinements = 20;
		dR = trapezoidal(sinphifun,0.0,S_Point,tol, max_refinements);
		//std::cout << dR << std::endl;
	}
	dR = 0.0;
	Rval = R0 + dR;
	Cov.setZero();
	Cov(0,0) = 1.;
	//Cov(1,1) = pow(Rval,2.0);
	Cov(1,1) = 1.0;
	AC = Cov.inverse();

	Proj.setZero();
	Proj(0,0) = sqrt(Cov(0,0));
	Proj(0,1) = Cov(0,1)/sqrt(Cov(0,0));
	Proj(1,1) = sqrt(Cov(1,1) - pow(Cov(0,1),2)/Cov(0,0));
	Proj(1,0) = 0.;

	Form2.setZero();
	Form2(0,0) = -dphifun(S_Point);
	Form2(1,1) = cos(phifun(S_Point))/Rval;

	Form2(0,0) = 0.0;
	Form2(1,1) = 1.0;



}

double Reference_Configuration::phifun(double s) {
	double offset = defmag;
	double offsetlim = 0.5;
	if (offset > offsetlim) {
		offset = offsetlim;
	}

	return 0.5*defmag*s*s + offset;
	//return -defmag * s + offset;
}

double Reference_Configuration::dphifun(double s) {
	return defmag * s;
	//return -defmag;
}

void Reference_Configuration::set_deformation_param(double lambda){
	defmag = lambda;
}

void Reference_Configuration::set_R0(double Rtemp){
	R0 = Rtemp;
}

Eigen::Matrix2d Reference_Configuration::get_Ac(){

	return Cov;
}

Eigen::Matrix2d Reference_Configuration::get_Bc(){

	return Form2;
}

Eigen::Matrix2d Reference_Configuration::get_AC(){
	return AC;
}


double Reference_Configuration::get_R(){
	return Rval;
}

double Reference_Configuration::get_phi(){
	return phifun(S_Point);
}

Eigen::Matrix2d Reference_Configuration::get_P(){
	return Proj;
}

double Reference_Configuration::get_dphi(){
	return dphifun(S_Point);
}

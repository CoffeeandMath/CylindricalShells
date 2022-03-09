#ifndef MATERIAL_CLASS_H_
#define MATERIAL_CLASS_H_

#include <fstream>
#include <iostream>
#include <string>
#include "reference_configuration.h"
#include "unsupported/Eigen/CXX11/Tensor"




class Material_Class
{
public:
	Material_Class();
	Material_Class(const double,const double,const Eigen::Matrix2d &);
	void Update_Values();

	double getQ2();
	Eigen::Matrix2d getdQ2dF();
	Eigen::Tensor<double,4> getddQ2ddF();

	void set_Params(const double,const double,const Eigen::Matrix2d &);


private:
	double Emod;
	double nu;
	Eigen::Matrix2d F;

	Eigen::Matrix2d Identity2D;
	Eigen::Tensor<double,4> Cofactor4d;
	Eigen::Tensor<double,4> Identity4d;
	Eigen::Tensor<double,4> Identity4d2;


	double Q2;
	Eigen::Matrix2d dQ2dF;
	Eigen::Tensor<double,4> ddQ2ddF;


	void Calc_Q2();
	void Calc_dQ2dF();
	void Calc_ddQ2ddF();

	double Det2D(Eigen::Matrix2d &);
	Eigen::Matrix2d Cofactor2D(Eigen::Matrix2d &);
	Eigen::Matrix2d outer_product(const Eigen::Vector2d &, const Eigen::Vector2d &);
	double Tr(Eigen::Matrix2d &);
	double Tr2(Eigen::Matrix2d &);

};



#endif /* MATERIAL_CLASS_H_ */

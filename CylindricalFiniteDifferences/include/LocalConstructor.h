/*
 * LocalConstructor.h
 *
 *  Created on: Jan 25, 2022
 *      Author: kevin
 */

#ifndef CYLINDRICALFINITEDIFFERENCES_SRC_LOCALCONSTRUCTOR_H_
#define CYLINDRICALFINITEDIFFERENCES_SRC_LOCALCONSTRUCTOR_H_
#include "Element.h"
#include "material_class.h"

class Local_Constructor {
public:
	Local_Constructor();
	Local_Constructor(Element *, Material_Class *, Material_Class *);
	Local_Constructor(double, Element *, Material_Class *, Material_Class *);
	virtual ~Local_Constructor();
	void calculate_gradient();
	void calculate_hessian();
	Eigen::VectorXd getGradient();
	Eigen::MatrixXd getHessian();
	void setEa(Eigen::Matrix2d &);
	void setBa(Eigen::Matrix2d &);
	void seth(double);
	std::vector<int> get_indices();

private:
	Element * elem = NULL;
	Material_Class * Q1_mat = NULL;
	Material_Class * Q2_mat = NULL;

	std::vector<int> indices;

	Eigen::Matrix2d Ea;
	Eigen::Matrix2d Ba;

	int ndofs = 0;

	double Emod = 1.;
	double nu = 0.4;

	double h = 0.01;
	double hsc = 12. * pow(h,2)/(1.- pow(nu,2) );
	Eigen::VectorXd DEloc;
	Eigen::MatrixXd DDEloc;


	double tensor_inner(Eigen::Matrix2d &, Eigen::Matrix2d &);
	double bilinear_inner(Eigen::Matrix2d &, Eigen::Tensor<double,4> &,Eigen::Matrix2d &);
};

#endif /* CYLINDRICALFINITEDIFFERENCES_SRC_LOCALCONSTRUCTOR_H_ */

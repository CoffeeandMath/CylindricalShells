/*
 * LocalConstructor.cpp
 *
 *  Created on: Jan 25, 2022
 *      Author: kevin
 */

#include "LocalConstructor.h"

Local_Constructor::Local_Constructor() {
	// TODO Auto-generated constructor stub

}
Local_Constructor::Local_Constructor(Element * elem_, Material_Class * Q1_mat_,Material_Class * Q2_mat_){
	elem = elem_;
	indices = elem->get_indices();
	Q1_mat = Q1_mat_;
	Q2_mat = Q2_mat_;

	ndofs = elem->get_indices().size();
	Eigen::VectorXd etemp(ndofs);
	DEloc = etemp;
	Eigen::MatrixXd ddetemp(ndofs,ndofs);
	DDEloc = ddetemp;
	Ea.setZero();
	Ba.setZero();
}

Local_Constructor::Local_Constructor(double h_, Element * elem_, Material_Class * Q1_mat_,Material_Class * Q2_mat_){
	elem = elem_;
	indices = elem->get_indices();
	Q1_mat = Q1_mat_;
	Q2_mat = Q2_mat_;
	ndofs = elem->get_indices().size();
	Eigen::VectorXd etemp(ndofs);
	DEloc = etemp;
	Eigen::MatrixXd ddetemp(ndofs,ndofs);
	DDEloc = ddetemp;
	Ea.setZero();
	Ba.setZero();
	seth(h_);
}

void Local_Constructor::calculate_gradient(){

	elem->calc();

	Q1_mat->set_Params(Emod, nu, elem->getE() - Ea);
	Q2_mat->set_Params(Emod, nu, elem->getB() - Ba);
	Eigen::Matrix2d dQ1dF = Q1_mat->getdQ2dF();
	Eigen::Matrix2d dQ2dF = Q2_mat->getdQ2dF();
	double J = elem->getl();
	DEloc.setZero();
	for (unsigned int i = 0; i < ndofs; i++){
		elem->calc_i_perturbation(i);



		Eigen::Matrix2d DF = elem->getDiE();


		DEloc[i] += tensor_inner(dQ1dF,DF) * J;

		DF = elem->getDiB();

		DEloc[i] += hsc * tensor_inner(dQ2dF,DF) * J;

	}
}

void Local_Constructor::calculate_hessian(){
	std::cout << "Got Here 1" << std::endl;
	elem->calc();

	Q1_mat->set_Params(Emod, nu, elem->getE());
	Q2_mat->set_Params(Emod, nu, elem->getB());
	Eigen::Matrix2d dQ1dF = Q1_mat->getdQ2dF();
	Eigen::Matrix2d dQ2dF = Q2_mat->getdQ2dF();
	Eigen::Tensor<double,4> ddQ1ddF = Q1_mat->getddQ2ddF();
	Eigen::Tensor<double,4> ddQ2ddF = Q2_mat->getddQ2ddF();
	double J = elem->getl();
	DEloc.setZero();
	DDEloc.setZero();
	for (unsigned int i = 0; i < ndofs; i++){
		elem->calc_i_perturbation(i);

		Eigen::Matrix2d DiE = elem->getDiE();
		DEloc[i] += tensor_inner(dQ1dF,DiE)*J;


		Eigen::Matrix2d DiB = elem->getDiB();


		DEloc[i] += hsc * tensor_inner(dQ2dF,DiB) * J;

		for (unsigned int j = 0; j < ndofs; j++){
			elem->calc_j_perturbation(j);
			Eigen::Matrix2d DjE = elem->getDjE();
			Eigen::Matrix2d DjB = elem->getDjB();

			Eigen::Matrix2d DDE = elem->getDDE();
			Eigen::Matrix2d DDB = elem->getDDB();



			DDEloc(i,j) += tensor_inner(dQ1dF,DDE) * J;
			DDEloc(i,j) += bilinear_inner(DiE,ddQ1ddF,DjE) * J;

			DDEloc(i,j) += hsc* tensor_inner(dQ2dF,DDB) * J;
			DDEloc(i,j) += hsc*bilinear_inner(DiB,ddQ2ddF,DjB)*J;

		}

	}



}

double Local_Constructor::tensor_inner(Eigen::Matrix2d & m1, Eigen::Matrix2d & m2){
	double sum = 0.;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			sum += m1(i,j)*m2(i,j);
		}
	}
	return sum;
}

double Local_Constructor::bilinear_inner(Eigen::Matrix2d & m1 , Eigen::Tensor<double,4> & C,Eigen::Matrix2d & m2){
	double sum = 0.;

	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					sum += m1(i,j)*C(i,j,k,l)*m2(k,l);
				}
			}
		}
	}

	return sum;
}

void Local_Constructor::setEa(Eigen::Matrix2d & Ea_){
	Ea = Ea_;
}
void Local_Constructor::setBa(Eigen::Matrix2d & Ba_){
	Ba = Ba_;
}

void Local_Constructor::seth(double h_){
	h = h_;
	hsc = 12. * pow(h,2)/(1.- pow(nu,2) );
}



Eigen::VectorXd Local_Constructor::getGradient(){

	return DEloc;

}
Eigen::MatrixXd Local_Constructor::getHessian(){

	return DDEloc;

}

std::vector<int> Local_Constructor::get_indices(){
	return indices;
}

Local_Constructor::~Local_Constructor() {
	// TODO Auto-generated destructor stub
}


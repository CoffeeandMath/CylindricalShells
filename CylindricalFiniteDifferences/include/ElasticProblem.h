/*
 * ElasticProblem.h
 *
 *  Created on: Jan 26, 2022
 *      Author: kevin
 */

#ifndef CYLINDRICALFINITEDIFFERENCES_SRC_ELASTICPROBLEM_H_
#define CYLINDRICALFINITEDIFFERENCES_SRC_ELASTICPROBLEM_H_
#include "LocalConstructor.h"
#include "DOFHandler.h"
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Core>
#include <math.h>
#include <chrono>
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SpMat;

class ElasticProblem {
public:
	ElasticProblem();
	virtual ~ElasticProblem();
	void run();
private:
	void set_geometry();
	void newton_raphson();
	void gradient_descent();
	void solve_path();
	void solve();
	void fill_vector(Eigen::VectorXd *,const std::vector<int> &,const Eigen::VectorXd &);
	void fill_matrix(std::vector<Eigen::Triplet<double>> *,const std::vector<int> &,const Eigen::MatrixXd &);
	void calc_E();
	void calc_dE();
	void calc_ddE();
	Eigen::VectorXd numerical_derivative();

	Eigen::ConjugateGradient<SpMat> solver;
	double tol = 1.e-10;


	std::vector<Node> Nodes;
	std::vector<Reference_Configuration> References;

	std::vector<Element> Elements;
	std::vector<Material_Class> Q1_Material;
	std::vector<Material_Class> Q2_Material;
	std::vector<Local_Constructor> LocalConstructors;
	std::vector<Eigen::Matrix2d> Ea;
	std::vector<Eigen::Matrix2d> Ba;


	DOFHandler DOFHan;

	int DofPerNode = 2;
	int ndofs = 0;

	Eigen::VectorXd solution;
	double E = 0.;
	Eigen::VectorXd dE;
	Eigen::MatrixXd ddE;
	std::vector<T> ddEtrips;

	Eigen::VectorXd linearsolve;

	double R0 = 1.0;
	double L = 1.0;
	int Nnodes = 100;
	double Emod = 1.e3;
	double nu = 0.4;
	double h = 0.01;
	double hsc = 12.*pow(h,2)/(1.0 - pow(nu,2));

	double tensor_inner(Eigen::Matrix2d &, Eigen::Matrix2d &);
	double bilinear_inner(Eigen::Matrix2d &, Eigen::Tensor<double,4> &,Eigen::Matrix2d &);

};

#endif /* CYLINDRICALFINITEDIFFERENCES_SRC_ELASTICPROBLEM_H_ */

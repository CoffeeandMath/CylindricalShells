/*
 * ElasticProblem.cpp
 *
 *  Created on: Jan 26, 2022
 *      Author: kevin
 */

#include "ElasticProblem.h"

ElasticProblem::ElasticProblem() {
	// TODO Auto-generated constructor stub

}

void ElasticProblem::solve_path(){
	newton_raphson();
}


void ElasticProblem::set_geometry(){
	R0 = 1.0;
	L = 1.0;
	Nnodes = 100;

	Eigen::VectorXd Si = Eigen::VectorXd::LinSpaced(Nnodes,0.0,L);


	//Setting nodal level parameters
	Nodes.clear();
	References.clear();
	Reference_Configuration tempref;

	Q1_Material.clear();
	Q2_Material.clear();
	for (int i = 0; i < Nnodes; i++){
		std::vector<int> indices = {2*i,2*i+1};
		Eigen::Vector2d values = {R0+0.1,Si[i]};

		Nodes.emplace_back(Si[i], indices, values);


		tempref.set_R0(R0);
		tempref.set_point(Si[i]);
		References.push_back(tempref);

	}
	Nodes[0].setfixed(0);
	Nodes[0].setfixed(1);
	Nodes[1].setfixed(0);
	Nodes[1].setfixed(1);
	Nodes[2].setfixed(0);
	Nodes[2].setfixed(1);
	Nodes[3].setfixed(0);
	Nodes[3].setfixed(1);

	//Setting element level parameters
	int Nelements = Nnodes - 4;
	Elements.clear();

	for (int i = 0; i < Nelements; i++){




		Elements.emplace_back(&(Nodes[i]),&(Nodes[i+1]),&(Nodes[i+2]),&(Nodes[i+3]),&(Nodes[i+4]), &(References[i+2]));
		Elements[i].calc();

		Q1_Material.emplace_back(Emod, nu, Eigen::Matrix2d::Identity());
		Q2_Material.emplace_back(Emod, nu, Eigen::Matrix2d::Identity());

	}


	for (int i = 0; i < Elements.size(); i++){

		Elements[i].calc();

		for (int ii = 0; ii < 10; ii++){
			Elements[i].calc_i_perturbation(ii);
			for (int ij = 0; ij < 10; ij++){
				Elements[i].calc_j_perturbation(ij);
			}
		}

	}

	Ea.resize(Elements.size());
	Ba.resize(Elements.size());

	for (int i = 0; i < Ea.size(); i++){
		Ea[i].setZero();
		Ba[i].setZero();
	}

	ndofs = Nnodes*DofPerNode;
	Eigen::VectorXd dEtemp(ndofs); dEtemp.setZero();
	dE = dEtemp;
	solution = dEtemp;

	linearsolve = dEtemp;

	Eigen::MatrixXd ddEtemp(ndofs,ndofs); ddEtemp.setZero();
	ddE = ddEtemp;

	DOFHan = DOFHandler(&Nodes,&Elements);
	solution = DOFHan.get_solution();
}





void ElasticProblem::fill_vector(Eigen::VectorXd * vect,const std::vector<int> & ind, const Eigen::VectorXd & vals){
	for (int i = 0; i < ind.size(); i++){
		(*vect)[ind[i]] += vals[i];
	}
}

void ElasticProblem::fill_matrix(std::vector<Eigen::Triplet<double>> * tripmatr,const std::vector<int> & ind, const Eigen::MatrixXd & vals){


	for (int i = 0; i < ind.size(); i++){
		for (int j = 0; j < ind.size(); j++){
			(*tripmatr).push_back(T(ind[i],ind[j],vals(i,j)));
		}
	}


}


void ElasticProblem::calc_E(){
	E = 0.;


	int cellindex = 0;
	for (auto & elem : Elements){

		elem.calc();

		Q1_Material[cellindex].set_Params(Emod, nu, elem.getE());
		Q2_Material[cellindex].set_Params(Emod, nu, elem.getB());
		//std::cout << elem.getB() << std::endl;

		Eigen::Matrix2d dQ1dF = Q1_Material[cellindex].getdQ2dF();
		Eigen::Matrix2d dQ2dF = Q2_Material[cellindex].getdQ2dF();



		//std::cout << E << std::endl;

		int nsize = elem.get_indices().size();
		double J = elem.getl();


		E += Q1_Material[cellindex].getQ2() * J;
		E += hsc *Q2_Material[cellindex].getQ2() * J;

		cellindex++;
	}
}


void ElasticProblem::calc_dE(){
	E = 0.;
	dE.setZero();

	Eigen::VectorXd DEloc(Elements[0].get_indices().size());
	int cellindex = 0;
	for (auto & elem : Elements){

		DEloc.setZero();
		elem.calc();

		Q1_Material[cellindex].set_Params(Emod, nu, elem.getE());
		Q2_Material[cellindex].set_Params(Emod, nu, elem.getB());
		//std::cout << elem.getB() << std::endl;

		Eigen::Matrix2d dQ1dF = Q1_Material[cellindex].getdQ2dF();
		Eigen::Matrix2d dQ2dF = Q2_Material[cellindex].getdQ2dF();



		//std::cout << E << std::endl;

		int nsize = elem.get_indices().size();
		double J = elem.getl();


		E += Q1_Material[cellindex].getQ2() * J;
		E += hsc *Q2_Material[cellindex].getQ2() * J;
		for (int i = 0; i < nsize; i++){

			elem.calc_i_perturbation(i);


			Eigen::Matrix2d DiE = elem.getDiE();
			Eigen::Matrix2d DiB = elem.getDiB();


			DEloc[i] += tensor_inner(dQ1dF,DiE)*J;

			DEloc[i] +=  hsc * tensor_inner(dQ2dF,DiB) * J;




		}

		fill_vector(&dE ,elem.get_indices(), DEloc);
		cellindex++;
	}


}

void ElasticProblem::calc_ddE(){
	E = 0.;
	dE.setZero();
	ddE.setZero();
	ddEtrips.clear(); ddEtrips.reserve(100*Elements.size());
	Eigen::VectorXd DEloc(Elements[0].get_indices().size());
	Eigen::MatrixXd DDEloc(Elements[0].get_indices().size(),Elements[0].get_indices().size());
	int cellindex = 0;
	for (auto & elem : Elements){

		DEloc.setZero();
		DDEloc.setZero();
		elem.calc();

		Q1_Material[cellindex].set_Params(Emod, nu, elem.getE());
		Q2_Material[cellindex].set_Params(Emod, nu, elem.getB());
		//std::cout << elem.getB() << std::endl;

		Eigen::Matrix2d dQ1dF = Q1_Material[cellindex].getdQ2dF();
		Eigen::Matrix2d dQ2dF = Q2_Material[cellindex].getdQ2dF();
		Eigen::Tensor<double,4> ddQ1ddF = Q1_Material[cellindex].getddQ2ddF();
		Eigen::Tensor<double,4> ddQ2ddF = Q2_Material[cellindex].getddQ2ddF();



		//std::cout << E << std::endl;

		int nsize = elem.get_indices().size();
		double J = elem.getl();


		E += Q1_Material[cellindex].getQ2() * J;
		E += hsc *Q2_Material[cellindex].getQ2() * J;
		for (int i = 0; i < nsize; i++){

			elem.calc_i_perturbation(i);


			Eigen::Matrix2d DiE = elem.getDiE();
			Eigen::Matrix2d DiB = elem.getDiB();


			DEloc[i] += tensor_inner(dQ1dF,DiE)*J;

			DEloc[i] +=  hsc * tensor_inner(dQ2dF,DiB) * J;

			for (int j = 0; j < nsize; j++){
				elem.calc_j_perturbation(j);

				Eigen::Matrix2d DjE = elem.getDjE();
				Eigen::Matrix2d DjB = elem.getDjB();

				Eigen::Matrix2d DDE = elem.getDDE();
				Eigen::Matrix2d DDB = elem.getDDB();



				DDEloc(i,j) += tensor_inner(dQ1dF,DDE) * J;
				DDEloc(i,j) += bilinear_inner(DiE,ddQ1ddF,DjE) * J;

				DDEloc(i,j) += hsc* tensor_inner(dQ2dF,DDB) * J;
				DDEloc(i,j) += hsc*bilinear_inner(DiB,ddQ2ddF,DjB)*J;


			}


		}

		fill_vector(&dE ,elem.get_indices(), DEloc);
		fill_matrix(&ddEtrips,elem.get_indices(), DDEloc);
		cellindex++;
	}

}



double ElasticProblem::tensor_inner(Eigen::Matrix2d & m1, Eigen::Matrix2d & m2){
	double sum = 0.;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			sum += m1(i,j)*m2(i,j);
		}
	}
	return sum;
}

double ElasticProblem::bilinear_inner(Eigen::Matrix2d & m1 , Eigen::Tensor<double,4> & C,Eigen::Matrix2d & m2){
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

void ElasticProblem::newton_raphson(){
	double stepsize = 2.0*tol;
	int cntr = 0;


	solution = DOFHan.get_solution();
	while (stepsize > tol){
		cntr++;
		calc_ddE();
		solve();
		stepsize = linearsolve.norm()/linearsolve.size();
		std::cout << stepsize << std::endl;
	}
}

void ElasticProblem::solve(){
	linearsolve.setZero();

	SpMat ddEsp(DOFHan.n_free_dofs(),DOFHan.n_free_dofs());
	std::vector<T> ddEtripsconst = DOFHan.condense(ddEtrips);
	ddEsp.setFromTriplets(ddEtripsconst.begin(), ddEtripsconst.end());
	ddEsp.makeCompressed();

	solver.analyzePattern(ddEsp);
	solver.factorize(ddEsp);
	if (solver.info() != Eigen::Success) {
		std::cout << "decomposition failed" << std::endl;
		return;
	}
	Eigen::VectorXd dEcond = DOFHan.condense(dE);
	Eigen::VectorXd linearsolveconst(DOFHan.n_free_dofs());

	linearsolveconst = -1.*solver.solve(dEcond);
	linearsolve = DOFHan.distribute(linearsolveconst);
	//linearsolve = -1e-4*DOFHan.apply(dE);

	//std::cout << linearsolve << std::endl;
	solution += linearsolve;
	DOFHan.apply_solution(solution);

}


void ElasticProblem::gradient_descent(){
	double stepsize = 2.0*tol;
	int cntr = 0;

	double maxpivot = -1.e-6;

	solution = DOFHan.get_solution();
	while (stepsize > tol){
		cntr++;
		std::cout << cntr << std::endl;
		//linearsolve = numerical_derivative();
		calc_dE();
		linearsolve = dE;
		linearsolve = DOFHan.apply(linearsolve);

		stepsize = fabs(maxpivot)*linearsolve.norm()/linearsolve.size();

		solution += maxpivot*linearsolve;
		std::cout << "stepsize = " << stepsize << ", tol = " << tol << std::endl;
		DOFHan.apply_solution(solution);
	}
}


Eigen::VectorXd ElasticProblem::numerical_derivative(){
	Eigen::VectorXd ndE(ndofs); ndE.setZero();
	Eigen::VectorXd pert(ndofs);
	Eigen::VectorXd origsol = solution;
	double dh = 1e-10;
	for (int i = 0; i < ndofs; i++){
		pert.setZero();
		pert[i] = dh/2.;

		Eigen::VectorXd xp = origsol + pert;
		Eigen::VectorXd xm = origsol - pert;

		DOFHan.apply_solution(xp);
		calc_E();
		double Ep = E;

		DOFHan.apply_solution(xm);
		calc_E();
		double Em = E;

		ndE[i] = (Ep - Em)/dh;

	}
	solution = origsol;
	DOFHan.apply_solution(solution);
	return ndE;

}


void ElasticProblem::run(){
	set_geometry();
	/*
	calc_ddE();
	Eigen::VectorXd dEan = dE;
	Eigen::VectorXd dEnum = numerical_derivative();

	std::cout << (dEan - dEnum).norm() << std::endl;
	 */

	solve_path();
	std::cout << "Rsol" << std::endl;
	for (int i = 0; i < Nodes.size(); i++){
		std::cout << Nodes[i].get(0) << std::endl;
	}


}


ElasticProblem::~ElasticProblem() {
	// TODO Auto-generated destructor stub
}



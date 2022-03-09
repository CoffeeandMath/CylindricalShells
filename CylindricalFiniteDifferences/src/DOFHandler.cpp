/*
 * DOFHandler.cpp
 *
 *  Created on: Jan 31, 2022
 *      Author: kevin
 */

#include "DOFHandler.h"

DOFHandler::DOFHandler() {
	// TODO Auto-generated constructor stub
	NodeP = NULL;
	ElemP = NULL;
}

DOFHandler::DOFHandler(std::vector<Node> * np, std::vector<Element> * ep) {
	// TODO Auto-generated constructor stub
	NodeP = np;
	ElemP = ep;
	ndofs = dofspernode*NodeP->size();

	Eigen::VectorXd soltemp(ndofs); soltemp.setZero();
	solution = soltemp;


	globaltonodenumber.resize(ndofs);
	globaltonodeindex.resize(ndofs);
	isdofconstrained.resize(ndofs);
	globaltoconstrained.resize(ndofs);
	nfree = ndofs;
	int freecntr = 0;
	for (int i = 0; i < NodeP->size(); i++){
		for (int j = 0; j < (*NodeP)[i].getindex().size(); j++){
			globaltonodenumber[(*NodeP)[i].getindex()[j]] = i;
			globaltonodeindex[(*NodeP)[i].getindex()[j]] = j;

			if ((*NodeP)[i].isfree(j)){
				isdofconstrained[i] = false;
				globaltoconstrained[(*NodeP)[i].getindex()[j]] = freecntr;
				freecntr++;
			} else {
				isdofconstrained[i] = true;
				globaltoconstrained[(*NodeP)[i].getindex()[j]] = -1;
				nfree--;
			}
		}
	}

	isdofconstrained.resize(ndofs);
}

Eigen::VectorXd DOFHandler::condense(Eigen::VectorXd & sol){
	Eigen::VectorXd solout(nfree); solout.setZero();
	int cntr = 0;
	for (int i = 0; i < sol.size(); i++){
		if (globaltoconstrained[i] > -1){
			solout[globaltoconstrained[i]] = sol[i];
		}
	}


	return solout;
}

std::vector<T> DOFHandler::condense(std::vector<T> & tripsin){
	std::vector<T> tripsout; tripsout.reserve(tripsin.size());


	for (int i = 0; i < tripsin.size(); i++){
		if (globaltoconstrained[tripsin[i].row()] > -1 && globaltoconstrained[tripsin[i].col()] > -1){
			tripsout.push_back(T(globaltoconstrained[tripsin[i].row()],globaltoconstrained[tripsin[i].col()], tripsin[i].value()));
		}
	}

	return tripsout;
}

Eigen::VectorXd DOFHandler::distribute(Eigen::VectorXd & sol){
	Eigen::VectorXd solout(ndofs); solout.setZero();

	for (int i = 0; i < ndofs; i++){
		if (globaltoconstrained[i] > -1){
			solout[i] = sol[globaltoconstrained[i]];
		}

	}


	return solout;

}

Eigen::VectorXd DOFHandler::apply(Eigen::VectorXd & sol){
	Eigen::VectorXd solout(ndofs); solout.setZero();
	for (int i = 0; i < ndofs; i++){
			if (globaltoconstrained[i] > -1){
				solout[i] = sol[i];
			}

		}

	return solout;
}


std::vector<T> DOFHandler::apply(std::vector<T> & tripsin){
	std::vector<T> tripsout(tripsin.size());

	for (int i = 0; i < tripsin.size(); i++){
		if (globaltoconstrained[tripsin[i].row()] > -1 && globaltoconstrained[tripsin[i].col()] > -1){
			tripsout[i] = tripsin[i];
		} else {
			tripsout[i] = T(tripsin[i].row(),tripsin[i].col(),0.);
		}
	}
	return tripsout;


}


int DOFHandler::n_free_dofs(){
	return nfree;
}


void DOFHandler::apply_solution(const Eigen::VectorXd & sol){
	solution = sol;

	for (int i = 0; i < solution.size(); i++){
		((*NodeP)[globaltonodenumber[i]]).set(globaltonodeindex[i],solution[i]);
	}
}

Eigen::VectorXd DOFHandler::get_solution(){
	Eigen::VectorXd solout(ndofs); solout.setZero();

	for (int i = 0; i < NodeP->size(); i++){
		for (int j = 0; j < ((*NodeP)[i]).getindex().size(); j++){
			solout[((*NodeP)[i].getindex())[j]] += ((*NodeP)[i]).get(j);
		}
	}


	return solout;
}
DOFHandler::~DOFHandler() {
	// TODO Auto-generated destructor stub
}


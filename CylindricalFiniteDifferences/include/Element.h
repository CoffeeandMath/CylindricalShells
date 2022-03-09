/*
 * Element.h
 *
 *  Created on: Jan 24, 2022
 *      Author: kevin
 */

#ifndef CYLINDRICALFINITEDIFFERENCES_SRC_ELEMENT_H_
#define CYLINDRICALFINITEDIFFERENCES_SRC_ELEMENT_H_
#include <vector>
#include "Node.h"
#include "reference_configuration.h"
#include "unsupported/Eigen/CXX11/Tensor"


using Eigen::Tensor;
using Eigen::TensorMap;
using Eigen::TensorRef;


class Element {
public:
	Element();
	Element(std::vector<Node*> & , Reference_Configuration *);
	Element(Node*,Node*,Node*,Node*,Node*, Reference_Configuration *);
	virtual ~Element();
	void calc();
	void update_weights();
	void calc_i_perturbation(int);
	void calc_j_perturbation(int);
	std::vector<int> getindex();
	Eigen::Vector2d getu();
	Eigen::Vector2d getdu();
	Eigen::Vector2d getddu();
	Eigen::Matrix2d getac();
	Eigen::Matrix2d getaC();
	Eigen::Matrix2d getE();
	Eigen::Matrix2d getB();
	Eigen::Matrix2d getDiE();
	Eigen::Matrix2d getDjE();
	Eigen::Matrix2d getDiB();
	Eigen::Matrix2d getDjB();
	Eigen::Matrix2d getDDE();
	Eigen::Matrix2d getDDB();
	double getl();
	std::vector<int> get_indices();

private:

	std::vector<Node*> NodeP;

	Reference_Configuration *Ref = NULL;

	double l = 1.;
	Node * n0p = NULL;
	Node * n1p = NULL;
	Node * n2p = NULL;
	Node * n3p = NULL;
	Node * n4p = NULL;

	Eigen::Vector2d ui;
	Eigen::Vector2d uip;
	Eigen::Vector2d uipp;

	Eigen::Matrix2d ac;
	Eigen::Matrix2d aC;

	Eigen::Vector3d f;

	Eigen::Matrix2d bc;

	Eigen::VectorXd u_weights;
	Eigen::VectorXd du_weights;
	Eigen::VectorXd ddu_weights;

	std::vector<int> index;



	Eigen::Vector3d a1;
	Eigen::Vector3d a2;

	Eigen::Vector3d n;

	Eigen::Matrix2d E;
	Eigen::Matrix2d B;

	void set_i_perturbation(int);


	Eigen::Matrix2d symmetric(Eigen::Matrix2d &);

	std::vector<Eigen::Vector2d> Uivec;

	Eigen::Vector2d Diui;
	Eigen::Vector2d Diuip;
	Eigen::Vector2d Diuipp;

	Eigen::Vector3d DiddxddS;
	Eigen::Vector3d DiddxdSdtheta;
	Eigen::Vector3d Diddxddtheta;

	Eigen::Vector3d Dif;

	Eigen::Vector3d Dia1;
	Eigen::Vector3d Dia2;

	Eigen::Matrix2d Diac;
	Eigen::Matrix2d DiaC;

	Eigen::Vector3d Din;

	Eigen::Matrix2d Dibc;

	Eigen::Matrix2d DiE;
	Eigen::Matrix2d DiB;

	void set_j_perturbation(int);


	std::vector<Eigen::Vector2d> Ujvec;

	Eigen::Vector2d Djui;
	Eigen::Vector2d Djuip;
	Eigen::Vector2d Djuipp;

	Eigen::Vector3d DjddxddS;
	Eigen::Vector3d DjddxdSdtheta;
	Eigen::Vector3d Djddxddtheta;

	Eigen::Vector3d Djf;

	Eigen::Vector3d Dja1;
	Eigen::Vector3d Dja2;

	Eigen::Matrix2d Djac;
	Eigen::Matrix2d DjaC;
	Eigen::Vector3d Djn;

	Eigen::Matrix2d Djbc;

	Eigen::Matrix2d DjE;
	Eigen::Matrix2d DjB;

	Eigen::Matrix2d DDac;
	Eigen::Matrix2d DDaC;

	Eigen::Vector3d DDf;
	Eigen::Vector3d DDn;

	Eigen::Matrix2d DDbc;

	Eigen::Matrix2d DDE;
	Eigen::Matrix2d DDB;



};

#endif /* CYLINDRICALFINITEDIFFERENCES_SRC_ELEMENT_H_ */

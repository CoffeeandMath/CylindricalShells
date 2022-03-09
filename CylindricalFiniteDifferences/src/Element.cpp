/*
 * Element.cpp
 *
 *  Created on: Jan 24, 2022
 *      Author: kevin
 */

#include "Element.h"

Element::Element() {
	// TODO Auto-generated constructor stub
}

Element::Element(std::vector<Node*> & NodeP_,Reference_Configuration * Ref_){
	NodeP.clear();
	NodeP.resize(NodeP_.size());
	for (int i = 0; i < NodeP_.size(); i++){
		NodeP[i] = NodeP_[i];
	}



	Ref = Ref_;
	index.clear();
	for (int i = 0; i < NodeP.size(); i++){
		for (int j = 0; j < (NodeP[i]->getindex()).size(); j++){
			index.push_back((NodeP[i]->getindex())[j]);

		}
	}



	l = (NodeP[1]->getS() - NodeP[0]->getS());

	Eigen::VectorXd etemp(5); etemp.setZero();
	u_weights = etemp;
	du_weights = etemp;
	ddu_weights = etemp;

	update_weights();
	calc();
}

Element::Element(Node * n0, Node * n1,Node * n2,Node * n3,Node * n4, Reference_Configuration * Ref_){
	n0p = n0;
	n1p = n1;
	n2p = n2;
	n3p = n3;
	n4p = n4;

	Ref = Ref_;

	index.clear();
	index.resize(10);
	index[0] = n0p->getindex()[0];index[1] = n0p->getindex()[1];
	index[2] = n1p->getindex()[0];index[3] = n1p->getindex()[1];
	index[4] = n2p->getindex()[0];index[5] = n2p->getindex()[1];
	index[6] = n3p->getindex()[0];index[7] = n3p->getindex()[1];
	index[8] = n4p->getindex()[0];index[9] = n4p->getindex()[1];


	l = (n1p->getS() - n0p->getS());

	Eigen::VectorXd etemp(5); etemp.setZero();
	u_weights = etemp;
	du_weights = etemp;
	ddu_weights = etemp;

	update_weights();
	calc();
}


void Element::calc(){
	Eigen::Vector2d zerovec; zerovec.setZero();
	ui = zerovec;
	uip = zerovec;
	uipp = zerovec;

	//std::cout << NodeP.size() << std::endl;
	update_weights();

	/*
	for (int i = 0; i < NodeP.size(); i++){
		std::cout << "Got Here 2" << std::endl;
		ui += u_weights[i]*(NodeP[i]->get());
		uip += du_weights[i]*(NodeP[i]->get());
		uipp += ddu_weights[i]*(NodeP[i]->get());
	}
	 */
	ui = u_weights[0] * (n0p->get())
			+ u_weights[1] * (n1p->get())
			+ u_weights[2] * (n2p->get())
			+ u_weights[3] * (n3p->get())
			+ u_weights[4] * (n4p->get());

	uip = du_weights[0] * (n0p->get())
			+ du_weights[1] * (n1p->get())
			+ du_weights[2] * (n2p->get())
			+ du_weights[3] * (n3p->get())
			+ du_weights[4] * (n4p->get());

	uipp = ddu_weights[0] * (n0p->get())
			+ ddu_weights[1] * (n1p->get())
			+ ddu_weights[2] * (n2p->get())
			+ ddu_weights[3] * (n3p->get())
			+ ddu_weights[4] * (n4p->get());


	a1[0] = uip[0];
	a1[1] = 0.;
	a1[2] = uip[1];

	a2[0] = 0.;
	a2[1] = ui[0];
	a2[2] = 0.;

	ac(0,0) = a1.dot(a1);
	ac(0,1) = a1.dot(a2);
	ac(1,0) = ac(0,1);
	ac(1,1) = a2.dot(a2);

	aC = ac.inverse();

	f = a1.cross(a2);
	n = f/f.norm();

	Eigen::Vector3d ddxddS = {uipp[0], 0., uipp[1]};
	Eigen::Vector3d ddxdSdtheta = {0., uip[0], 0.};
	Eigen::Vector3d ddxddtheta = {-1.*ui[0], 0., 0.};
	bc(0,0) = n.dot(ddxddS);
	bc(0,1) = n.dot(ddxdSdtheta);
	bc(1,0) = bc(0,1);
	bc(1,1) = n.dot(ddxddtheta);

	E = 0.5*(Ref->get_P())*(Ref->get_AC() * ac * Ref->get_AC() - Ref->get_AC())*(Ref->get_P().transpose());

	B = Ref->get_P()*aC*bc*aC*(Ref->get_P().transpose());

}

void Element::update_weights(){


	u_weights[0] = 0.;
	u_weights[1] = 0.;
	u_weights[2] = 1.;
	u_weights[3] = 0.;
	u_weights[4] = 0.;


	du_weights[0] = 0.;
	du_weights[1] = -1./(2.*l);
	du_weights[2] = 0.;
	du_weights[3] = 1./(2.*l);
	du_weights[4] = 0.;

	ddu_weights[0] = 1./pow(2.*l,2);
	ddu_weights[1] = 0.;
	ddu_weights[2] = -2./pow(2.*l,2);
	ddu_weights[3] = 0.;
	ddu_weights[4] = 1./pow(2.*l,2);
}

Eigen::Vector2d Element::getu(){
	return ui;
}

Eigen::Vector2d Element::getdu(){
	return uip;
}

Eigen::Vector2d Element::getddu(){
	return uipp;
}


Eigen::Matrix2d Element::getac(){
	return ac;
}

Eigen::Matrix2d Element::getaC(){
	return aC;
}


void Element::set_i_perturbation(int i){

	Uivec.resize(5);
	for (int k = 0; k < Uivec.size(); k++){
		Uivec[k].setZero();
	}

	Uivec[i/2][i%2] = 1.0;


}

void Element::calc_i_perturbation(int i){
	set_i_perturbation(i);


	Diui.setZero();
	Diuip.setZero();
	Diuipp.setZero();

	for (int k = 0; k < Uivec.size(); k++){
		Diui += u_weights[k]*Uivec[k];
		Diuip += du_weights[k]*Uivec[k];
		Diuipp += ddu_weights[k]*Uivec[k];
	}


	Dia1[0] = Diuip[0];
	Dia1[1] = 0.;
	Dia1[2] = Diuip[1];

	Dia2[0] = 0.;
	Dia2[1] = Diui[0];
	Dia2[2] = 0.;

	Diac(0,0) = 2.0*a1.dot(Dia1);
	Diac(0,1) = a1.dot(Dia2) + a2.dot(Dia1);
	Diac(1,0) = Diac(0,1);
	Diac(1,1) = 2.0*a2.dot(Dia2);

	DiaC = -1.0*aC*Diac*aC;

	Dif = Dia1.cross(a2) + a1.cross(Dia2);
	Din = (Eigen::Matrix3d::Identity() - n*n.transpose())*(Dif) / f.norm();

	Eigen::Vector3d ddxddS = {uipp[0], 0., uipp[1]};
	Eigen::Vector3d ddxdSdtheta = {0., uip[0], 0.};
	Eigen::Vector3d ddxddtheta = {-1.*ui[0], 0., 0.};

	DiddxddS[0] = Diuipp[0]; DiddxddS[1] =  0.; DiddxddS[2] = Diuipp[1];
	DiddxdSdtheta[0] = 0.; DiddxdSdtheta[1] = Diuip[0]; DiddxdSdtheta[2] =  0.;
	Diddxddtheta[0] = -1.*Diui[0]; Diddxddtheta[1] = 0.; Diddxddtheta[2] = 0.;

	Eigen::Matrix2d Dib1;
	Dib1(0,0) = Din.dot(ddxddS);
	Dib1(0,1) = Din.dot(ddxdSdtheta);
	Dib1(1,0) = Dib1(0,1);
	Dib1(1,1) = Din.dot(ddxddtheta);

	Eigen::Matrix2d Dib2;
	Dib2(0,0) = n.dot(DiddxddS);
	Dib2(0,1) = n.dot(DiddxdSdtheta);
	Dib2(1,0) = Dib2(0,1);
	Dib2(1,1) = n.dot(Diddxddtheta);

	Dibc = Dib1 + Dib2;

	DiE = 0.5 *  (Ref->get_P()) * (Ref->get_AC()) * Diac * (Ref->get_AC()) * (Ref->get_P().transpose());
	Eigen::Matrix2d DDtemp = aC*bc*DiaC;
	DiB = Ref->get_P() * (aC * Dibc * aC + 2. * symmetric(DDtemp)) * (Ref->get_P().transpose());


}

void Element::calc_j_perturbation(int j){
	set_j_perturbation(j);

	Djui.setZero();
	Djuip.setZero();
	Djuipp.setZero();

	for (int k = 0; k < Ujvec.size(); k++){
		Djui += u_weights[k]*Ujvec[k];
		Djuip += du_weights[k]*Ujvec[k];
		Djuipp += ddu_weights[k]*Ujvec[k];
	}




	Dja1[0] = Djuip[0];
	Dja1[1] = 0.;
	Dja1[2] = Djuip[1];

	Dja2[0] = 0.;
	Dja2[1] = Djui[0];
	Dja2[2] = 0.;

	Djac(0,0) = 2.0*a1.dot(Dja1);
	Djac(0,1) = a1.dot(Dja2) + a2.dot(Dja1);
	Djac(1,0) = Djac(0,1);
	Djac(1,1) = 2.0*a2.dot(Dja2);

	DjaC = -1.0*aC*Djac*aC;

	Djf = Dja1.cross(a2) + a1.cross(Dja2);
	Djn = (Eigen::Matrix3d::Identity() - n*n.transpose())*(Djf) / f.norm();

	Eigen::Vector3d ddxddS = {uipp[0], 0., uipp[1]};
	Eigen::Vector3d ddxdSdtheta = {0., uip[0], 0.};
	Eigen::Vector3d ddxddtheta = {-1.*ui[0], 0., 0.};

	DjddxddS[0] = Djuipp[0]; DjddxddS[1] = 0.; DjddxddS[2] = Djuipp[1];
	DjddxdSdtheta[0] = 0.; DjddxdSdtheta[1] = Djuip[0]; DjddxdSdtheta[2] = 0.;
	Djddxddtheta[0] = -1.*Djui[0]; Djddxddtheta[1] = 0.; Djddxddtheta[2] = 0.;

	Eigen::Matrix2d Djb1;
	Djb1(0,0) = Djn.dot(ddxddS);
	Djb1(0,1) = Djn.dot(ddxdSdtheta);
	Djb1(1,0) = Djb1(0,1);
	Djb1(1,1) = Djn.dot(ddxddtheta);

	Eigen::Matrix2d Djb2;
	Djb2(0,0) = n.dot(DjddxddS);
	Djb2(0,1) = n.dot(DjddxdSdtheta);
	Djb2(1,0) = Djb2(0,1);
	Djb2(1,1) = n.dot(Djddxddtheta);

	Djbc = Djb1 + Djb2;

	DjE = 0.5 *  Ref->get_P() * Ref->get_AC() * Djac * Ref->get_AC() * (Ref->get_P().transpose());
	Eigen::Matrix2d DDtemp = aC*bc*DjaC;
	DjB = Ref->get_P() * (aC * Djbc * aC + 2. * symmetric(DDtemp)) * (Ref->get_P().transpose());
	// Make sure i is already updated

	DDac(0,0) = 2.0*Dja1.dot(Dia1);
	DDac(0,1) = Dja1.dot(Dia2) + Dja2.dot(Dia1);
	DDac(1,0) = DDac(0,1);
	DDac(1,1) = 2.0*Dja2.dot(Dia2);


	DDtemp = aC*Diac*DjaC;
	DDaC = -aC*DDac*aC - 2.*symmetric(DDtemp);


	DDf = Dja1.cross(Dia2) + Dia1.cross(Dja2);

	DDn = -(n.dot(Djf))*Dif / f.norm()
							- (Djn*n.transpose() + n*Djn.transpose())*Dif / f.norm()
							+(Eigen::Matrix3d::Identity() - n*n.transpose())*DDf / f.norm();

	Eigen::Matrix2d DDb1;
	DDb1(0,0) = DDn.dot(ddxddS);
	DDb1(0,1) = DDn.dot(ddxdSdtheta);
	DDb1(1,0) = DDb1(0,1);
	DDb1(1,1) = DDn.dot(ddxddtheta);

	Eigen::Matrix2d DDb2;
	DDb2(0,0) = Din.dot(DjddxddS);
	DDb2(0,1) = Din.dot(DjddxdSdtheta);
	DDb2(1,0) = DDb2(0,1);
	DDb2(1,1) = Din.dot(Djddxddtheta);

	Eigen::Matrix2d DDb3;
	DDb3(0,0) = Djn.dot(DiddxddS);
	DDb3(0,1) = Djn.dot(DiddxdSdtheta);
	DDb3(1,0) = DDb3(0,1);
	DDb3(1,1) = Djn.dot(Diddxddtheta);

	DDbc = DDb1 + DDb2 + DDb3;


	DDE = 0.5*(Ref->get_P())*(Ref->get_AC() )* DDac * Ref->get_AC()*(Ref->get_P().transpose());

	DDtemp = DjaC*bc*DiaC + aC*(Dibc*DjaC + Djbc*DiaC + bc * DDaC);

	DDB = Ref->get_P() * (aC * DDbc * aC + 2.*symmetric(DDtemp)) * (Ref->get_P().transpose());

}


void Element::set_j_perturbation(int j){
	Ujvec.resize(5);
	for (int k = 0; k < Ujvec.size(); k++){
		Ujvec[k].setZero();
	}

	Ujvec[j/2][j%2] = 1.0;

}

Eigen::Matrix2d Element::symmetric(Eigen::Matrix2d & F){
	return 0.5*(F + F.transpose());
}


Eigen::Matrix2d Element::getE(){
	return E;
}
Eigen::Matrix2d Element::getB(){
	return B;
}
Eigen::Matrix2d Element::getDiE(){
	return DiE;
}
Eigen::Matrix2d Element::getDjE(){
	return DjE;
}
Eigen::Matrix2d Element::getDiB(){
	return DiB;
}
Eigen::Matrix2d Element::getDjB(){
	return DjB;
}
Eigen::Matrix2d Element::getDDE(){
	return DDE;
}
Eigen::Matrix2d Element::getDDB(){
	return DDB;
}

std::vector<int> Element::get_indices(){
	return index;
}

double Element::getl(){
	return l;
}
Element::~Element() {
	// TODO Auto-generated destructor stub
}


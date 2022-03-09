
#ifndef MATERIAL_CLASS_CC_
#define MATERIAL_CLASS_CC_
#include "material_class.h"







Material_Class::Material_Class() {
	Emod = 1.0;

	nu = 0.0;
	Cofactor4d = Eigen::Tensor<double, 4> (2,2,2,2); Cofactor4d.setZero();
	Identity4d = Eigen::Tensor<double, 4> (2,2,2,2); Identity4d.setZero();
	Identity4d2 = Eigen::Tensor<double, 4> (2,2,2,2); Identity4d2.setZero();
	ddQ2ddF = Eigen::Tensor<double, 4> (2,2,2,2); ddQ2ddF.setZero();

	//Setting values for the cofactor and identity 4d matrices

	for (unsigned int i = 0; i < 2; i++){
		Eigen::Vector2d ei;
		ei[i] = 1.0;
		Identity2D(i,i) = 1.0;
		for (unsigned int j = 0; j < 2; j++){
			Eigen::Vector2d ej;
			ej[j] = 1.0;


			for (unsigned int k = 0; k < 2; k++) {
				Eigen::Vector2d ek;
				ek[k] = 1.0;
				for (unsigned int l = 0; l < 2; l++){
					if (i==j && k==l){
						Identity4d(i,j,k,l) = 1.0;
					}

					if (i==k && j==l){
						Identity4d2(i,j,k,l) = 1.0;
					}
					Eigen::Vector2d el;
					el[l] = 1.0;

					Eigen::Matrix2d ekotimesel = outer_product(ek,el);
					Eigen::Matrix2d cofacekel = Cofactor2D(ekotimesel);
					Cofactor4d(i,j,k,l) = ei.dot(cofacekel*ej);


				}
			}
		}
	}
}

Material_Class::Material_Class(const double Emodtemp, const double nutemp, const Eigen::Matrix2d & Ftemp) {
	Emod = Emodtemp;
	nu = nutemp;
	F = Ftemp;

	Cofactor4d = Eigen::Tensor<double, 4> (2,2,2,2); Cofactor4d.setZero();
	Identity4d = Eigen::Tensor<double, 4> (2,2,2,2); Identity4d.setZero();
	Identity4d2 = Eigen::Tensor<double, 4> (2,2,2,2); Identity4d2.setZero();
	ddQ2ddF = Eigen::Tensor<double, 4> (2,2,2,2); ddQ2ddF.setZero();


	//std::cout << F.norm_square() << std::endl;

	for (unsigned int i = 0; i < 2; i++){
		Eigen::Vector2d ei; ei.setZero();
		ei[i] = 1.0;
		Identity2D(i,i) = 1.0;
		for (unsigned int j = 0; j < 2; j++){
			Eigen::Vector2d ej; ej.setZero();
			ej[j] = 1.0;


			for (unsigned int k = 0; k < 2; k++) {
				Eigen::Vector2d ek; ek.setZero();
				ek[k] = 1.0;
				for (unsigned int l = 0; l < 2; l++){
					if (i==j && k==l){
						Identity4d(i,j,k,l) = 1.0;
					}
					Eigen::Vector2d el; el.setZero();
					el[l] = 1.0;

					Eigen::Matrix2d ekotimesel = outer_product(ek,el);
					Eigen::Matrix2d cofacekel = Cofactor2D(ekotimesel);
					Cofactor4d(i,j,k,l) = ei.dot(cofacekel*ej);


				}
			}
		}
	}


	Calc_Q2();
	Calc_dQ2dF();
	Calc_ddQ2ddF();


}

void Material_Class::set_Params(const double Emodtemp, const double nutemp, const Eigen::Matrix2d & Ftemp) {
	Emod = Emodtemp;
	nu = nutemp;
	F = Ftemp;
	//std::cout << F.norm_square() << std::endl;


	Calc_Q2();
	Calc_dQ2dF();
	Calc_ddQ2ddF();


}

void Material_Class::Calc_Q2() {
	//Q2 = 0.5*Emod*Tr2(F);
	Q2 = 0.;
	Q2 = 0.5*Emod*(pow(Tr(F),2.0) - 2.0*(1.0 - nu)*Det2D(F));
}


void Material_Class::Calc_dQ2dF() {
	dQ2dF.setZero();
	//dQ2dF = Emod*F;
	dQ2dF = Emod*(Tr(F)*Identity2D - (1.0 - nu)*Cofactor2D(F));

}

void Material_Class::Calc_ddQ2ddF() {
	ddQ2ddF.setZero();
	//ddQ2ddF = Emod*(Identity4d2 );
	ddQ2ddF = Emod*(Identity4d - (1.0 - nu)*Cofactor4d);
}

double Material_Class::getQ2(){
	return Q2;
}

Eigen::Matrix2d Material_Class::getdQ2dF(){
	return dQ2dF;
}

Eigen::Tensor<double,4> Material_Class::getddQ2ddF(){
	return ddQ2ddF;
}


double Material_Class::Det2D(Eigen::Matrix2d & H) {
	return H(0,0)*H(1,1) - H(0,1)*H(1,0);
}


Eigen::Matrix2d Material_Class::Cofactor2D(Eigen::Matrix2d & H) {
	Eigen::Matrix2d Cofac; Cofac.setZero();
	Cofac(0,0) = H(1,1);
	Cofac(1,0) = -1.0*H(0,1);
	Cofac(0,1) = -1.0*H(1,0);
	Cofac(1,1) = H(0,0);
	return Cofac;
}

Eigen::Matrix2d Material_Class::outer_product(const Eigen::Vector2d & v1, const Eigen::Vector2d & v2) {
	Eigen::Matrix2d Tout; Tout.setZero();
	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			Tout(i,j) = v1[i]*v2[j];
		}
	}
	return Tout;
}

double Material_Class::Tr(Eigen::Matrix2d & H) {
	double sum = 0.0;

	for (unsigned int i = 0; i < 2; i++){
		sum += H(i,i);
	}
	return sum;
}

double Material_Class::Tr2(Eigen::Matrix2d & H){
	double sum = 0.0;
	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++){
			sum += pow(H(i,j),2.0);
		}
	}
	return sum;
}


#endif // MATERIAL_CLASS_CC_

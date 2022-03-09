/*
 * vector.cpp
 *
 *  Created on: Nov 10, 2021
 *      Author: kevin
 */

#include "../include/Vector.h"
namespace HAZEL{

Vector::Vector(){
	values = NULL;

}


Vector::Vector(int Nset){
	resize(Nset);


}


Vector::Vector(const Vector& vold){
	resize(vold.N);
	std::copy(vold.values, vold.values+N, values);
}



void Vector::resize(int Nset){
	values = new double [Nset];
	N = Nset;
}


void Vector::set(int loc, double val){
	values[loc] = val;
}


void Vector::ptrset(int N_,double* vals_){
	resize(N_);
	values = vals_;
}


double Vector::value(int loc){
	return values[loc];
}


int Vector::size(){
	return N;
}


double Vector::norm(){
	/*
	double sqrsum = 0.0;
	for (int i = 0; i < N; i++){
		sqrsum += pow(*(values+i),2);
	}
	return sqrt(sqrsum);
	*/
	return dnrm2_(&N, values, &spacing);
}


void Vector::setzero(){
	for (int i = 0; i < N; i++){
		*(values+i) = 0.0;
	}
}

double * Vector::getptr(){
	return values;
}




double& Vector::operator[](unsigned int i){
	return values[i];
}


const double& Vector::operator[](unsigned int i)const{
	return values[i];
}


void Vector::normalize(){
	double vnorminv = 1.0/norm();
	for (int i = 0; i < N; i++){
		*(values+i) *= vnorminv;
	}
}


void Vector::scale(double d){

	dscal_(&N, &d, values, &spacing);
}


void Vector::set_equal(const Vector & rhs){

	 resize(rhs.N);


	dcopy_(&N,rhs.values, &spacing,values,&spacing);
}


void Vector::Print(){
	std::cout << "-----------------------------" << std::endl;
	for (int i = 0; i < N; i++){
		std::cout << *(values+i) << std::endl;
	}
	std::cout << "-----------------------------" << std::endl;

}


Vector::~Vector() {
	// TODO Auto-generated destructor stub

	if(values){ // True if values is not a null pointer
		delete[] values;
	}


}
}


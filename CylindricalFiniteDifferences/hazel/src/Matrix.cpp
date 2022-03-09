/*
 * Matrix.cpp
 *
 *  Created on: Nov 11, 2021
 *      Author: kevin
 */

#include "../include/Matrix.h"

namespace HAZEL{
Matrix::Matrix() {
	// TODO Auto-generated constructor stub
	values = NULL;

}

Matrix::Matrix(int N_, int M_){
	resize(N_,M_);

}
Matrix::Matrix(const Matrix& vold){
	resize(vold.N,vold.M);
	std::copy(vold.values, vold.values+NM, values);
}

void Matrix::resize(int N_,int M_){
	N = N_;
	M = M_;
	NM = N*M;
	values = new double [NM];
}


void Matrix::set(int i_, int j_, double val_){
	*(values + index(i_,j_)) = val_;

}

void Matrix::set(int i_, double val_){
	*(values + i_) = val_;
}


double Matrix::value(int i_, int j_){
	return values[index(i_,j_)];
}

void Matrix::setzero(){
	for (int i = 0; i < NM; i++){
		*(values+i) = 0.0;
	}
}







void Matrix::set_equal(const Matrix & rhs){

	resize(rhs.N,rhs.M);


	dcopy_(&NM,rhs.values, &spacing,values,&spacing);
}


void Matrix::scale(double d){

	dscal_(&NM, &d, values, &spacing);
}

int Matrix::index(int i_, int j_){
	return i_+ N*j_;
}


int Matrix::rows(){
	return N;
}


int Matrix::cols(){
	return M;
}


double * Matrix::operator[](unsigned int i){

	return (values + N*i);
}


const double * Matrix::operator[](unsigned int i)const{

	return (values + N*i);
}

double Matrix::operator()(unsigned int i_, unsigned int j_){
	return *(values + index(i_,j_));
}



void Matrix::Print(){
	std::cout << "-----------------------------" << std::endl;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			std::cout << *(values + index(i,j)) << ", ";
		}
		std::cout<< std::endl;
	}
	std::cout << "-----------------------------" << std::endl;
}

Matrix::~Matrix() {
	// TODO Auto-generated destructor stub

	if (values){
		delete[] values;
	}
}

}

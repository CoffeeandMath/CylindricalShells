/*
 * Matrix.h
 *
 *  Created on: Nov 11, 2021
 *      Author: kevin
 */

#ifndef HAZEL_MATRIX_H_
#define HAZEL_MATRIX_H_

#include <cblas.h>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <math.h>
#include "Vector.h"


namespace HAZEL{

extern "C"{
extern void dscal_(int *, double *, double *, int *);
}

extern "C"{
extern void dcopy_(int *, double *, int *, double *, int *);
}
extern "C"{
extern void dgemv_(char*, int*, int*, double *, double *, int *, double *, int *, double *, double *, int *);
}

class Matrix {
public:
	Matrix();
	Matrix(int, int);
	Matrix(const Matrix&);
	void resize(int,int);
	void set(int, int, double);
	double value(int,int);
	inline Matrix transpose();
	void Print();
	int rows();
	int cols();
	inline Matrix & add(Matrix &);
	void setzero();
	void scale(double);
	void set_equal(const Matrix &);
	template<typename T>
	inline T vmult(T &);
	template<typename T>
	inline void vmult(T &, T *);


	inline Matrix & operator= (const Matrix &);
	double * operator[](unsigned int);
	const double * operator[](unsigned int)const;
	double operator()(unsigned int, unsigned int);

	virtual ~Matrix();

private:
	int N = 0;
	int M = 0;
	int NM = 0;
	double * values = NULL;;
	int spacing = 1;

	int index(int, int);
	void set(int, double);
};



inline Matrix & Matrix::operator= (const Matrix & rhs)
{
	if (rhs.N != N || rhs.M != M){
		resize(rhs.N,rhs.M);
	}


	std::copy(rhs.values, rhs.values+NM, values);


	return *this;
}

inline Matrix & Matrix::add(Matrix & m1){

	for (int i = 0; i < NM; i++){
		values[i] += *(m1.values+i);
		//vout.set(i,*(values+i) + *(m1.values+i));
	}

	return *this;

}

template<typename T>
inline T Matrix::vmult(T & vect){
	T vout;
	vout.resize(N);
	char ntrans = 'N';
	double one = 1.0;
	int oneint = 1;
	double beta = 0.0;
	dgemv_(&ntrans, &N, &M, &one, values, &N, vect.getptr(), &spacing, &beta, vout.getptr(), &spacing);
	return vout;
}

template<typename T>
inline void Matrix::vmult(T & vect, T * voutptr){

	voutptr->resize(N);
	char ntrans = 'T';
	double one = 1.0;
	int oneint = 1;
	double beta = 0.0;
	dgemv_(&ntrans, &N, &M, &one, values, &N, vect.getptr(), &spacing, &beta, voutptr->getptr(), &spacing);
}

inline Matrix Matrix::transpose(){
	Matrix Mout(M,N);

	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			Mout.set(j,i,values[index(i,j)]);
		}
	}

	return Mout;
}


}

#endif /* HAZEL_MATRIX_H_ */

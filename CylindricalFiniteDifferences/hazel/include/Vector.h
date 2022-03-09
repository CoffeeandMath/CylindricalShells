/*
 * vector.h
 *
 *  Created on: Nov 10, 2021
 *      Author: kevin
 */

#ifndef HAZEL_VECTOR_H_
#define HAZEL_VECTOR_H_


#include <cblas.h>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <math.h>
#include "Matrix.h"







namespace HAZEL{

extern "C" {
extern double ddot_(int *, double*, int *, double*, int *);
}

extern "C"{
extern double dnrm2_(int *, double *, int *);
}

extern "C"{
extern void dscal_(int *, double *, double *, int *);
}

extern "C"{
extern void dcopy_(int *, double *, int *, double *, int *);
}
class Vector {
public:
	Vector();

	Vector(int);
	Vector(const Vector&);
	void resize(int);
	void set(int, double);
	void ptrset(int,double*);
	double value(int);
	int size();
	void Print();
	inline Vector & add(Vector &);
	inline Vector & operator= (const Vector &);
	inline Vector operator+ (const Vector &);
	inline void operator+= (const Vector &);
	inline Vector operator- (const Vector &);
	inline void operator-= (const Vector &);
	inline Vector operator* (double);
	inline double operator* (const Vector &);


	double& operator[](unsigned int);
	const double& operator[](unsigned int)const;
	double norm();
	void normalize();
	void setzero();
	void scale(double);
	void set_equal(const Vector &);
	double * getptr();



	virtual ~Vector();

private:
	int N = 0;
	int spacing = 1;


	double * values = NULL;



};


inline Vector & Vector::add(Vector & v1){

	for (int i = 0; i < N; i++){
		values[i] += *(v1.values +i);
	}
	return *this;
}

inline Vector & Vector::operator= (const Vector & rhs)
{
	if (rhs.N != N){
		values = new double[ rhs.N];
		N = rhs.N;
	}


	std::copy(rhs.values, rhs.values+N, values);


	return *this;
}


inline Vector Vector::operator+ (const Vector & rhs){

	Vector w(rhs.N);
	for (int i = 0; i < N; i++) {
		w.set(i, *(values+i) + *(rhs.values+i));
	}
	return w;

}


inline void Vector::operator+= (const Vector & rhs){
	for (int i = 0; i < N; i++){
		*(values+i) += *(rhs.values+i);
	}
}


inline void Vector::operator-= (const Vector & rhs){
	for (int i = 0; i < N; i++){
		*(values+i) -= *(rhs.values+i);
	}
}


inline Vector Vector::operator- (const Vector & rhs){

	Vector w(rhs.N);
	for (int i = 0; i < N; i++) {
		w.set(i, *(values+i) - *(rhs.values+i));
	}
	return w;

}


inline Vector Vector::operator* (double d){

	Vector w;
	w.set_equal(*this);
	w.scale(d);



	return w;
}

inline double Vector::operator* (const Vector &v){
	/*
	double sum = 0.0;
	for (int i = 0; i < N; i++){
		sum += *(values+i) * (*(v.values+i));
	}
	return sum;
	 */
	return ddot_(&N, values, &spacing, v.values, &spacing);
}


inline Vector operator*(double d,Vector &v){
	return v*d;
}

}


#endif /* HAZEL_VECTOR_H_ */

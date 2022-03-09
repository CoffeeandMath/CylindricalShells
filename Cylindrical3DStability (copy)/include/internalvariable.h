/*
 * internalvariable.h
 *
 *  Created on: Aug 19, 2021
 *      Author: kevin
 */

#ifndef SRC_INTERNALVARIABLE_H_
#define SRC_INTERNALVARIABLE_H_
#include <vector>


template<typename T>
class internal_variable {
public:
	internal_variable(const int,const int);
	internal_variable();
	void set_value(const int,const int,const T& );
	void set_all_values(const T&);
	T get_value(const int, const int) const;
	T* get_pointer(const int, const int);
	void resize(const int, const int);
	int get_n_cells() const;
	int get_n_quad() const;


private:

	std::vector<T> var;
	int n_cells;
	int n_quad;

};
template<typename T>
internal_variable<T>::internal_variable(const int ncells, const int nquad)
: n_cells(ncells), n_quad(n_quad)
{
	var.resize(n_cells*n_quad);

}
template<typename T>
internal_variable<T>::internal_variable()
: n_cells(-1), n_quad(-1)
{

}
template<typename T>
void internal_variable<T>::set_value(const int cell_index, const int q_index, const T& vloc){

	var[cell_index + q_index * n_quad] = vloc;
}
template<typename T>
void internal_variable<T>::set_all_values(const T& vloc){
	for (unsigned int i = 0; i < n_cells; i++){
		for (unsigned int j = 0; j < n_quad; j++){
			var[i + j * n_quad] = vloc;
		}
	}
}
template<typename T>
T internal_variable<T>::get_value(const int cell_index, const int q_index) const{
	return var[cell_index + q_index * n_quad];
}

template<typename T>
T* internal_variable<T>::get_pointer(const int cell_index, const int q_index){
	return &var[cell_index + q_index * n_quad];
}

template<typename T>
void internal_variable<T>::resize(const int ncell, const int nquad ){
	if ((n_cells != ncell) || (nquad != nquad)){
		n_cells = ncell;
		n_quad = nquad;
		var.resize(ncell*nquad);
	}

}

template<typename T>
int internal_variable<T>::get_n_cells() const {
	return n_cells;
}

template<typename T>
int internal_variable<T>::get_n_quad() const {
	return n_quad;
}


#endif /* SRC_INTERNALVARIABLE_H_ */

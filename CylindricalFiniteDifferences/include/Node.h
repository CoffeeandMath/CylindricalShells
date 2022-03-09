/*
 * Node.h
 *
 *  Created on: Jan 24, 2022
 *      Author: kevin
 */

#ifndef CYLINDRICALFINITEDIFFERENCES_SRC_NODE_H_
#define CYLINDRICALFINITEDIFFERENCES_SRC_NODE_H_
#include <vector>
#include <Eigen/Dense>
#include <iostream>

class Node {
public:
	Node();
	Node(double,std::vector<int> &);
	Node(double,std::vector<int> &, Eigen::Vector2d &);
	virtual ~Node();
	void set(Eigen::Vector2d &);
	void set(int,double);
	void setS(double);
	void setfixed(int);
	bool isfree(int);


	Eigen::Vector2d get();
	double get(int);
	double getS();

	std::vector<int> getindex();

private:
	std::vector<bool> is_free;
	Eigen::Vector2d value;
	std::vector<int> index;
	double S = 0.;
};

#endif /* CYLINDRICALFINITEDIFFERENCES_SRC_NODE_H_ */

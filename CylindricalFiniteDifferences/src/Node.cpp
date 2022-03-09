/*
 * Node.cpp
 *
 *  Created on: Jan 24, 2022
 *      Author: kevin
 */

#include "Node.h"

Node::Node() {
	// TODO Auto-generated constructor stub

}

Node::Node(double S_, std::vector<int> & index_){
	S = S_;
	index = index_;
	value.setZero();

	is_free.resize(index_.size());
	for (int i = 0; i < is_free.size(); i++){
		is_free[i] = true;
	}
}

Node::Node(double S_,std::vector<int> & index_, Eigen::Vector2d & values_){
	S = S_;
	index=  index_;
	value = values_;

	is_free.resize(index_.size());
	for (int i = 0; i < is_free.size(); i++){
		is_free[i] = true;
	}


}
void Node::setfixed(int k){
	is_free[k] = false;
}


void Node::set(Eigen::Vector2d & value_){

	value = value_;

}


void Node::set(int i, double value_){
	value[i] = value_;
}

void Node::setS(double S_){

	S = S_;
}


Eigen::Vector2d Node::get(){
	return value;
}

double Node::get(int i){
	return value[i];
}

bool Node::isfree(int j){
	return is_free[j];
}



double Node::getS(){
	return S;
}

std::vector<int> Node::getindex(){
	return index;
}

Node::~Node() {
	// TODO Auto-generated destructor stub
}


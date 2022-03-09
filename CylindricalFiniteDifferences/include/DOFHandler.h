/*
 * DOFHandler.h
 *
 *  Created on: Jan 31, 2022
 *      Author: kevin
 */

#ifndef SRC_DOFHANDLER_H_
#define SRC_DOFHANDLER_H_
#include "Node.h"
#include "Element.h"
#include <Eigen/Sparse>


typedef Eigen::Triplet<double> T;
class DOFHandler {
public:
	DOFHandler();
	DOFHandler(std::vector<Node> *, std::vector<Element> *);

	void apply_solution(const Eigen::VectorXd &);
	Eigen::VectorXd get_solution();
	Eigen::VectorXd condense(Eigen::VectorXd &);
	std::vector<T> condense(std::vector<T> &);
	Eigen::VectorXd distribute(Eigen::VectorXd &);
	Eigen::VectorXd apply(Eigen::VectorXd &);
	std::vector<T> apply(std::vector<T> &);
	int n_free_dofs();
	virtual ~DOFHandler();

private:
	int ndofs = 0;
	int dofspernode = 2;
	std::vector<Node> * NodeP;
	std::vector<Element> * ElemP;

	Eigen::VectorXd solution;
	std::vector<int> globaltonodenumber;
	std::vector<int> globaltonodeindex;
	std::vector<bool> isdofconstrained;
	std::vector<int> globaltoconstrained;
	int nfree = 0;

};

#endif /* SRC_DOFHANDLER_H_ */

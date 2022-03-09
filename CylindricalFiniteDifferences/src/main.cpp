#include <iostream>
#include <Eigen/Dense>
#include "Node.h"
#include "Element.h"
#include "ElasticProblem.h"

using namespace Eigen;

int main()
{

	ElasticProblem SymmetricProblem;
	SymmetricProblem.run();

}

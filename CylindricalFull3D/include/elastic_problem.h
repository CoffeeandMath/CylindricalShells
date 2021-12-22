#ifndef ELASTIC_PROBLEM_H_
#define ELASTIC_PROBLEM_H_

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_direct.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <filesystem>
#include <thread>
#include <future>
#include <chrono>
#include <deal.II/fe/fe_dgq.h>


#include <deal.II/lac/lapack_full_matrix.h>
#include "read_file.h"
#include "internalvariable.h"



//namespace fs = std::filesystem;

#include "material_class.h"
#include "reference_configuration.h"

#define DIM 2
#define pi 3.1415926535897932384626433832795028841971
#define shapefunctiondeg 2
namespace Step4{
using namespace dealii;



class ElasticProblem
{
public:
	ElasticProblem();
	void run();

private:
	void make_grid();
	void setup_system();
	void assemble_system();
	void setup_constraints();
	void initialize_configuration();
	void output_results(const unsigned int) const;
	void newton_raphson();
	void solve();
	void solve_path();
	void calc_displacement();

	Triangulation<DIM> triangulation;
	FESystem<DIM>          fe;
	DoFHandler<DIM>    dof_handler;


	std::vector<bool> create_bool_vector(int,int);


	int quadegadd = 1;
	int refinelevel = 5;
	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	SparseMatrix<double> constraint_matrix;

	AffineConstraints<double> constraints;


	Vector<double> solution;
	Vector<double> system_rhs;
	Vector<double> linearsolve;
	Vector<double> prev_solution;
	Vector<double> displacement;


	std::vector<Material_Class> Material_Vector_InPlane;
	std::vector<Material_Class> Material_Vector_Bending;
	internal_variable<Tensor<2,2>> Applied_Bending;
	internal_variable<Tensor<2,2>> Applied_Strain;

	Tensor<2,2> inverse2x2(const Tensor<2,2> &);
	Tensor<1,3> Normalize(const Tensor<1,3> &);
	Tensor<1,3> Cross(const Tensor<1,3> & , const Tensor<1,3> &);
	Tensor<2,2> Transpose(const Tensor<2,2> &);
	template<int d>
	Tensor<2,d> Symmetric(const Tensor<2,d> &);
	template<int d>
	Tensor<2,d> Outer_Product(const Tensor<1,d> & , const Tensor<1,d>&);
	double Tensor_Inner(const Tensor<2,2> &, const Tensor<2,2> &);
	double BilinearProduct(const Tensor<2,2> &, const Tensor<4,2> &, const Tensor<2,2> &);
	std::vector<double> linspace(double, double, int);


	double h = 0.01;
	double Emodv = 1.e3;
	double nu = 0.3;
	double hsc = pow(h,2)/(12.*(1.-nu*nu));

	double mu = 1.0;


	double tol = 1.e-10;

	double xdrag = 0.e-10;
};
}



#endif /* ELASTIC_PROBLEM_H_ */

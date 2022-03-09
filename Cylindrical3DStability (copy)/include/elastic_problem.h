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

#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include "read_file.h"
#include "internalvariable.h"



//namespace fs = std::filesystem;

#include "material_class.h"
#include "reference_configuration.h"

#define DIM 1
#define pi 3.1415926535897932384626433832795028841971

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
	void setup_constraints();
	std::vector<bool> create_bool_vector(int,int);
	std::vector<double> linspace(double, double, int);
	void initialize_reference_config();
	void update_internal_metrics();
	void assemble_system();
	void assemble_constraint_system();
	void construct_reduced_mappings();
	Vector<double> calculate_stability(int);
	void save_stability(std::vector<std::vector<Vector<double>>>);
	void solve_path();
	void load_external_params();
	static void calc_MMult(Eigen::MatrixXd *,const Eigen::MatrixXd *, const Eigen::MatrixXd *);
	static void calc_MAdd(Eigen::MatrixXd *,const Eigen::MatrixXd *, const Eigen::MatrixXd *,const double a, const double b);
	void calc_Kx2(const Eigen::MatrixXd &);
	void writeToCSVfile(std::string , Eigen::MatrixXd);

	Tensor<2,2> inverse2x2(const Tensor<2,2> &);
	template<int d>
	Tensor<2,d> Symmetric(const Tensor<2,d> &);
	Tensor<1,3> cross_product(const Tensor<1,3> &, const Tensor<1,3> &);
	template<int T>
	Tensor<2,T> outer_product(const Tensor<1,T> &, const Tensor<1,T> &);
	Tensor<1,3> normalize(const Tensor<1,3> &);
	double Tensor_Inner(const Tensor<2,2> &, const Tensor<2,2> &);
	double BilinearProduct(const Tensor<2,2> &, const Tensor<4,2> &, const Tensor<2,2> &);
	void savestate(Vector<double> &, int);
	Triangulation<DIM> triangulation;
	FESystem<DIM>          fe;
	DoFHandler<DIM>    dof_handler;


	void load_state(unsigned int);
	template<typename T>
	T max_val(std::vector<T> &);
	Triangulation<DIM> x0triangulation;
	FESystem<DIM> x0fe;
	DoFHandler<DIM> x0dof_handler;
	Vector<double> x0solution;
	Vector<double> R0solution;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	SparseMatrix<double> constraint_matrix;

	AffineConstraints<double> constraints;


	Vector<double> solution;
	std::vector<std::vector<Reference_Configuration>> Reference_Configuration_Vec;

	std::vector<Material_Class> Material_Vector_InPlane;
	std::vector<Material_Class> Material_Vector_Bending;

	internal_variable<Tensor<2,2>> inplane_a;
	internal_variable<Tensor<2,2>> bending_a;
	Tensor<2,2> current_inplane_a;
	Tensor<2,2> current_bending_a;


	int quadegadd = 4;
	int refinelevel = 7;
	double h = 0.02;
	double nval = 1.0;
	double nu = 0.3;
	double Emodv = 1.0e-3;

	std::vector<int> x_global_to_reduced;
	std::vector<int> xi_global_to_reduced;
	std::vector<int> x_reduced_to_global;
	std::vector<int> xi_reduced_to_global;

	unsigned int nxdofs = 0;
	unsigned int nxidofs = 0;
	int timestep = 0;

	double defstep = 0.0;

	Eigen::MatrixXd XiProj;
	Eigen::MatrixXd XiProjtrans;
	Eigen::MatrixXd Ktilde;

	Eigen::MatrixXd Kx1;
	Eigen::MatrixXd Kx2;
	LAPACKFullMatrix<double> KtildeLA;
	Eigen::MatrixXd Evals;

	read_file rf;
	some_data dat;
};
}



#endif /* ELASTIC_PROBLEM_H_ */

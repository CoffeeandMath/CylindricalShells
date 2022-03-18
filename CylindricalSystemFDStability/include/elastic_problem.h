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
#include <deal.II/fe/fe_dgq.h>
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
#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <filesystem>
#include <Eigen/Dense>
#include <deal.II/lac/lapack_full_matrix.h>
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
	void assemble_system();
	void solve();
	void output_results() const;
	void output_data_csv();
	void output_data_csv_iterative(std::string,int);
	void setup_constraints();
	void newton_raphson();
	void solve_path();
	void initialize_reference_config();
	void update_applied_strains();
	void update_internal_metrics();
	void construct_reduced_mappings();
	void sparse_matrix_to_csv(SparseMatrix<double> &, std::string);
	template<typename T>
	void vector_to_csv(std::vector<T> & , std::string );
	void assemble_constraint_system();
	void extract_and_sort_values();

	Tensor<2,2> calcDDf(const Tensor<1,2> &, const Tensor<1,2> &, const Tensor<1,2> &, int, int,int);
	Tensor<2,6> calcDDf_Fourier(const Tensor<1,2> &, const Tensor<1,2> &, const Tensor<1,2> &, int, int,int,double);
	void calculate_fd_stability();
	void calculate_fd_stability2();
	void calculate_fd_fourier_stability(double);
	void calculate_fd_fourier_stability2(double);
	void output_stability_matrices(int);
	void save_current_state(unsigned int, bool);
	void exportdata(std::string);
	void flip_solution();
	double Tensor_Inner(const Tensor<2,2> &, const Tensor<2,2> &);
	double BilinearProduct(const Tensor<2,2> &, const Tensor<4,2> &, const Tensor<2,2> &);
	Tensor<2,2> transpose2D(Tensor<2,2> &);
	Tensor<2,6> transpose6D(Tensor<2,6> &);
	Tensor<2,2> inverse2D(Tensor<2,2> &);
	Tensor<2,6> inverse6D(Tensor<2,6> &);
	template<int d>
	Tensor<2,d> Symmetric(const Tensor<2,d> & );
	Tensor<1,3> cross_product(const Tensor<1,3> &, const Tensor<1,3> &);
	Vector<double> calculate_stability();
	void save_eigenvalues(std::vector<Vector<double>>);
	void save_matrix(const LAPACKFullMatrix<double> &);


	Triangulation<DIM> triangulation;
	FESystem<DIM>          fe;
	DoFHandler<DIM>    dof_handler;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;
	SparseMatrix<double> constraint_matrix;

	AffineConstraints<double> constraints;
	std::vector<Material_Class> Material_Vector_InPlane;
	std::vector<Material_Class> Material_Vector_Bending;
	std::vector<std::vector<Reference_Configuration>> Reference_Configuration_Vec;
	std::vector<Reference_Configuration> Nodal_Reference_Configuration_Vec;
	std::vector<std::vector<Tensor<2,2>>> epsilon_a;
	std::vector<std::vector<Tensor<2,2>>> b_a;
	std::vector<std::vector<double>> fr;
	std::vector<std::vector<double>> fz;

	Vector<double> solution;
	Vector<double> prev_solution;
	Vector<double> linearsolve;
	Vector<double> system_rhs;

	std::vector<std::pair<double,double>> rvals;
	std::vector<std::pair<double,double>> zvals;
	std::vector<std::pair<double,double>> xirvals;
	std::vector<std::pair<double,double>> xizvals;


	int quadegadd = 10;
	int solveiteration = 0;


	double tol = 1e-15;
	double h = .01;
	double z0 = 0.;
	double r0 = .10;
	double Smax = 1.0;
	int refinelevel = 4;


	double Emodv = 1.0e-1;
	double homog = 0.000;
	double dhomog = 0.000;
	double defmag = 0.0;
	double defmaginplane = 0.0;
	double defmagbending = 0.0;
	double mu = 1000.00 * Emodv;
	double nu = 0.3;

	std::vector<double> linspace(double, double, int);
	Tensor<2,2> inplane_dir;
	Tensor<2,2> bending_dir;
	internal_variable<Tensor<2,2>> inplane_a;
	internal_variable<Tensor<2,2>> bending_a;

	std::vector<int> x_global_to_reduced;
	std::vector<int> xi_global_to_reduced;
	std::vector<int> x_reduced_to_global;
	std::vector<int> xi_reduced_to_global;


	Eigen::MatrixXd XiProj;
	Eigen::MatrixXd XiProjtrans;
	Eigen::MatrixXd Ktilde;

	Eigen::MatrixXd Kx1;
	Eigen::MatrixXd Kx2;
	LAPACKFullMatrix<double> KtildeLA;
	Eigen::MatrixXd Evals;

	unsigned int nxdofs = 0;
	unsigned int nxidofs = 0;
	int N2max = 0;

	int fmodein = 1;

	enum BoundaryConditions{
		Real=0,Testcase,Slopecase,TaylorExpansioncase, TaylorExpansioncase2
	};

};
}



#endif /* ELASTIC_PROBLEM_H_ */

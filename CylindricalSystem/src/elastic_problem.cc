#ifndef ELASTIC_PROBLEM_CC_
#define ELASTIC_PROBLEM_CC_

#include "../include/elastic_problem.h"

#include <deal.II/base/iterator_range.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>
#include <deal.II/dofs/dof_accessor.templates.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_q1.h>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>  // std::setprecision()



namespace Step4
{
using namespace dealii;



ElasticProblem::ElasticProblem()
: fe(FE_Q<DIM>(2), 6)
, dof_handler(triangulation){}


// @sect4{Step4::make_grid}

// Grid creation is something inherently dimension dependent. However, as long
// as the domains are sufficiently similar in 2D or 3D, the library can
// abstract for you. In our case, we would like to again solve on the square
// $[-1,1]\times [-1,1]$ in 2D, or on the cube $[-1,1] \times [-1,1] \times
// [-1,1]$ in 3D; both can be termed GridGenerator::hyper_cube(), so we may
// use the same function in whatever dimension we are. Of course, the
// functions that create a hypercube in two and three dimensions are very much
// different, but that is something you need not care about. Let the library
// handle the difficult things.


void ElasticProblem::solve_path(){

	h = 0.005;
	homog = 0.0;
	dhomog = 0.0;
	//r0 = 1.0;

	bool StabilityCalcs = false;
	initialize_reference_config();

	double defmagmin = 0.00;
	//double defmagmax = 0.6*pi;
	double defmagmax = 4.;
	int Nmax = 400;
	N2max = 800;

	std::vector<double> defmagvec = linspace(defmagmin,defmagmax,Nmax);

	int cntr = 0;
	for (double defmagtemp :defmagvec){
		cntr++;
		defmag=defmagtemp;
		initialize_reference_config();
		update_applied_strains();
		std::cout << "Solve Iteration: " << cntr << "---------------------------" << std::endl;
		newton_raphson();
		//output_data_csv();
	}
	homog = 0.0000;
	dhomog = 0.0000;
	double defmagbendingmin = 0.0;
	double defmagbendingmax = 0.0*2.5*pi;

	std::vector<double> defmagbendingvec = linspace(defmagbendingmin,defmagbendingmax,N2max);
	std::vector<double> defmaginplanevec = linspace(0.0,0.5,N2max);
	inplane_dir[0][0] = 1.0;
	inplane_dir[1][1] = -.30;

	bending_dir[0][0] = -1.0;
	bending_dir[1][1] = 0.0;
	cntr = 0;
	//fs::create_directory("solutions");
	//fs::remove_all("stabilitymatrices");
	//fs::create_directory("stabilitymatrices");



	vector_to_csv(defmagbendingvec, "forward_defmag.csv");

	assemble_constraint_system();
	assemble_system();
	construct_reduced_mappings();
	vector_to_csv(x_global_to_reduced, "x_global_to_reduced.csv");
	vector_to_csv(x_reduced_to_global, "x_reduced_to_global.csv");
	vector_to_csv(xi_global_to_reduced, "xi_global_to_reduced.csv");
	vector_to_csv(xi_reduced_to_global, "xi_reduced_to_global.csv");

	std::vector<Vector<double>> evalsall(N2max);


	for (unsigned int cnt = 0; cnt < N2max; cnt++){
		defmaginplane = defmaginplanevec[cnt];
		defmagbending = defmagbendingvec[cnt];
		std::cout << defmagbending << std::endl;
		initialize_reference_config();
		update_internal_metrics();

		std::cout << "Solve Iteration: " << cnt << "---------------------------" << std::endl;
		newton_raphson();

		if(StabilityCalcs){

			assemble_constraint_system();
			assemble_system();
			construct_reduced_mappings();
			output_stability_matrices(cnt);
			evalsall[cnt] = calculate_stability();
			save_eigenvalues(evalsall);
		}
		output_data_csv_iterative("solutions",cnt);
		save_current_state(cnt, (cnt==0));

	}
	/*
	for (double defmagbendingtemp : defmagbendingvec){
		cntr++;
		defmagbending = defmagbendingtemp;
		initialize_reference_config();
		update_internal_metrics();


		std::cout << "Solve Iteration: " << cntr << "---------------------------" << std::endl;
		newton_raphson();

		if(StabilityCalcs){

			assemble_constraint_system();
			assemble_system();
			construct_reduced_mappings();
			output_stability_matrices(cntr);
			evalsall[cntr-1] = calculate_stability();
			save_eigenvalues(evalsall);
		}
		output_data_csv_iterative("solutions",cntr);
		save_current_state(cntr, (cntr==1));
	}
	 */

	/*

	double defmag3start = defmag2max;
	double defmag3end = defmag2min;
	int N3max = N2max;
	std::vector<double> defmag3vec = linspace(defmag3start,defmag3end,N3max);
	vector_to_csv(defmag3vec, "reverse_defmag.csv");

	for (double defmag3temp : defmag3vec) {
		cntr++;
		defmag2 = defmag3temp;
		initialize_reference_config();
		update_internal_metrics();

		std::cout << "Solve Iteration: " << cntr << "---------------------------" << std::endl;
		newton_raphson();

		if(StabilityCalcs){

			assemble_constraint_system();
			assemble_system();
			construct_reduced_mappings();
			output_stability_matrices(cntr);
			evalsall[cntr-1] = calculate_stability();
			save_eigenvalues(evalsall);
		}
		output_data_csv_iterative("solutions",cntr);
		save_current_state(cntr, (cntr==1));
	}
	 */
	exportdata("params.dat");


}

void ElasticProblem::make_grid()
{
	GridGenerator::hyper_cube(triangulation, 0.0, Smax);
	triangulation.refine_global(refinelevel);

	std::cout << "   Number of active cells: " << triangulation.n_active_cells()
																																																																					<< std::endl << "   Total number of cells: "
																																																																					<< triangulation.n_cells() << std::endl;
}

// @sect4{Step4::setup_system}

// This function looks exactly like in the previous example, although it
// performs actions that in their details are quite different if
// <code>dim</code> happens to be 3. The only significant difference from a
// user's perspective is the number of cells resulting, which is much higher
// in three than in two space dimensions!


void ElasticProblem::initialize_reference_config(){


	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	// We wanted to have a non-constant right hand side, so we use an object of
	// the class declared above to generate the necessary data. Since this right
	// hand side object is only used locally in the present function, we declare
	// it here as a local variable:

	// Compared to the previous example, in order to evaluate the non-constant
	// right hand side function we now also need the quadrature points on the
	// cell we are presently on (previously, we only required values and
	// gradients of the shape function from the FEValues object, as well as the
	// quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
	// do for us by also giving it the #update_quadrature_points flag:
	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	// We then again define the same abbreviation as in the previous program.
	// The value of this variable of course depends on the dimension which we
	// are presently using, but the FiniteElement class does all the necessary
	// work for you and you don't have to care about the dimension dependent
	// parts:
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();



	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();


		Reference_Configuration_Vec[cell_index].resize(n_q_points);
		epsilon_a[cell_index].resize(n_q_points);
		b_a[cell_index].resize(n_q_points);
		fr[cell_index].resize(n_q_points);
		fz[cell_index].resize(n_q_points);
		for (const unsigned int q_index : fe_values.quadrature_point_indices()){
			const auto &x_q = fe_values.quadrature_point(q_index);
			Reference_Configuration_Vec[cell_index][q_index].set_deformation_param(defmag);
			Reference_Configuration_Vec[cell_index][q_index].set_R0(r0);
			Reference_Configuration_Vec[cell_index][q_index].set_point(x_q[0]);

			fr[cell_index][q_index] = 0.0;
			fz[cell_index][q_index] = 0.0;
		}


	}
}

void ElasticProblem::update_internal_metrics(){

	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	// We wanted to have a non-constant right hand side, so we use an object of
	// the class declared above to generate the necessary data. Since this right
	// hand side object is only used locally in the present function, we declare
	// it here as a local variable:

	// Compared to the previous example, in order to evaluate the non-constant
	// right hand side function we now also need the quadrature points on the
	// cell we are presently on (previously, we only required values and
	// gradients of the shape function from the FEValues object, as well as the
	// quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
	// do for us by also giving it the #update_quadrature_points flag:
	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	// We then again define the same abbreviation as in the previous program.
	// The value of this variable of course depends on the dimension which we
	// are presently using, but the FiniteElement class does all the necessary
	// work for you and you don't have to care about the dimension dependent
	// parts:
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();



	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();

		for (const unsigned int q_index : fe_values.quadrature_point_indices()){
			const auto &x_q = fe_values.quadrature_point(q_index);

			Tensor<2,2> epsilontemp;
			Tensor<2,2> btemp;
			double Rval = Reference_Configuration_Vec[cell_index][q_index].get_R();

			inplane_a.set_value(cell_index, q_index, defmaginplane*inplane_dir);
			bending_a.set_value(cell_index,q_index, defmagbending*bending_dir);
		}


	}
}



void ElasticProblem::setup_system()
{
	dof_handler.distribute_dofs(fe);
	solution.reinit(dof_handler.n_dofs());
	prev_solution.reinit(dof_handler.n_dofs());
	linearsolve.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
	Reference_Configuration_Vec.resize(triangulation.n_active_cells());
	Material_Vector_InPlane.resize(triangulation.n_active_cells());
	Material_Vector_Bending.resize(triangulation.n_active_cells());
	epsilon_a.resize(triangulation.n_active_cells());
	b_a.resize(triangulation.n_active_cells());
	fr.resize(triangulation.n_active_cells());
	fz.resize(triangulation.n_active_cells());
	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()<< std::endl;


	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);
	inplane_a.resize(triangulation.n_active_cells(), quadrature_formula.size());
	bending_a.resize(triangulation.n_active_cells(), quadrature_formula.size());

	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler,
			dsp,
			constraints,
			/*keep_constrained_dofs = */ false);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);
	constraint_matrix.reinit(sparsity_pattern);

	std::ofstream out("sparsity_pattern.svg");
	sparsity_pattern.print_svg(out);




}

void ElasticProblem::update_applied_strains(){

	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	// We wanted to have a non-constant right hand side, so we use an object of
	// the class declared above to generate the necessary data. Since this right
	// hand side object is only used locally in the present function, we declare
	// it here as a local variable:

	// Compared to the previous example, in order to evaluate the non-constant
	// right hand side function we now also need the quadrature points on the
	// cell we are presently on (previously, we only required values and
	// gradients of the shape function from the FEValues object, as well as the
	// quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
	// do for us by also giving it the #update_quadrature_points flag:
	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	// We then again define the same abbreviation as in the previous program.
	// The value of this variable of course depends on the dimension which we
	// are presently using, but the FiniteElement class does all the necessary
	// work for you and you don't have to care about the dimension dependent
	// parts:
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();



	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();



		for (const unsigned int q_index : fe_values.quadrature_point_indices()){
			const auto &x_q = fe_values.quadrature_point(q_index);


			fr[cell_index][q_index] = 0.00*defmag;
			fz[cell_index][q_index] = 0.*defmag;
		}

	}
}

// @sect4{Step4::assemble_system}

// Unlike in the previous example, we would now like to use a non-constant
// right hand side function and non-zero boundary values. Both are tasks that
// are readily achieved with only a few new lines of code in the assemblage of
// the matrix and right hand side.
//
// More interesting, though, is the way we assemble matrix and right hand side
// vector DIMension independently: there is simply no difference to the
// two-dimensional case. Since the important objects used in this function
// (quadrature formula, FEValues) depend on the dimension by way of a template
// parameter as well, they can take care of setting up properly everything for
// the dimension for which this function is compiled. By declaring all classes
// which might depend on the dimension using a template parameter, the library
// can make nearly all work for you and you don't have to care about most
// things.

void ElasticProblem::assemble_system()
{
	system_matrix = 0;
	system_rhs    = 0;
	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	// We wanted to have a non-constant right hand side, so we use an object of
	// the class declared above to generate the necessary data. Since this right
	// hand side object is only used locally in the present function, we declare
	// it here as a local variable:

	// Compared to the previous example, in order to evaluate the non-constant
	// right hand side function we now also need the quadrature points on the
	// cell we are presently on (previously, we only required values and
	// gradients of the shape function from the FEValues object, as well as the
	// quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
	// do for us by also giving it the #update_quadrature_points flag:
	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	// We then again define the same abbreviation as in the previous program.
	// The value of this variable of course depends on the dimension which we
	// are presently using, but the FiniteElement class does all the necessary
	// work for you and you don't have to care about the dimension dependent
	// parts:
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);
	std::vector<double> r_q(n_q_points);
	std::vector<double> z_q(n_q_points);
	std::vector<double> lambda_r_q(n_q_points);
	std::vector<double> lambda_z_q(n_q_points);
	std::vector<double> xi_r_q(n_q_points);
	std::vector<double> xi_z_q(n_q_points);

	std::vector<double> r_old_q(n_q_points);
	std::vector<double> z_old_q(n_q_points);

	std::vector<Tensor<1,DIM>> dr_q(n_q_points);
	std::vector<Tensor<1,DIM>> dz_q(n_q_points);

	std::vector<Tensor<1,DIM>> dxi_r_q(n_q_points);
	std::vector<Tensor<1,DIM>> dxi_z_q(n_q_points);

	std::vector<Tensor<1,DIM>> dr_old_q(n_q_points);
	std::vector<Tensor<1,DIM>> dz_old_q(n_q_points);


	std::vector<Tensor<2,DIM>> ddr_q(n_q_points);
	std::vector<Tensor<2,DIM>> ddz_q(n_q_points);



	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	// Next, we again have to loop over all cells and assemble local
	// contributions.  Note, that a cell is a quadrilateral in two space
	// dimensions, but a hexahedron in 3D. In fact, the
	// <code>active_cell_iterator</code> data type is something different,
	// depending on the dimension we are in, but to the outside world they look
	// alike and you will probably never see a difference. In any case, the real
	// type is hidden by using `auto`:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();
		cell_matrix = 0;
		cell_rhs    = 0;

		const FEValuesExtractors::Scalar r(0);
		const FEValuesExtractors::Scalar z(1);
		const FEValuesExtractors::Scalar lambda_r(2);
		const FEValuesExtractors::Scalar lambda_z(3);
		const FEValuesExtractors::Scalar xi_r(4);
		const FEValuesExtractors::Scalar xi_z(5);



		fe_values[r].get_function_values(solution,r_q);
		fe_values[z].get_function_values(solution,z_q);

		fe_values[lambda_r].get_function_values(solution,lambda_r_q);
		fe_values[lambda_z].get_function_values(solution, lambda_z_q);
		fe_values[xi_r].get_function_values(solution, xi_r_q);
		fe_values[xi_z].get_function_values(solution, xi_z_q);




		fe_values[r].get_function_values(prev_solution,r_old_q);
		fe_values[z].get_function_values(prev_solution,z_old_q);

		fe_values[r].get_function_gradients(solution,dr_q);
		fe_values[z].get_function_gradients(solution,dz_q);

		fe_values[xi_r].get_function_gradients(solution,dxi_r_q);
		fe_values[xi_z].get_function_gradients(solution,dxi_z_q);


		// Now we have to assemble the local matrix and right hand side. This is
		// done exactly like in the previous example, but now we revert the
		// order of the loops (which we can safely do since they are independent
		// of each other) and merge the loops for the local matrix and the local
		// vector as far as possible to make things a bit faster.
		//
		// Assembling the right hand side presents the only significant
		// difference to how we did things in step-3: Instead of using a
		// constant right hand side with value 1, we use the object representing
		// the right hand side and evaluate it at the quadrature points:

		for (const unsigned int q_index : fe_values.quadrature_point_indices()){








			Tensor<2,2> CovariantMetric;
			Tensor<2,2> Covariant2Form;

			double hsc = 12.0*pow(h,2.0)/(1.0 - nu*nu);

			double Rq = Reference_Configuration_Vec[cell_index][q_index].get_R();

			double rq = r_q[q_index];
			double rp = dr_q[q_index][0];
			double rpp = dxi_r_q[q_index][0];
			double zp = dz_q[q_index][0];
			double zpp = dxi_z_q[q_index][0];

			const auto &x_q = fe_values.quadrature_point(q_index);

			CovariantMetric[0][0] = 0.5*(pow(rp,2.0) + pow(zp,2.0) - 1.0 );
			CovariantMetric[1][1] = 0.5*(pow(rq/Rq,2.0) - 1.0);

			Covariant2Form[0][0] = (-(rpp*zp) + rp*zpp)/pow(pow(rp,2) + pow(zp,2),2.5);
			Covariant2Form[1][1] = (pow(Rq,2)*zp)/(pow(rq,3)*sqrt(pow(rp,2) + pow(zp,2)));

			Tensor<2,2> InPlane = CovariantMetric  - inplane_a.get_value(cell_index, q_index);
			Tensor<2,2> Bending = Covariant2Form - Reference_Configuration_Vec[cell_index][q_index].get_Covariant_2Form() - bending_a.get_value(cell_index, q_index);

			Material_Vector_InPlane[cell_index].set_Params(Emodv, nu, InPlane);
			Material_Vector_Bending[cell_index].set_Params(Emodv, nu, Bending);


			Tensor<2,4> ddb11mat;
			ddb11mat[0][0] = (5.*(-6.*pow(rp,2)*rpp*zp + rpp*pow(zp,3) + 4.*pow(rp,3)*zpp -
					3.*rp*pow(zp,2)*zpp))/pow(pow(rp,2) + pow(zp,2),4.5);
			ddb11mat[0][1] = (5.*rp*zp)/pow(pow(rp,2) + pow(zp,2),3.5);
			ddb11mat[1][0] = ddb11mat[0][1];
			ddb11mat[0][2] = (5.*(pow(rp,3)*rpp - 6.*rp*rpp*pow(zp,2) + 6.*pow(rp,2)*zp*zpp - pow(zp,3)*zpp))/
					pow(pow(rp,2) + pow(zp,2),4.5);
			ddb11mat[2][0] = ddb11mat[0][2];
			ddb11mat[0][3] = (-4.*pow(rp,2) + pow(zp,2))/pow(pow(rp,2) + pow(zp,2),3.5);
			ddb11mat[3][0] = ddb11mat[0][3];
			ddb11mat[1][2] = (-pow(rp,2) + 4.*pow(zp,2))/pow(pow(rp,2) + pow(zp,2),3.5);
			ddb11mat[2][1] = ddb11mat[1][2];
			ddb11mat[2][2] = (-5.*(-3.*pow(rp,2)*rpp*zp + 4.*rpp*pow(zp,3) + pow(rp,3)*zpp -
					6.*rp*pow(zp,2)*zpp))/pow(pow(rp,2) + pow(zp,2),4.5);
			ddb11mat[2][3] = (-5.*rp*zp)/pow(pow(rp,2) + pow(zp,2),3.5);
			ddb11mat[3][2] = ddb11mat[2][3];

			Tensor<2,3> ddb22mat;
			ddb22mat[0][0] = (12.*pow(Rq,2)*zp)/(pow(rq,5)*sqrt(pow(rp,2) + pow(zp,2)));
			ddb22mat[0][1] = (3.*rp*pow(Rq,2)*zp)/(pow(rq,4)*pow(pow(rp,2) + pow(zp,2),1.5));
			ddb22mat[1][0] = ddb22mat[0][1];
			ddb22mat[0][2] = (-3.*pow(rp,2)*pow(Rq,2))/(pow(rq,4)*pow(pow(rp,2) + pow(zp,2),1.5));
			ddb22mat[2][0] = ddb22mat[0][2];
			ddb22mat[1][1] = -((pow(Rq,2)*zp*(-2*pow(rp,2) + pow(zp,2)))/
					(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),2.5)));
			ddb22mat[1][2] = -((pow(Rq,2)*(pow(rp,3) - 2.*rp*pow(zp,2)))/
					(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),2.5)));
			ddb22mat[2][1] = ddb22mat[1][2];
			ddb22mat[2][2] = (-3.*pow(rp,2)*pow(Rq,2)*zp)/(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),2.5));
			for (const unsigned int i : fe_values.dof_indices())
			{



				const double R_i_q = fe_values[r].value(i,q_index);
				const double Z_i_q = fe_values[z].value(i,q_index);

				const double Lambda_r_i_q = fe_values[lambda_r].value(i,q_index);
				const double Lambda_z_i_q = fe_values[lambda_z].value(i,q_index);

				const double Xi_r_i_q = fe_values[xi_r].value(i,q_index);
				const double Xi_z_i_q = fe_values[xi_z].value(i,q_index);

				const Tensor<1,DIM> dR_i_q = fe_values[r].gradient(i,q_index);
				const Tensor<1,DIM> dZ_i_q = fe_values[z].gradient(i,q_index);

				const Tensor<1,DIM> dXi_r_i_q = fe_values[xi_r].gradient(i,q_index);
				const Tensor<1,DIM> dXi_z_i_q = fe_values[xi_z].gradient(i,q_index);

				double Di_r = R_i_q;
				double Di_rp = dR_i_q[0];
				double Di_rpp = dXi_r_i_q[0];
				double Di_z = Z_i_q;
				double Di_zp = dZ_i_q[0];
				double Di_zpp = dXi_z_i_q[0];



				Tensor<2,2> d_CovariantMetric_i_q;
				d_CovariantMetric_i_q[0][0] = rp*Di_rp + zp*Di_zp;
				d_CovariantMetric_i_q[1][1] = rq*Di_r/pow(Rq,2);

				Tensor<2,2> d_Covariant2Form_i_q;
				d_Covariant2Form_i_q[0][0] = (5.*rp*rpp*zp - 4.*pow(rp,2)*zpp + pow(zp,2)*zpp)/pow(pow(rp,2) + pow(zp,2),3.5) * Di_rp
						-(zp/pow(pow(rp,2) + pow(zp,2),2.5)) * Di_rpp
						+(-(pow(rp,2)*rpp) + 4.*rpp*pow(zp,2) - 5.*rp*zp*zpp)/pow(pow(rp,2) + pow(zp,2),3.5) * Di_zp
						+rp/pow(pow(rp,2) + pow(zp,2),2.5) * Di_zpp;

				d_Covariant2Form_i_q[1][1] = (-3.*pow(Rq,2)*zp)/(pow(rq,4)*sqrt(pow(rp,2) + pow(zp,2)))*Di_r
						-((rp*pow(Rq,2)*zp)/(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),1.5))) * Di_rp
						+(pow(rp,2)*pow(Rq,2))/(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),1.5)) * Di_zp;


				for (const unsigned int j : fe_values.dof_indices()) {

					const double R_j_q = fe_values[r].value(j,q_index);
					const double Z_j_q = fe_values[z].value(j,q_index);

					const double Lambda_r_j_q = fe_values[lambda_r].value(j,q_index);
					const double Lambda_z_j_q = fe_values[lambda_z].value(j,q_index);

					const double Xi_r_j_q = fe_values[xi_r].value(j,q_index);
					const double Xi_z_j_q = fe_values[xi_z].value(j,q_index);

					const Tensor<1,DIM> dR_j_q = fe_values[r].gradient(j,q_index);
					const Tensor<1,DIM> dZ_j_q = fe_values[z].gradient(j,q_index);

					const Tensor<1,DIM> dXi_r_j_q = fe_values[xi_r].gradient(j,q_index);
					const Tensor<1,DIM> dXi_z_j_q = fe_values[xi_z].gradient(j,q_index);

					double Dj_r = R_j_q;
					double Dj_rp = dR_j_q[0];
					double Dj_rpp = dXi_r_j_q[0];
					double Dj_z = Z_j_q;
					double Dj_zp = dZ_j_q[0];
					double Dj_zpp = dXi_z_j_q[0];


					Tensor<2,2> d_CovariantMetric_j_q;
					d_CovariantMetric_j_q[0][0] = rp*Dj_rp + zp*Dj_zp;
					d_CovariantMetric_j_q[1][1] = rq*Dj_r/pow(Rq,2);

					Tensor<2,2> d_Covariant2Form_j_q;
					d_Covariant2Form_j_q[0][0] = (5.*rp*rpp*zp - 4.*pow(rp,2)*zpp + pow(zp,2)*zpp)/pow(pow(rp,2) + pow(zp,2),3.5) * Dj_rp
							-(zp/pow(pow(rp,2) + pow(zp,2),2.5)) * Dj_rpp
							+(-(pow(rp,2)*rpp) + 4.*rpp*pow(zp,2) - 5.*rp*zp*zpp)/pow(pow(rp,2) + pow(zp,2),3.5) * Dj_zp
							+rp/pow(pow(rp,2) + pow(zp,2),2.5) * Dj_zpp;

					d_Covariant2Form_j_q[1][1] = (-3.*pow(Rq,2)*zp)/(pow(rq,4)*sqrt(pow(rp,2) + pow(zp,2)))*Dj_r
							-((rp*pow(Rq,2)*zp)/(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),1.5))) * Dj_rp
							+(pow(rp,2)*pow(Rq,2))/(pow(rq,3)*pow(pow(rp,2) + pow(zp,2),1.5)) * Dj_zp;




					Tensor<2,2> dd_CovariantMetric_ij_q;
					dd_CovariantMetric_ij_q[0][0] = Di_rp*Dj_rp + Di_zp*Dj_zp;
					dd_CovariantMetric_ij_q[1][1] = Di_r*Dj_r/pow(Rq,2);



					Tensor<1,4> di;
					di[0] = Di_rp;
					di[1] = Di_rpp;
					di[2] = Di_zp;
					di[3] = Di_zpp;

					Tensor<1,4> dj;
					dj[0] = Dj_rp;
					dj[1] = Dj_rpp;
					dj[2] = Dj_zp;
					dj[3] = Dj_zpp;

					Tensor<1,3> d2i;
					d2i[0] = Di_r;
					d2i[1] = Di_rp;
					d2i[2] = Di_zp;

					Tensor<1,3> d2j;
					d2j[0] = Dj_r;
					d2j[1] = Dj_rp;
					d2j[2] = Dj_zp;



					Tensor<2,2> dd_Covariant2Form_ij_q;
					dd_Covariant2Form_ij_q[0][0] = dj*ddb11mat*di;
					dd_Covariant2Form_ij_q[1][1] = d2j*ddb22mat*d2i;
					/*
					cell_matrix(i,j) += (dR_i_q[0]*dR_j_q[0] + dZ_i_q[0]*dZ_j_q[0])*fe_values.JxW(q_index);
					cell_matrix(i,j) += 100.0*(R_i_q*R_j_q + Z_i_q*Z_j_q)*fe_values.JxW(q_index);
					 */
					///*

					//Contribution of the in plane stretching
					cell_matrix(i,j) += Rq*( BilinearProduct(d_CovariantMetric_i_q,Material_Vector_InPlane[cell_index].getddQ2ddF(),d_CovariantMetric_j_q) )*fe_values.JxW(q_index);
					cell_matrix(i,j) += Rq*(Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),dd_CovariantMetric_ij_q))*fe_values.JxW(q_index);

					//Contribution of the bending component
					cell_matrix(i,j) += hsc*Rq*( BilinearProduct(d_Covariant2Form_i_q,Material_Vector_Bending[cell_index].getddQ2ddF(),d_Covariant2Form_j_q) )*fe_values.JxW(q_index);
					cell_matrix(i,j) += hsc*Rq*(Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),dd_Covariant2Form_ij_q))*fe_values.JxW(q_index);



					//Homogenization/drag term
					cell_matrix(i,j) += homog*Rq*(R_i_q*R_j_q + Z_i_q*Z_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += dhomog*Rq*(dR_i_q[0]*dR_j_q[0] + dZ_i_q[0]*dZ_j_q[0])*fe_values.JxW(q_index);



					// Augmented Lagrangian terms
					cell_matrix(i,j) -= Lambda_r_i_q*(dR_j_q[0] - Xi_r_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) -= Lambda_r_j_q*(dR_i_q[0] - Xi_r_i_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) -= Lambda_z_i_q*(dZ_j_q[0] - Xi_z_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) -= Lambda_z_j_q*(dZ_i_q[0] - Xi_z_i_q)*fe_values.JxW(q_index);

					cell_matrix(i,j) += mu*(dR_i_q[0] - Xi_r_i_q)*(dR_j_q[0] - Xi_r_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += mu*(dZ_i_q[0] - Xi_z_i_q)*(dZ_j_q[0] - Xi_z_j_q)*fe_values.JxW(q_index);

					//*/
				}


				/*
				cell_rhs(i) += ((dr_q[q_index][0]*dR_i_q[0] + dz_q[q_index][0]*dZ_i_q[0]) * // phi_i(x_q)
						fe_values.JxW(q_index));            // dx

				cell_rhs(i) += 100.0*((r_q[q_index]-0.5)*R_i_q + (z_q[q_index] - 0.75)*Z_i_q)*fe_values.JxW(q_index);
				 */

				// /*

				cell_rhs(i) += (Rq*(Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),d_CovariantMetric_i_q)))*fe_values.JxW(q_index);

				cell_rhs(i) += hsc*(Rq*(Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),d_Covariant2Form_i_q)))*fe_values.JxW(q_index);


				cell_rhs(i) += homog*Rq*((r_q[q_index] - r_old_q[q_index])*R_i_q + (z_q[q_index] - z_old_q[q_index])*Z_i_q)*fe_values.JxW(q_index);
				cell_rhs(i) += dhomog*Rq*((dr_q[q_index][0] -dr_old_q[q_index][0])*dR_i_q[0] + (dz_q[q_index][0] -dz_old_q[q_index][0])*dZ_i_q[0])*fe_values.JxW(q_index);

				cell_rhs(i) -= Rq*(fr[cell_index][q_index]*R_i_q + fz[cell_index][q_index]*Z_i_q)*fe_values.JxW(q_index);


				// Augmented Lagruangian terms

				cell_rhs(i) -= lambda_r_q[q_index]*(dR_i_q[0] - Xi_r_i_q)*fe_values.JxW(q_index);
				cell_rhs(i) -= Lambda_r_i_q*(dr_q[q_index][0] - xi_r_q[q_index])*fe_values.JxW(q_index);
				cell_rhs(i) -= lambda_z_q[q_index]*(dZ_i_q[0] - Xi_z_i_q)*fe_values.JxW(q_index);
				cell_rhs(i) -= Lambda_z_i_q*(dz_q[q_index][0] - xi_z_q[q_index])*fe_values.JxW(q_index);

				cell_rhs(i) += mu*(dR_i_q[0] - Xi_r_i_q)*(dr_q[q_index][0] - xi_r_q[q_index])*fe_values.JxW(q_index);
				cell_rhs(i) += mu*(dZ_i_q[0] - Xi_z_i_q)*(dz_q[q_index][0] - xi_z_q[q_index])*fe_values.JxW(q_index);

				//cell_rhs(i) += hscinv*homog*Z_i_q*(z_q[q_index] - x_q[0])*fe_values.JxW(q_index);
				// */




			}
		}

		cell->get_dof_indices(local_dof_indices);

		//constraints.distribute_local_to_global( cell_matrix, cell_rhs, local_dof_indices,system_matrix,system_rhs);

		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < local_dof_indices.size(); i++) {

			system_rhs[local_dof_indices[i]] += cell_rhs[i];
			for (unsigned int j = 0; j < local_dof_indices.size(); j++) {
				system_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[i][j]);
			}
		}
	}

	constraints.condense(system_matrix);
	constraints.condense(system_rhs);
}


void ElasticProblem::assemble_constraint_system()
{
	constraint_matrix = 0;

	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	// We wanted to have a non-constant right hand side, so we use an object of
	// the class declared above to generate the necessary data. Since this right
	// hand side object is only used locally in the present function, we declare
	// it here as a local variable:

	// Compared to the previous example, in order to evaluate the non-constant
	// right hand side function we now also need the quadrature points on the
	// cell we are presently on (previously, we only required values and
	// gradients of the shape function from the FEValues object, as well as the
	// quadrature weights, FEValues::JxW() ). We can tell the FEValues object to
	// do for us by also giving it the #update_quadrature_points flag:
	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	// We then again define the same abbreviation as in the previous program.
	// The value of this variable of course depends on the dimension which we
	// are presently using, but the FiniteElement class does all the necessary
	// work for you and you don't have to care about the dimension dependent
	// parts:
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);
	std::vector<double> r_q(n_q_points);
	std::vector<double> z_q(n_q_points);
	std::vector<double> lambda_r_q(n_q_points);
	std::vector<double> lambda_z_q(n_q_points);
	std::vector<double> xi_r_q(n_q_points);
	std::vector<double> xi_z_q(n_q_points);

	std::vector<double> r_old_q(n_q_points);
	std::vector<double> z_old_q(n_q_points);

	std::vector<Tensor<1,DIM>> dr_q(n_q_points);
	std::vector<Tensor<1,DIM>> dz_q(n_q_points);

	std::vector<Tensor<1,DIM>> dxi_r_q(n_q_points);
	std::vector<Tensor<1,DIM>> dxi_z_q(n_q_points);

	std::vector<Tensor<1,DIM>> dr_old_q(n_q_points);
	std::vector<Tensor<1,DIM>> dz_old_q(n_q_points);


	std::vector<Tensor<2,DIM>> ddr_q(n_q_points);
	std::vector<Tensor<2,DIM>> ddz_q(n_q_points);



	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	// Next, we again have to loop over all cells and assemble local
	// contributions.  Note, that a cell is a quadrilateral in two space
	// dimensions, but a hexahedron in 3D. In fact, the
	// <code>active_cell_iterator</code> data type is something different,
	// depending on the dimension we are in, but to the outside world they look
	// alike and you will probably never see a difference. In any case, the real
	// type is hidden by using `auto`:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();
		cell_matrix = 0;
		cell_rhs    = 0;

		const FEValuesExtractors::Scalar r(0);
		const FEValuesExtractors::Scalar z(1);

		const FEValuesExtractors::Scalar xi_r(4);
		const FEValuesExtractors::Scalar xi_z(5);



		fe_values[r].get_function_values(solution,r_q);
		fe_values[z].get_function_values(solution,z_q);


		fe_values[xi_r].get_function_values(solution, xi_r_q);
		fe_values[xi_z].get_function_values(solution, xi_z_q);


		fe_values[r].get_function_gradients(solution,dr_q);
		fe_values[z].get_function_gradients(solution,dz_q);

		fe_values[xi_r].get_function_gradients(solution,dxi_r_q);
		fe_values[xi_z].get_function_gradients(solution,dxi_z_q);


		// Now we have to assemble the local matrix and right hand side. This is
		// done exactly like in the previous example, but now we revert the
		// order of the loops (which we can safely do since they are independent
		// of each other) and merge the loops for the local matrix and the local
		// vector as far as possible to make things a bit faster.
		//
		// Assembling the right hand side presents the only significant
		// difference to how we did things in step-3: Instead of using a
		// constant right hand side with value 1, we use the object representing
		// the right hand side and evaluate it at the quadrature points:

		for (const unsigned int q_index : fe_values.quadrature_point_indices()){


			double R_ref_q = Reference_Configuration_Vec[cell_index][q_index].get_R();


			const auto &x_q = fe_values.quadrature_point(q_index);


			for (const unsigned int i : fe_values.dof_indices())
			{
				const double R_i_q = fe_values[r].value(i,q_index);
				const double Z_i_q = fe_values[z].value(i,q_index);


				const double Xi_r_i_q = fe_values[xi_r].value(i,q_index);
				const double Xi_z_i_q = fe_values[xi_z].value(i,q_index);

				const Tensor<1,DIM> dR_i_q = fe_values[r].gradient(i,q_index);
				const Tensor<1,DIM> dZ_i_q = fe_values[z].gradient(i,q_index);

				const Tensor<1,DIM> dXi_r_i_q = fe_values[xi_r].gradient(i,q_index);
				const Tensor<1,DIM> dXi_z_i_q = fe_values[xi_z].gradient(i,q_index);




				Tensor<2,2> d_CovariantMetric_i_q;



				for (const unsigned int j : fe_values.dof_indices()) {

					const double R_j_q = fe_values[r].value(j,q_index);
					const double Z_j_q = fe_values[z].value(j,q_index);


					const double Xi_r_j_q = fe_values[xi_r].value(j,q_index);
					const double Xi_z_j_q = fe_values[xi_z].value(j,q_index);

					const Tensor<1,DIM> dR_j_q = fe_values[r].gradient(j,q_index);
					const Tensor<1,DIM> dZ_j_q = fe_values[z].gradient(j,q_index);

					const Tensor<1,DIM> dXi_r_j_q = fe_values[xi_r].gradient(j,q_index);
					const Tensor<1,DIM> dXi_z_j_q = fe_values[xi_z].gradient(j,q_index);





					//cell_matrix(i,j) += mu*(dR_i_q[0] - Xi_r_i_q)*(dR_j_q[0] - Xi_r_j_q)*fe_values.JxW(q_index);
					//cell_matrix(i,j) += mu*(dZ_i_q[0] - Xi_z_i_q)*(dZ_j_q[0] - Xi_z_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += 100000.0*R_ref_q*(dR_i_q[0] - Xi_r_i_q)*(Xi_r_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += 100000.0*R_ref_q*(dZ_i_q[0] - Xi_z_i_q)*(Xi_z_j_q)*fe_values.JxW(q_index);
					//*/
				}






			}
		}

		cell->get_dof_indices(local_dof_indices);

		//constraints.distribute_local_to_global( cell_matrix, cell_rhs, local_dof_indices,system_matrix,system_rhs);

		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < local_dof_indices.size(); i++) {
			for (unsigned int j = 0; j < local_dof_indices.size(); j++) {
				constraint_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[i][j]);
			}
		}
	}

	constraints.condense(constraint_matrix);
}

// @sect4{Step4::solve}




// I want to write a new function to deal with the boundary conditions as constraints

void ElasticProblem::flip_solution(){


	//DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	//VectorTools::interpolate_boundary_values(dof_handler,
	//		1,
	//		Functions::ZeroFunction<DIM>(DIM+5),
	//		constraints);


	const int ndofs = dof_handler.n_dofs();
	// Constraint stuff
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);

	std::vector<bool> xir_components = {false,false,false,false,true,false};
	ComponentMask xir_mask(r_components);

	std::vector<bool> xiz_components = {false,false,false,false,false,true};
	ComponentMask xiz_mask(z_components);


	std::vector<bool> is_r_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> is_z_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);

	std::vector<bool> is_xir_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xir_mask, is_xir_comp);

	std::vector<bool> is_xiz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xiz_mask, is_xiz_comp);

	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);


	for (unsigned int i = 0; i < ndofs; i++) {

		///*

		if (is_z_comp[i] && is_xiz_comp[i]){
			solution[i] = -1.0*solution[i];
		}


		//*/


	}
	std::cout<< "Flipped Up-Down" << std::endl;



}
void ElasticProblem::setup_constraints(){
	constraints.clear();

	//DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	//VectorTools::interpolate_boundary_values(dof_handler,
	//		1,
	//		Functions::ZeroFunction<DIM>(DIM+5),
	//		constraints);


	const int ndofs = dof_handler.n_dofs();
	// Constraint stuff
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);


	std::vector<bool> is_r_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> is_z_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);


	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);


	for (unsigned int i = 0; i < ndofs; i++) {

		///*


		if (fabs(support_points[i][0]) < 1.0e-8) {
			std::cout << support_points[i] << std::endl;
			if (is_r_comp[i]) {
				//constraints.add_line(i);
				solution[i] = r0;
			} else if (is_z_comp[i]) {
				constraints.add_line(i);
				std::cout << i << std::endl;
				solution[i] = z0;
			}



			std::cout << "Constraint index: " << i << std::endl;
		}
		//*/


	}
	std::cout<< "Filled out constraints" << std::endl;

	constraints.close();


	for (unsigned int i = 0; i < ndofs; i++) {
		if (is_r_comp[i]) {
			solution[i] = r0;
		} else if (is_z_comp[i]) {
			solution[i] = support_points[i][0];
		}
	}
	prev_solution = solution;

}

void ElasticProblem::output_data_csv( ){
	//constraints.clear();

	//DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	//VectorTools::interpolate_boundary_values(dof_handler,
	//		1,
	//		Functions::ZeroFunction<DIM>(DIM+5),
	//		constraints);


	const int ndofs = dof_handler.n_dofs();
	// Constraint stuff
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);


	std::vector<bool> is_r_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> is_z_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);


	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	std::vector<double> rSvals;
	std::vector<double> rvals;

	std::vector<double> zSvals;
	std::vector<double> zvals;

	for (unsigned int i = 0; i < ndofs; i++) {

		if (is_r_comp[i]){
			rSvals.push_back(support_points[i][0]);
			rvals.push_back(solution[i]);
		} else if (is_z_comp[i]) {
			zSvals.push_back(support_points[i][0]);
			zvals.push_back(solution[i]);
		}
	}


	std::string rname = "r_values.csv";
	std::ofstream rfile;
	rfile.open(rname);
	rfile << "S_values,r_values\n";
	for (unsigned int i = 0; i < rvals.size(); i++){
		rfile << rSvals[i] << "," << rvals[i] << "\n";
	}
	rfile.close();

	std::string zname = "z_values.csv";
	std::ofstream zfile;
	zfile.open(zname);
	zfile << "S_values,z_values\n";
	for (unsigned int i = 0; i < zvals.size(); i++){
		zfile << zSvals[i] << "," << zvals[i] << "\n";
	}
	zfile.close();




}


void ElasticProblem::output_data_csv_iterative(std::string foldername,int iter){
	//constraints.clear();

	//DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	//VectorTools::interpolate_boundary_values(dof_handler,
	//		1,
	//		Functions::ZeroFunction<DIM>(DIM+5),
	//		constraints);


	const int ndofs = dof_handler.n_dofs();
	// Constraint stuff
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);

	std::vector<bool> xi_r_components = {false,false,false,false,true,false};
	ComponentMask xi_r_mask(xi_r_components);

	std::vector<bool> xi_z_components = {false,false,false,false,false,true};
	ComponentMask xi_z_mask(xi_z_components);


	std::vector<bool> is_r_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> is_z_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);

	std::vector<bool> is_xi_r_comp(ndofs,false);
	DoFTools::extract_dofs(dof_handler, xi_r_mask, is_xi_r_comp);

	std::vector<bool> is_xi_z_comp(ndofs,false);
	DoFTools::extract_dofs(dof_handler, xi_z_mask, is_xi_z_comp);



	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	std::vector<double> rSvals;
	std::vector<double> rvals;

	std::vector<double> zSvals;
	std::vector<double> zvals;

	std::vector<double> xiRSvals;
	std::vector<double> xiRvals;

	std::vector<double> xiZSvals;
	std::vector<double> xiZvals;

	for (unsigned int i = 0; i < ndofs; i++) {

		if (is_r_comp[i]){
			rSvals.push_back(support_points[i][0]);
			rvals.push_back(solution[i]);
		} else if (is_z_comp[i]) {
			zSvals.push_back(support_points[i][0]);
			zvals.push_back(solution[i]);
		} else if (is_xi_r_comp[i]) {
			xiRSvals.push_back(support_points[i][0]);
			xiRvals.push_back(solution[i]);
		} else if (is_xi_z_comp[i]){
			xiZSvals.push_back(support_points[i][0]);
			xiZvals.push_back(solution[i]);
		}
	}


	std::string rname = foldername + "/r_values_" + std::to_string(iter) + ".csv";
	std::ofstream rfile;
	rfile.open(rname);
	rfile << "S_values,r_values\n";
	for (unsigned int i = 0; i < rvals.size(); i++){
		rfile << rSvals[i] << "," << rvals[i] << "\n";
	}
	rfile.close();

	std::string zname = foldername + "/z_values_" + std::to_string(iter) + ".csv";
	std::ofstream zfile;
	zfile.open(zname);
	zfile << "S_values,z_values\n";
	for (unsigned int i = 0; i < zvals.size(); i++){
		zfile << zSvals[i] << "," << zvals[i] << "\n";
	}
	zfile.close();

	std::string xirname = foldername + "/xir_values_" + std::to_string(iter) + ".csv";
	std::ofstream xirfile;
	xirfile.open(xirname);
	xirfile << "S_values,xir_values\n";
	for (unsigned int i = 0; i < xiRvals.size(); i++){
		xirfile << xiRSvals[i] << "," << xiRvals[i] << "\n";
	}
	xirfile.close();

	std::string xizname = foldername + "/xiz_values_" + std::to_string(iter) + ".csv";
	std::ofstream xizfile;
	xizfile.open(xizname);
	xizfile << "S_values,xiz_values\n";
	for (unsigned int i = 0; i < xiZvals.size(); i++){
		xizfile << xiZSvals[i] << "," << xiZvals[i] << "\n";
	}
	xizfile.close();


}

void ElasticProblem::exportdata(std::string outname){
	std::ofstream rfile;
	rfile.open(outname);

	rfile << h << '\n';
	rfile << refinelevel << '\n';
	rfile << N2max << '\n';
	rfile << nu << '\n';
	rfile << defmagbending << '\n';
	rfile << bending_dir[0][0] << '\n';
	rfile << bending_dir[0][1] << '\n';
	rfile << bending_dir[1][1] << '\n';
	rfile << defmaginplane << '\n';
	rfile << inplane_dir[0][0] << '\n';
	rfile << inplane_dir[0][1] << '\n';
	rfile << inplane_dir[1][1] << '\n';
	rfile.close();
}

// Solving the linear system of equations is something that looks almost
// identical in most programs. In particular, it is dimension independent, so
// this function is copied verbatim from the previous example.
void ElasticProblem::solve()
{



	constraints.distribute(linearsolve);
	linearsolve = 0;
	SparseDirectUMFPACK a_direct;
	a_direct.initialize(system_matrix);
	a_direct.vmult(linearsolve,system_rhs);

	constraints.distribute(linearsolve);
	solution.add(-1.0,linearsolve);


}




void ElasticProblem::newton_raphson() {
	double stepsize = 2.0*tol;
	int cntr = 0;
	prev_solution = solution;
	while (stepsize > tol) {
		//prev_solution = solution;
		cntr++;

		assemble_system();

		solve();


		stepsize = sqrt(linearsolve.norm_sqr())/linearsolve.size();
		//residual = step_direction;

		//std::cout << "Iteration: " << cntr << std::endl;
		//std::cout << "Step Size: " << stepsize<< std::endl;
	}
	std::cout << "Number of Steps: " << cntr << std::endl;
}


double ElasticProblem::Tensor_Inner(const Tensor<2,2> & H1, const Tensor<2,2> & H2) {
	double sum = 0.0;
	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			sum += H1[i][j]*H2[i][j];
		}
	}
	return sum;
}

double ElasticProblem::BilinearProduct(const Tensor<2,2> & F1, const Tensor<4,2> & C, const Tensor<2,2> & F2) {
	double sum = 0.0;

	for (unsigned int i = 0; i < 2; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			for (unsigned int k = 0; k < 2; k++) {
				for (unsigned int l = 0; l < 2; l++) {
					sum += F1[i][j]*C[i][j][k][l]*F2[k][l];
				}
			}
		}
	}
	return sum;
}


void ElasticProblem::construct_reduced_mappings(){


	x_global_to_reduced.clear();
	x_reduced_to_global.clear();
	xi_global_to_reduced.clear();
	xi_reduced_to_global.clear();


	const int ndofs = dof_handler.n_dofs();

	x_global_to_reduced.resize(ndofs);
	xi_global_to_reduced.resize(ndofs);
	// Constraint stuff
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);

	std::vector<bool> xi_r_components = {false,false,false,false,true,false};
	ComponentMask xi_r_mask(xi_r_components);

	std::vector<bool> xi_z_components = {false,false,false,false,false,true};
	ComponentMask xi_z_mask(xi_z_components);


	std::vector<bool> is_r_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> is_z_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);

	std::vector<bool> is_xi_r_comp(ndofs,false);
	DoFTools::extract_dofs(dof_handler, xi_r_mask, is_xi_r_comp);

	std::vector<bool> is_xi_z_comp(ndofs,false);
	DoFTools::extract_dofs(dof_handler, xi_z_mask, is_xi_z_comp);


	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	unsigned int xcntr = 0;
	unsigned int xicntr = 0;
	nxdofs = 0;
	nxidofs = 0;

	for (unsigned int i = 0; i < ndofs; i++) {

		if (is_r_comp[i] || is_z_comp[i] ){
			x_reduced_to_global.push_back(i);
			x_global_to_reduced[i] = xcntr;
			xcntr++;
			nxdofs++;
		} else {
			x_global_to_reduced[i] = -1;
		}

		if (is_xi_r_comp[i] || is_xi_z_comp[i]){
			xi_reduced_to_global.push_back(i);
			xi_global_to_reduced[i] = xicntr;
			xicntr++;
			nxidofs++;
		} else {
			xi_global_to_reduced[i] = -1;
		}

	}


	KtildeLA.reinit(nxdofs,nxdofs);
	KtildeLA = 0.0;


	Eigen::MatrixXd Kxx(nxdofs,nxdofs); Kxx.setZero();
	Eigen::MatrixXd Kxxi(nxdofs,nxidofs); Kxxi.setZero();
	Eigen::MatrixXd Kxixi(nxidofs,nxidofs); Kxixi.setZero();
	Eigen::MatrixXd Kxix(nxidofs,nxdofs); Kxix.setZero();
	Eigen::MatrixXd Cxxi(nxdofs,nxidofs); Cxxi.setZero();
	Eigen::MatrixXd Cxixi(nxidofs,nxidofs); Cxixi.setZero();

	std::cout << "Constructing Reduced Matrices" << std::endl;
	for (unsigned int row = 0; row < system_matrix.m(); ++row)
	{
		const typename SparseMatrix<double>::const_iterator end_row = system_matrix.end(row);
		int global_row = row;
		for (typename SparseMatrix<double>::const_iterator entry = system_matrix.begin(row);
				entry != end_row; ++entry)
		{

			int global_column = entry->column();
			if(fabs(entry->value()) > 1e-15 ){

				//       system_matrix_petsc.set(row, entry->column(),entry->value());
				if(x_global_to_reduced[global_row] > -1 && x_global_to_reduced[global_column] > -1) {
					Kxx(x_global_to_reduced[global_row] ,x_global_to_reduced[global_column]) += entry->value();
				} else if (x_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){
					Kxxi(x_global_to_reduced[global_row] ,xi_global_to_reduced[global_column]) += entry->value();
				} else if (xi_global_to_reduced[global_row] > -1 && x_global_to_reduced[global_column] > -1){
					Kxix(xi_global_to_reduced[global_row] ,x_global_to_reduced[global_column]) += entry->value();
				} else if (xi_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){
					Kxixi(xi_global_to_reduced[global_row] ,xi_global_to_reduced[global_column]) += entry->value();
				}
			}

		}
	}

	if (XiProj.cols() < 1){
		std::cout << "Constructing XiR" << std::endl;
		for (unsigned int row = 0; row < constraint_matrix.m(); ++row)
		{
			const typename SparseMatrix<double>::const_iterator end_row = constraint_matrix.end(row);
			int global_row = row;
			for (typename SparseMatrix<double>::const_iterator entry = constraint_matrix.begin(row);
					entry != end_row; ++entry)
			{

				int global_column = entry->column();
				if(fabs(entry->value()) > 1e-15 ){
					//       system_matrix_petsc.set(row, entry->column(),entry->value());

					if (x_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){

						Cxxi(x_global_to_reduced[global_row] ,xi_global_to_reduced[global_column]) += entry->value();

					} else if (xi_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){
						Cxixi(xi_global_to_reduced[global_row] ,xi_global_to_reduced[global_column]) += entry->value();
					}
				}

			}
		}

		Eigen::MatrixXd Cxix = Cxxi.transpose();


		XiProj = -1.0*Cxixi.inverse()*Cxix;

		for (unsigned int i = 0; i < XiProj.rows(); i++){
			for (unsigned int j = 0; j < XiProj.cols(); j++){
				if (fabs(XiProj(i,j)) < 1e-15){
					XiProj(i,j) = 0.0;
				}
			}
		}
		XiProjtrans = XiProj.transpose();

	}
	std::cout << "Calculating Ktilde" << std::endl;

	Eigen::MatrixXd Ktildeasym = Kxx + Kxxi*XiProj + XiProjtrans*Kxix + XiProjtrans*Kxixi*XiProj;

	Ktilde.setZero();
	Ktilde = 0.5*(Ktildeasym + Ktildeasym.transpose());


	for (unsigned int i = 0; i < Ktilde.rows(); i++) {
		for (unsigned int j = 0; j < Ktilde.cols(); j++){
			if (fabs(Ktilde(i,j)) < 1e-15){
				Ktilde(i,j) = 0.0;
			} else {
				KtildeLA(i,j) = Ktilde(i,j);

			}
		}
	}



}

Vector<double> ElasticProblem::calculate_stability(){
	Vector<double> eigenvalues;
	FullMatrix<double> eigenvectors;

	double tol = 1.0e-10;
	Vector<double> evalsout(10);
	KtildeLA.compute_eigenvalues_symmetric(-1.0*Emodv, 1000.0*h*Emodv, tol, eigenvalues, eigenvectors);

	for (unsigned int i = 0; i < 10; i++){
		std::cout << eigenvalues[i] << std::endl;
		evalsout[i] = eigenvalues[i]/Emodv;
	}
	//std::cout << "Smallest Eigenvalue : " << eigenvalues[0] << std::endl;


	return evalsout;


}

void ElasticProblem::save_eigenvalues(std::vector<Vector<double>> evout){
	std::string strout = "evals.csv";
	std::ofstream rfile;
	rfile.open(strout);

	for (unsigned int i = 0; i < evout.size(); i++){
		for (unsigned int j = 0; j < evout[i].size(); j++){
			rfile << std::setprecision(17) << evout[i][j];
			if (j < (evout[i].size() - 1)){
				rfile << ',';
			}
		}
		rfile << '\n';
	}

	rfile.close();
}


void ElasticProblem::sparse_matrix_to_csv(SparseMatrix<double> & SPMat, std::string strout){
	std::ofstream rfile;
	rfile.open(strout);


	for (unsigned int row = 0; row < SPMat.m(); ++row)
	{
		const typename SparseMatrix<double>::const_iterator end_row = SPMat.end(row);
		for (typename SparseMatrix<double>::const_iterator entry = SPMat.begin(row);
				entry != end_row; ++entry)
		{
			if(fabs(entry->value()) > 1e-12 ){
				//       system_matrix_petsc.set(row, entry->column(),entry->value());

				rfile << std::setprecision(17) << row << "," << entry->column() << "," << entry->value() << '\n';
			}

		}
	}
	rfile.close();
}

template<typename T>
void ElasticProblem::vector_to_csv(std::vector<T> & vec, std::string strout){
	std::ofstream rfile;
	rfile.open(strout);

	for (unsigned int i = 0; i < vec.size(); ++i){
		rfile << vec[i] << '\n';
	}


	rfile.close();
}



void ElasticProblem::output_stability_matrices(int cntr){
	std::string n1 = "stabilitymatrices/x_global_to_reduced" + std::to_string(cntr) + ".csv";
	std::string n2 = "stabilitymatrices/x_reduced_to_global" + std::to_string(cntr) + ".csv";
	std::string n3 = "stabilitymatrices/xi_global_to_reduced" + std::to_string(cntr) + ".csv";
	std::string n4 = "stabilitymatrices/xi_reduced_to_global" + std::to_string(cntr) + ".csv";
	std::string n5 = "stabilitymatrices/spout" + std::to_string(cntr) + ".csv";
	std::string n6 = "stabilitymatrices/spconstraint" + std::to_string(cntr) + ".csv";
	vector_to_csv(x_global_to_reduced, n1);
	vector_to_csv(x_reduced_to_global, n2);
	vector_to_csv(xi_global_to_reduced, n3);
	vector_to_csv(xi_reduced_to_global, n4);
	sparse_matrix_to_csv(system_matrix, n5);
	sparse_matrix_to_csv(constraint_matrix,n6);
}



std::vector<double> ElasticProblem::linspace(double start_in, double end_in, int num_in)
{

	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (end - start) / (num - 1);

	for(int i=0; i < num-1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
	// are exactly the same as the input
	return linspaced;
}
// @sect4{Step4::output_results}

// This function also does what the respective one did in step-3. No changes
// here for dimension independence either.
//
// Since the program will run both 2d and 3d versions of the Laplace solver,
// we use the dimension in the filename to generate distinct filenames for
// each run (in a better program, one would check whether <code>dim</code> can
// have other values than 2 or 3, but we neglect this here for the sake of
// brevity).

void ElasticProblem::output_results() const
{
	DataOut<DIM> data_out;

	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, "solution");

	data_out.build_patches();

	std::ofstream output(DIM == 2 ? "solution-2d.vtk" : "solution-1d.vtk");
	data_out.write_vtk(output);
}

void ElasticProblem::save_current_state(unsigned int indx, bool firstTime){


	std::vector<Point<DIM>> support_points(dof_handler.n_dofs());
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
	std::vector<bool> r_components = {true,false,false,false,false,false};
	ComponentMask r_mask(r_components);

	std::vector<bool> is_r_comp(dof_handler.n_dofs(), false);
	DoFTools::extract_dofs(dof_handler, r_mask, is_r_comp);

	std::vector<bool> z_components = {false,true,false,false,false,false};
	ComponentMask z_mask(z_components);

	std::vector<bool> is_z_comp(dof_handler.n_dofs(), false);
	DoFTools::extract_dofs(dof_handler, z_mask, is_z_comp);

	std::vector<bool> lr_components = {false,false,true,false,false,false};
	ComponentMask lr_mask(lr_components);

	std::vector<bool> is_lr_comp(dof_handler.n_dofs(), false);
	DoFTools::extract_dofs(dof_handler, lr_mask, is_lr_comp);

	if (firstTime==true){
		std::string dof_file_name = "de/dof_handler.dofdat";
		std::ofstream dof_out_u(dof_file_name);
		boost::archive::text_oarchive dof_ar_u(dof_out_u);
		dof_handler.save(dof_ar_u,1);

		std::string triang_file_name = "de/triang.tridat";
		std::ofstream triang_out_u(triang_file_name);
		boost::archive::text_oarchive triang_ar_u(triang_out_u);
		triangulation.save(triang_ar_u,1);

		Vector<double> R0vec;
		R0vec.reinit(solution.size());


		for (unsigned int i = 0; i < R0vec.size(); i++){

			if (is_r_comp[i]) {
				Reference_Configuration Refobj;
				Refobj.set_deformation_param(defmag);
				Refobj.set_R0(r0);
				Refobj.set_point(support_points[i][0]);

				R0vec[i] = Refobj.get_R();

			} else if (is_z_comp[i]) {
				Reference_Configuration Refobj;
				Refobj.set_deformation_param(defmag);
				Refobj.set_R0(r0);
				Refobj.set_point(support_points[i][0]);

				R0vec[i] = Refobj.get_phi();

			} else if (is_lr_comp[i]){
				Reference_Configuration Refobj;
				Refobj.set_deformation_param(defmag);
				Refobj.set_R0(r0);
				Refobj.set_point(support_points[i][0]);

				R0vec[i] = Refobj.get_dphi();
			}

		}


		std::string R0_file_name = "de/R0.soldat";
		std::ofstream R0_out_u(R0_file_name);
		boost::archive::text_oarchive R0_ar_u(R0_out_u);
		R0vec.save(R0_ar_u,1);

	}


	std::string sol_file_name = "de/solution_" + std::to_string(indx) + ".soldat";
	std::ofstream sol_out_u(sol_file_name);
	boost::archive::text_oarchive sol_ar_u(sol_out_u);
	solution.save(sol_ar_u,1);



}


// @sect4{Step4::run}

// This is the function which has the top-level control over everything. Apart
// from one line of additional output, it is the same as for the previous
// example.

void ElasticProblem::run()
{
	std::cout << "Solving problem in " << DIM << " space dimensions."
			<< std::endl;

	make_grid();

	setup_system();
	setup_constraints();
	initialize_reference_config();
	//assemble_system();
	//solve();
	//newton_raphson();
	solve_path();
	output_results();
	output_data_csv();
}


}

#endif // ELASTIC_PROBLEM_CC_

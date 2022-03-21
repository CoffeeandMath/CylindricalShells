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
#include<bits/stdc++.h>




namespace Step4
{
using namespace dealii;



ElasticProblem::ElasticProblem()
: fe(FESystem<DIM>(FE_Q<DIM>(2), 2),1,FESystem<DIM>(FE_DGQ<DIM>(0), 2),1,FESystem<DIM>(FE_Q<DIM>(2), 2),1)
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

	h = 0.01;
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
		solveiteration = cnt;
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

		if (fabs(fmodein) < 1e-3){
			calculate_fd_stability();
		} else {
			calculate_fd_fourier_stability((double) fmodein);
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


	Nodal_Reference_Configuration_Vec.resize(rvals.size());
	extract_and_sort_values();
	for (int i = 0; i < rvals.size(); i++){
		Nodal_Reference_Configuration_Vec[i].set_deformation_param(defmag);
		Nodal_Reference_Configuration_Vec[i].set_R0(r0);
		Nodal_Reference_Configuration_Vec[i].set_point(rvals[i].first );
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

	extract_and_sort_values();




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


Tensor<2,2> ElasticProblem::calcDDf(const Tensor<1,2> & u, const Tensor<1,2> & up, const Tensor<1,2> & upp, int dir1, int dir2,int index){

	Tensor<2,2> ddf;
	double R_ref_q = Nodal_Reference_Configuration_Vec[index].get_R();


	Tensor<2,2> CovariantMetric;
	Tensor<2,2> Covariant2Form;

	double hsc = 12.0*pow(h,2.0)/(1.0 - nu*nu);

	double Rq = R_ref_q;

	double rq = u[0];
	double rp = up[0];
	double rpp = upp[0];
	double zp = up[1];
	double zpp = upp[1];



	CovariantMetric[0][0] = 0.5*(pow(rp,2.0) + pow(zp,2.0) - 1.0 );
	CovariantMetric[1][1] = 0.5*(pow(rq/Rq,2.0) - 1.0);

	Covariant2Form[0][0] = (-(rpp*zp) + rp*zpp)/pow(pow(rp,2) + pow(zp,2),2.5);
	Covariant2Form[1][1] = (pow(Rq,2)*zp)/(pow(rq,3)*sqrt(pow(rp,2) + pow(zp,2)));

	Tensor<2,2> InPlane = CovariantMetric  - inplane_a.get_value(10, 1);
	Tensor<2,2> Bending = Covariant2Form - Nodal_Reference_Configuration_Vec[index].get_Covariant_2Form() - bending_a.get_value(10, 1);

	Material_Vector_InPlane[0].set_Params(Emodv, nu, InPlane);
	Material_Vector_Bending[0].set_Params(Emodv, nu, Bending);


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



	for (int i = 0; i < 2; i++){
		Tensor<1,2> Du_i;
		Tensor<1,2> Dup_i;
		Tensor<1,2> Dupp_i;
		if (dir1 == 0) {
			Du_i[i] = 1.;
		} else if (dir1 == 1) {
			Dup_i[i] = 1.;
		} else if (dir1 == 2){
			Dupp_i[i] = 1.;
		}





		double Di_r = Du_i[0];
		double Di_rp = Dup_i[0];
		double Di_rpp = Dupp_i[0];
		double Di_z = Du_i[1];
		double Di_zp = Dup_i[1];
		double Di_zpp = Dupp_i[1];



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




		for (int j = 0; j < 2; j++){
			Tensor<1,2> Du_j;
			Tensor<1,2> Dup_j;
			Tensor<1,2> Dupp_j;
			if (dir2 == 0) {
				Du_j[j] = 1.;
			} else if (dir2 == 1) {
				Dup_j[j] = 1.;
			} else if (dir2 == 2){
				Dupp_j[j] = 1.;
			}




			double Dj_r = Du_j[0];
			double Dj_rp = Dup_j[0];
			double Dj_rpp = Dupp_j[0];
			double Dj_z = Du_j[1];
			double Dj_zp = Dup_j[1];
			double Dj_zpp = Dupp_j[1];


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
			dd_Covariant2Form_ij_q[0][0] = dj*(ddb11mat*di);
			dd_Covariant2Form_ij_q[1][1] = d2j*(ddb22mat*d2i);

			//Contribution of the in plane stretching
			double val = 0.;
			val += Rq*( BilinearProduct(d_CovariantMetric_i_q,Material_Vector_InPlane[0].getddQ2ddF(),d_CovariantMetric_j_q) );
			val += Rq*(Tensor_Inner(Material_Vector_InPlane[0].getdQ2dF(),dd_CovariantMetric_ij_q));

			//Contribution of the bending component
			val += hsc*Rq*( BilinearProduct(d_Covariant2Form_i_q,Material_Vector_Bending[0].getddQ2ddF(),d_Covariant2Form_j_q) );
			val += hsc*Rq*(Tensor_Inner(Material_Vector_Bending[0].getdQ2dF(),dd_Covariant2Form_ij_q));
			ddf[i][j] = val;
		}


	}


	return ddf;


}

void ElasticProblem::calculate_fd_stability(){
	extract_and_sort_values();


	const int Ndof = 2;
	std::vector<double> Si(rvals.size());
	std::vector<Tensor<1,2>> u(rvals.size());
	std::vector<Tensor<1,2>> up(rvals.size());
	std::vector<Tensor<1,2>> upp(rvals.size());

	std::vector<Tensor<2,Ndof>> f00(rvals.size());
	std::vector<Tensor<2,Ndof>> f01(rvals.size());
	std::vector<Tensor<2,Ndof>> f02(rvals.size());
	std::vector<Tensor<2,Ndof>> f11(rvals.size());
	std::vector<Tensor<2,Ndof>> f12(rvals.size());
	std::vector<Tensor<2,Ndof>> f22(rvals.size());


	for (int i = 0; i < rvals.size(); i++){
		Si[i] = rvals[i].first;
		u[i][0] = rvals[i].second;
		u[i][1] = zvals[i].second;
		up[i][0] = xirvals[i].second;
		up[i][1] = xizvals[i].second;


	}



	double dl = Si[1] - Si[0];
	//upp[0] = (up[1] - up[0])/dl;
	//upp[0] = -11.*up[0]/(6. * dl) + 3.*up[1]/dl - 1.5 * up[2]/dl + up[3]/(3.*dl);

	for (int i = 1; i < (rvals.size()-1); i++){
		upp[i] = (up[i+1] - up[i-1])/(2.*dl);

	}




	upp[0] = 2.5*upp[1] - 2.*upp[2] + 0.5*upp[3];
	upp[rvals.size()-1] = 2.5*upp[rvals.size()-2] - 2.*upp[rvals.size()-3] + 0.5*upp[rvals.size()-4];
	/*
	for (int i = 0; i < upp.size(); i++){
		std::cout << upp[i] << std::endl;
	}
	 */

	for (int i = 0; i < rvals.size(); i++){
		f00[i] = calcDDf(u[i], up[i], upp[i], 0, 0, i);
		f01[i] = calcDDf(u[i], up[i], upp[i], 0, 1, i);
		f02[i] = calcDDf(u[i], up[i], upp[i], 0, 2, i);
		f11[i] = calcDDf(u[i], up[i], upp[i], 1, 1, i);
		f12[i] = calcDDf(u[i], up[i], upp[i], 1, 2, i);
		f22[i] = calcDDf(u[i], up[i], upp[i], 2, 2, i);

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;


	}
	Tensor<2,Ndof> Zero;
	int zeropad = 1;
	for (int i = 0; i < zeropad; i++){
		f00[i] = Zero;
		f01[i] = Zero;
		f02[i] = Zero;
		f11[i] = Zero;
		f12[i] = Zero;
		f22[i] = Zero;


		f00[rvals.size() - i - 1] = Zero;
		f01[rvals.size() - i - 1] = Zero;
		f02[rvals.size() - i - 1] = Zero;
		f11[rvals.size() - i - 1] = Zero;
		f12[rvals.size() - i - 1] = Zero;
		f22[rvals.size() - i - 1] = Zero;

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;
	}



	Tensor<2,Ndof> Iden;
	for (int i = 0; i < Ndof; i++){
		Iden[i][i] = 1.;
	}
	// We will also need f11 defined on the half nodes
	std::vector<Tensor<2,Ndof>> f11ext(Ndof*rvals.size()-1);

	for (int i = 0; i < f11ext.size(); i++){
		if (i%2==0){
			f11ext[i] = f11[i/2];
		} else {

			f11ext[i] = calcDDf(0.5*(u[i/2]+u[i/2 + 1]),0.5*(up[i/2]+up[i/2 + 1]),0.5*(upp[i/2]+upp[i/2 + 1]),1,1,i/2);
		}


	}

	for (int i = 0; i < 2*zeropad; i++){
		f11ext[i] = Zero;
		f11ext[f11ext.size() - i - 1] = Zero;
	}




	//I now want to construct the local hessians. I will deal with boundary nodes later
	// the structure will be a vector of tensors at every node, where the vector is the local
	// constructed ones

	std::vector<std::vector<Tensor<2,Ndof>>> DDf22wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> DDf21wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> DDf20w(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df12wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df11wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df10w(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f02wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f01wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f00w(rvals.size());

	std::vector<std::vector<Tensor<2,Ndof>>> localhessian(rvals.size());

	//initializing the stencil
	int stencilsize = 5;
	for (int i = 0; i < rvals.size(); i++){
		DDf22wpp[i].resize(stencilsize);
		DDf21wp[i].resize(stencilsize);
		DDf20w[i].resize(stencilsize);
		Df12wpp[i].resize(stencilsize);
		Df11wp[i].resize(stencilsize);
		Df10w[i].resize(stencilsize);
		f02wpp[i].resize(stencilsize);
		f01wp[i].resize(stencilsize);
		f00w[i].resize(stencilsize);

		localhessian[i].resize(stencilsize);
	}


	//setting the values for interior nodes
	for (int i = 1; i < (rvals.size() - 1); i++){

		DDf22wpp[i][0] = f22[i-1]/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][2] = (f22[i+1] + 4.*f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22[i+1] + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22[i+1]/pow(dl,4);


		DDf21wp[i][0] = -0.5*(transpose(f12[i-1])/pow(dl,3));
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f12[i-1] - f12[i+1];
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(transpose(f12[i+1]))/pow(dl,3);


		DDf20w[i][1] = transpose(f02[i-1])/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = transpose(f02[i+1])/pow(dl,2);


		Df12wpp[i][0] = -0.5*f12[i-1]/pow(dl,3);
		Df12wpp[i][1] = f12[i-1]/pow(dl,3);
		Df12wpp[i][2] = 0.5*(f12[i+1]-f12[i-1])/pow(dl,3);
		Df12wpp[i][3] = -f12[i+1]/pow(dl,3);
		Df12wpp[i][4] = 0.5*f12[i+1]/pow(dl,3);

		Df11wp[i][1] = f11ext[2*i-1]/pow(dl,2);
		Df11wp[i][2] = -(f11ext[2*i-1] + f11ext[2*i+1])/pow(dl,2);
		Df11wp[i][3] = f11ext[2*i+1]/pow(dl,2);





		Df10w[i][1] = -0.5*transpose(f01[i-1])/dl;
		Df10w[i][3] = 0.5*transpose(f01[i+1])/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];

	}


	//Dealing with boundary conditions. for now I will use both bcs from integration by parts
	//beginning with 0 boundary

	Tensor<2,Ndof> f20m1 = 0.*(2.*transpose(f02[0]) - transpose(f02[1]));
	Tensor<2,Ndof> f21m1 = 0.*(2.*transpose(f12[0]) - transpose(f12[1]));
	Tensor<2,Ndof> f11m1 = 0.*(2.*f11[0] - f11[1]);
	Tensor<2,Ndof> f21m1_2 = 0.*(0.5*(f21m1 + transpose(f12[0])));
	Tensor<2,Ndof> f21p1_2 = 0.*calcDDf(0.5*(u[0]+u[1]),0.5*(up[0]+up[1]),0.5*(upp[0]+upp[1]),2,1,1);
	Tensor<2,Ndof> f22m1 = 0.*(2.*f22[0] - f22[1]);
	Tensor<2,Ndof> f01m1 = 0.*(2.*f01[0] - f01[1]);



	std::vector<Tensor<2,Ndof>> Df21wp0(3);
	Df21wp0[0] = f21m1_2/pow(dl,2);
	Df21wp0[1] = -1.0*(f21m1_2 + f21p1_2)/pow(dl,2);
	Df21wp0[2] = f21p1_2/pow(dl,2);

	std::vector<Tensor<2,Ndof>> Df20w0(3);
	Df20w0[0] = -.5*f20m1/dl;
	Df20w0[2] = .5*transpose(f02[1])/dl;

	std::vector<Tensor<2,Ndof>> Df22wpp0(5);
	Df22wpp0[0] = -0.5*f22m1/pow(dl,3);
	Df22wpp0[1] = f22m1/pow(dl,3);
	Df22wpp0[2] = 0.5*(f22[1] - f22m1)/pow(dl,3);
	Df22wpp0[3] = -f22[1]/pow(dl,3);
	Df22wpp0[4] = 0.5*f22[1]/pow(dl,3);






	{
		int i = 0;

		DDf22wpp[i][0] = f22m1/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22m1)/pow(dl,4);
		DDf22wpp[i][2] = (f22[i+1] + 4.*f22[i] + f22m1)/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22[i+1] + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22[i+1]/pow(dl,4);


		DDf21wp[i][0] = -0.5*f21m1/pow(dl,3);
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f21m1 - f12[i+1];
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(transpose(f12[i+1]))/pow(dl,3);


		DDf20w[i][1] = f20m1/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = transpose(f02[i+1])/pow(dl,2);


		Df12wpp[i][0] = -0.5*f21m1/pow(dl,3);
		Df12wpp[i][1] = f21m1/pow(dl,3);
		Df12wpp[i][2] = 0.5*(f12[i+1]-f21m1)/pow(dl,3);
		Df12wpp[i][3] = -f12[i+1]/pow(dl,3);
		Df12wpp[i][4] = 0.5*f12[i+1]/pow(dl,3);


		Df11wp[i][1] = f11m1/pow(dl,2);
		Df11wp[i][2] = -(f11m1 + f11ext[2*i+1])/pow(dl,2);
		Df11wp[i][3] = f11ext[2*i+1]/pow(dl,2);


		Df10w[i][1] = -0.5*transpose(f01m1)/dl;
		Df10w[i][3] = 0.5*transpose(f01[i+1])/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];
	}



	//Going to the L boundary boundary
	int Nend = rvals.size();


	Tensor<2,Ndof> f20Np1 = 0.*(2.*transpose(f02[Nend-1]) - transpose(f02[Nend-2]));
	Tensor<2,Ndof> f21Np1 = 0.*(2.*transpose(f12[Nend-1]) - transpose(f12[Nend-2]));
	Tensor<2,Ndof> f11Np1 = 0.*(2.*f11[Nend-1] - f11[Nend-2]);
	Tensor<2,Ndof> f21Np1_2 = 0.*(0.5*(f21Np1 + transpose(f12[Nend-1])));
	Tensor<2,Ndof> f21Nm1_2 = 0.*calcDDf(0.5*(u[Nend-2]+u[Nend-1]),0.5*(up[Nend-2]+up[Nend-1]),0.5*(upp[Nend-2]+upp[Nend-1]),2,1,Nend-1);
	Tensor<2,Ndof> f22Np1 = 0.*(2.*f22[Nend-1] - f22[Nend-2]);
	Tensor<2,Ndof> f01Np1 = 0.*(2.*f01[Nend-1] - f01[Nend-2]);



	//finding el equations for the N boundaries
	{
		int i = Nend-1;

		DDf22wpp[i][0] = f22[i-1]/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][2] = (f22Np1 + 4.*f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22Np1 + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22Np1/pow(dl,4);


		//this one has a problem
		DDf21wp[i][0] = -0.5*(transpose(f12[i-1])/pow(dl,3));
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f12[i-1] - transpose(f21Np1);
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(f21Np1)/pow(dl,3);



		DDf20w[i][1] = transpose(f02[i-1])/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = f20Np1/pow(dl,2);


		Df12wpp[i][0] = -0.5*f12[i-1]/pow(dl,3);
		Df12wpp[i][1] = f12[i-1]/pow(dl,3);
		Df12wpp[i][2] = 0.5*(transpose(f21Np1)-f12[i-1])/pow(dl,3);
		Df12wpp[i][3] = -transpose(f21Np1)/pow(dl,3);
		Df12wpp[i][4] = 0.5*transpose(f21Np1)/pow(dl,3);

		Df11wp[i][1] = f11ext[2*i-1]/pow(dl,2);
		Tensor<2,Ndof> f11Np1_2 = 0.0*0.5*(f11[i] + f11Np1);
		Df11wp[i][2] = -(f11ext[2*i-1] + f11Np1_2)/pow(dl,2);
		Df11wp[i][3] = f11Np1_2/pow(dl,2);


		Df10w[i][1] = -0.5*transpose(f01[i-1])/dl;
		Df10w[i][3] = 0.5*transpose(f01Np1)/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];


	}




	for (int i = 0; i < localhessian.size(); i++){
		for (int j = 0; j < stencilsize; j++){
			localhessian[i][j] = DDf22wpp[i][j] + DDf21wp[i][j] + DDf20w[i][j] - Df12wpp[i][j] - Df11wp[i][j] - Df10w[i][j] + f02wpp[i][j] + f01wp[i][j] + f00w[i][j];

		}
	}




	LAPACKFullMatrix<double> FDHessian;


	FDHessian.reinit(Ndof*(rvals.size()+4),Ndof*(rvals.size()+4));
	FDHessian = 0.;



	for (int i = 0; i < (rvals.size()); i++){
		for (int j = 0; j < stencilsize; j++){
			for (int k = 0; k < Ndof; k++){
				for (int l = 0; l < Ndof; l++){
					FDHessian(Ndof*(i+2)+k,Ndof*(i+j)+l) = localhessian[i][j][k][l];
				}
			}
		}
	}

	save_matrix(FDHessian);


}

void ElasticProblem::calculate_fd_stability2(){
	extract_and_sort_values();


	const int Ndof = 2;
	std::vector<double> Si(rvals.size());
	std::vector<Tensor<1,2>> u(rvals.size());
	std::vector<Tensor<1,2>> up(rvals.size());
	std::vector<Tensor<1,2>> upp(rvals.size());

	std::vector<Tensor<2,Ndof>> f00(rvals.size());
	std::vector<Tensor<2,Ndof>> f01(rvals.size());
	std::vector<Tensor<2,Ndof>> f02(rvals.size());
	std::vector<Tensor<2,Ndof>> f11(rvals.size());
	std::vector<Tensor<2,Ndof>> f12(rvals.size());
	std::vector<Tensor<2,Ndof>> f22(rvals.size());


	for (int i = 0; i < rvals.size(); i++){
		Si[i] = rvals[i].first;
		u[i][0] = rvals[i].second;
		u[i][1] = zvals[i].second;
		up[i][0] = xirvals[i].second;
		up[i][1] = xizvals[i].second;


	}

	std::vector<Tensor<1,2>> utilde(rvals.size()-1);
	std::vector<Tensor<1,2>> utilde_p(rvals.size()-1);
	std::vector<Tensor<1,2>> utilde_pp(rvals.size()-1);

	std::vector<Tensor<2,Ndof>> f00tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f01tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f02tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f11tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f12tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f22tilde(rvals.size()-1);

	double dl = Si[1] - Si[0];
	//upp[0] = (up[1] - up[0])/dl;
	//upp[0] = -11.*up[0]/(6. * dl) + 3.*up[1]/dl - 1.5 * up[2]/dl + up[3]/(3.*dl);

	for (int i = 1; i < (rvals.size()-1); i++){
		upp[i] = (up[i+1] - up[i-1])/(2.*dl);

	}




	upp[0] = 2.5*upp[1] - 2.*upp[2] + 0.5*upp[3];
	upp[rvals.size()-1] = 2.5*upp[rvals.size()-2] - 2.*upp[rvals.size()-3] + 0.5*upp[rvals.size()-4];
	/*
	for (int i = 0; i < upp.size(); i++){
		std::cout << upp[i] << std::endl;
	}
	 */

	for (int i = 0; i < rvals.size(); i++){
		f00[i] = calcDDf(u[i], up[i], upp[i], 0, 0, i);
		f01[i] = calcDDf(u[i], up[i], upp[i], 0, 1, i);
		f02[i] = calcDDf(u[i], up[i], upp[i], 0, 2, i);
		f11[i] = calcDDf(u[i], up[i], upp[i], 1, 1, i);
		f12[i] = calcDDf(u[i], up[i], upp[i], 1, 2, i);
		f22[i] = calcDDf(u[i], up[i], upp[i], 2, 2, i);

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;


	}

	//setting half values

	for (int i = 0; i < utilde.size(); i++){
		utilde[i] = 0.5*(u[i] + u[i+1]);
		utilde_p[i] = 0.5*(up[i] + up[i+1]);
		utilde_pp[i] = 0.5*(upp[i] + upp[i+1]);

		f00tilde[i] = 0.5*(f00[i]+f00[i+1]);
		f01tilde[i] = 0.5*(f01[i]+f01[i+1]);
		f02tilde[i] = 0.5*(f02[i]+f02[i+1]);
		f11tilde[i] = 0.5*(f11[i]+f11[i+1]);
		f12tilde[i] = 0.5*(f12[i]+f12[i+1]);
		f22tilde[i] = 0.5*(f22[i]+f22[i+1]);

	}




	Tensor<2,Ndof> Zero;


	Tensor<2,Ndof> Iden;
	for (int i = 0; i < Ndof; i++){
		Iden[i][i] = 1.;
	}
	// We will also need f11 defined on the half nodes


	int Nend = rvals.size();


	std::vector<std::vector<Tensor<2,Ndof>>> Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> D2i(rvals.size()-1);

	int StencilSize2 = 4;
	for (int i = 0; i < (rvals.size()-1); i++){
		Gi[i].resize(StencilSize2);
		Di[i].resize(StencilSize2);
		D2i[i].resize(StencilSize2);
	}

	for (int i = 0; i < Gi.size(); i++){
		Gi[i][0] = Zero;
		Gi[i][1] = 0.0*Iden;
		Gi[i][2] = 1.0*Iden;
		Gi[i][3] = Zero;

		Di[i][0] = Zero;
		Di[i][1] = -Iden/dl;
		Di[i][2] = Iden/dl;
		Di[i][3] = Zero;

		if (i == 0) {
			D2i[i][0] = Zero;
			D2i[i][1] = Iden/pow(dl,2);
			D2i[i][2] = -2.*Iden/pow(dl,2);
			D2i[i][3] = Iden/pow(dl,2);

		} else if (i == (Gi.size()-1)) {

			D2i[i][0] = Iden/pow(dl,2);
			D2i[i][1] = -2.*Iden/pow(dl,2);
			D2i[i][2] = Iden/pow(dl,2);
			D2i[i][3] = Zero;

		} else {
			/*
			D2i[i][0] = Iden/(2.*pow(dl,2));
			D2i[i][1] = -Iden/(2.*pow(dl,2));
			D2i[i][2] = -Iden/(2.*pow(dl,2));
			D2i[i][3] = Iden/(2.*pow(dl,2));
			 */

			D2i[i][0] = Iden/(pow(dl,2));
			D2i[i][1] = -2.*Iden/pow(dl,2);
			D2i[i][2] = Iden/pow(dl,2);
			D2i[i][3] = Zero;
		}
	}


	std::vector<std::vector<Tensor<2,Ndof>>> f22D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f21Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f20Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f12D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f11Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f10Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f02D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f01Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f00Gi(rvals.size()-1);

	for (int i = 0; i < (rvals.size()-1); i++){
		f22D2i[i].resize(StencilSize2);
		f21Di[i].resize(StencilSize2);
		f20Gi[i].resize(StencilSize2);
		f12D2i[i].resize(StencilSize2);
		f11Di[i].resize(StencilSize2);
		f10Gi[i].resize(StencilSize2);
		f02D2i[i].resize(StencilSize2);
		f01Di[i].resize(StencilSize2);
		f00Gi[i].resize(StencilSize2);

		for (int j = 0; j < StencilSize2; j++){
			f22D2i[i][j] = f22tilde[i]*D2i[i][j];
			f21Di[i][j] = transpose(f12tilde[i])*Di[i][j];
			f20Gi[i][j] = transpose(f02tilde[i])*Gi[i][j];

			f12D2i[i][j] = f12tilde[i]*D2i[i][j];
			f11Di[i][j] = f11tilde[i]*Di[i][j];
			f10Gi[i][j] = transpose(f01tilde[i])*Gi[i][j];

			f02D2i[i][j]= f02tilde[i]*D2i[i][j];
			f01Di[i][j] = f01tilde[i]*Di[i][j];
			f00Gi[i][j] = f00tilde[i]*Gi[i][j];

		}
	}


	LAPACKFullMatrix<double> G; G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); G = 0.;
	LAPACKFullMatrix<double> D; D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); D = 0.;
	LAPACKFullMatrix<double> D2; D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); D2 = 0.;
	LAPACKFullMatrix<double> f22D2; f22D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f22D2 = 0.;
	LAPACKFullMatrix<double> f21D; f21D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f21D = 0.;
	LAPACKFullMatrix<double> f20G; f20G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f20G = 0.;
	LAPACKFullMatrix<double> f12D2; f12D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f12D2 = 0.;
	LAPACKFullMatrix<double> f11D; f11D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f11D = 0.;
	LAPACKFullMatrix<double> f10G; f10G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f10G = 0.;
	LAPACKFullMatrix<double> f02D2; f02D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f02D2 = 0.;
	LAPACKFullMatrix<double> f01D; f01D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f01D = 0.;
	LAPACKFullMatrix<double> f00G; f00G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f00G = 0.;


	for (int i = 0; i < Di.size(); i++){
		if (i==0){
			for (int j = 1; j < StencilSize2; j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];

					}
				}
			}

		} else if (i == (Di.size()-1)){
			for (int j = 0; j < StencilSize2-1; j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];
					}
				}
			}
		} else {
			for (int j = 0; j < (StencilSize2); j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];
					}
				}
			}

		}
	}

	LAPACKFullMatrix<double> Hout; Hout.reinit(Ndof*rvals.size(),Ndof*rvals.size()); Hout = 0.;
	LAPACKFullMatrix<double> F2; F2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F2 = 0.;
	LAPACKFullMatrix<double> F1; F1.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F1 = 0.;
	LAPACKFullMatrix<double> F0; F0.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F0 = 0.;

	for (int i = 0; i < f22D2.m(); i++){
		for (int j = 0; j < f22D2.n(); j++){
			F2(i,j) = f22D2(i,j) + f21D(i,j) + f20G(i,j);
			F1(i,j) = f12D2(i,j) + f11D(i,j) + f10G(i,j);
			F0(i,j) = f02D2(i,j) + f01D(i,j) + f00G(i,j);
		}
	}


	D2.Tmmult(Hout,F2,true);
	D.Tmmult(Hout,F1,true);
	G.Tmmult(Hout,F0,true);


	save_matrix(Hout);


}


Tensor<2,6> ElasticProblem::calcDDf_Fourier(const Tensor<1,2> & u, const Tensor<1,2> & up, const Tensor<1,2> & upp, int dir1, int dir2,int index,double fmode){
	Tensor<2,6> ddf;
	Tensor<1,3> ez;
	ez[2] = 1.0;
	Tensor<2,3> Iden3D;
	Iden3D[0][0] = 1.0;
	Iden3D[1][1] = 1.0;
	Iden3D[2][2] = 1.0;
	double R_ref_q = Nodal_Reference_Configuration_Vec[index].get_R();



	double hsc = 12.0*pow(h,2.0)/(1.0 - nu*nu);

	double Rq = R_ref_q;

	double rq = u[0];
	double rp = up[0];
	double rpp = upp[0];
	double zp = up[1];
	double zpp = upp[1];





	// Deformed configuration quantities
	Tensor<1,3> x0_q;
	x0_q[0] = u[0]; x0_q[2] = u[1];
	Tensor<1,3> dx0dS_q;
	dx0dS_q[0] = up[0]; dx0dS_q[2] = up[1];
	Tensor<1,3> dx0dtheta_q;
	dx0dtheta_q[1] = u[0];

	Tensor<1,3> f_q = cross_product(dx0dS_q,dx0dtheta_q);
	Tensor<1,3> n0_q = f_q/f_q.norm();


	Tensor<1,3> ddx0ddS_q;
	ddx0ddS_q[0] = upp[0]; ddx0ddS_q[2] = upp[1];
	Tensor<1,3> ddx0dSdtheta_q;
	ddx0dSdtheta_q[1] = up[0];
	Tensor<1,3> ddx0ddtheta_q;
	ddx0ddtheta_q[0] = -u[0];


	Tensor<2,2> a0c;
	a0c[0][0] = dx0dS_q*dx0dS_q;
	a0c[0][1] = dx0dtheta_q*dx0dS_q;
	a0c[1][0] = a0c[0][1];
	a0c[1][1] = dx0dtheta_q*dx0dtheta_q;

	Tensor<2,2> a0C = inverse2D(a0c);

	Tensor<2,2> Ac;
	Ac[0][0] = 1.0;
	Ac[1][1] = Rq*Rq;

	Tensor<2,2> AC;
	AC[0][0] = 1.0;
	AC[1][1] = 1.0/(Rq*Rq);

	Tensor<2,2> Pt;
	Pt[0][0] = 1.0;
	Pt[1][1] = Rq;

	Tensor<2,2> CovariantMetric;

	CovariantMetric[0][0] = dx0dS_q*dx0dS_q;
	CovariantMetric[0][1] = dx0dS_q*dx0dtheta_q;
	CovariantMetric[1][0] = CovariantMetric[0][1];
	CovariantMetric[1][1] = dx0dtheta_q*dx0dtheta_q;


	Tensor<2,2> ProjectedMetric = 0.5*Pt*(AC*CovariantMetric*AC - AC)*Pt;
	const double stretch_q = sqrt(pow(up[0],2.0) + pow(up[1],2.0));


	Tensor<2,2> Covariant2Form;

	Covariant2Form[0][0] = n0_q*ddx0ddS_q;
	Covariant2Form[0][1] = n0_q*ddx0dSdtheta_q;
	Covariant2Form[1][0] = Covariant2Form[0][1];
	Covariant2Form[1][1] = n0_q*ddx0ddtheta_q;


	Tensor<2,2> Projected2Form;
	Projected2Form = Pt*a0C*Covariant2Form*a0C*Pt;

	Tensor<2,2> Reference2Form = Nodal_Reference_Configuration_Vec[index].get_Covariant_2Form();


	Tensor<2,2> InPlane = ProjectedMetric - inplane_a.get_value(0,0);
	Tensor<2,2> Bending = Projected2Form - Reference2Form - bending_a.get_value(0, 0);


	Material_Vector_InPlane[0].set_Params(Emodv,nu, InPlane);
	Material_Vector_Bending[0].set_Params(Emodv, nu, Bending);

	for (int i = 0; i < 6; i++){
		Tensor<1,3> Di_v, Di_w, Di_vp, Di_wp, Di_vpp, Di_wpp;

		if (dir1 == 0){
			if (i < 3){
				Di_v[i] = 1.;
			} else {
				Di_w[i%3] = 1.;
			}

		} else if (dir1 == 1){
			if (i < 3){
				Di_vp[i] = 1.;
			} else {
				Di_wp[i%3] = 1.;
			}

		} else if (dir1 == 2){
			if (i < 3){
				Di_vpp[i] = 1.;
			} else {
				Di_wpp[i%3] = 1.;
			}
		}

		const Tensor<1,3> dUdS_cos_i_q = Di_vp;
		const Tensor<1,3> dUdS_sin_i_q = Di_wp;
		const Tensor<1,3> dUdtheta_cos_i_q = cross_product(ez,Di_v) + fmode*Di_w;
		const Tensor<1,3> dUdtheta_sin_i_q = cross_product(ez,Di_w) - fmode*Di_v;

		const Tensor<1,3> ddUddS_cos_i_q = Di_vpp;
		const Tensor<1,3> ddUddS_sin_i_q = Di_wpp;
		const Tensor<1,3> ddUdSdtheta_cos_i_q = cross_product(ez,Di_vp) + fmode*Di_wp;
		const Tensor<1,3> ddUdSdtheta_sin_i_q = cross_product(ez,Di_wp) - fmode*Di_vp;
		const Tensor<1,3> ddUddtheta_cos_i_q = cross_product(ez,dUdtheta_cos_i_q) + fmode*cross_product(ez,Di_w) - fmode*fmode*Di_v;
		const Tensor<1,3> ddUddtheta_sin_i_q = cross_product(ez,dUdtheta_sin_i_q) - fmode*cross_product(ez,Di_v) - fmode*fmode*Di_w;

		Tensor<2,2> Da_cos_i_q;
		Da_cos_i_q[0][0] = 2.0 * dx0dS_q*dUdS_cos_i_q;
		Da_cos_i_q[0][1] = (dx0dS_q*dUdtheta_cos_i_q + dx0dtheta_q*dUdS_cos_i_q);
		Da_cos_i_q[1][0] = Da_cos_i_q[0][1];
		Da_cos_i_q[1][1] = 2.0*dx0dtheta_q*dUdtheta_cos_i_q;

		Tensor<2,2> DaC_cos_i_q = -a0C*Da_cos_i_q*a0C;

		Tensor<2,2> Da_sin_i_q;
		Da_sin_i_q[0][0] = 2.0 * dx0dS_q*dUdS_sin_i_q;
		Da_sin_i_q[0][1] = (dx0dS_q*dUdtheta_sin_i_q + dx0dtheta_q*dUdS_sin_i_q);
		Da_sin_i_q[1][0] = Da_sin_i_q[0][1];
		Da_sin_i_q[1][1] = 2.0*dx0dtheta_q*dUdtheta_sin_i_q;

		Tensor<2,2> DaC_sin_i_q = -a0C*Da_sin_i_q*a0C;


		Tensor<2,2> Db1_cos_i_q;
		Db1_cos_i_q[0][0] = n0_q*ddUddS_cos_i_q;
		Db1_cos_i_q[0][1] = n0_q*ddUdSdtheta_cos_i_q;
		Db1_cos_i_q[1][0] = Db1_cos_i_q[0][1];
		Db1_cos_i_q[1][1] = n0_q*ddUddtheta_cos_i_q;

		Tensor<2,2> Db1_sin_i_q;
		Db1_sin_i_q[0][0] = n0_q*ddUddS_sin_i_q;
		Db1_sin_i_q[0][1] = n0_q*ddUdSdtheta_sin_i_q;
		Db1_sin_i_q[1][0] = Db1_sin_i_q[0][1];
		Db1_sin_i_q[1][1] = n0_q*ddUddtheta_sin_i_q;

		Tensor<1,3> df_cos_i_q = (cross_product(dx0dS_q,dUdtheta_cos_i_q) - cross_product(dx0dtheta_q,dUdS_cos_i_q));
		Tensor<2,3> scaledn0proj;
		double fnorm_q = (cross_product(dx0dS_q,dx0dtheta_q).norm());
		scaledn0proj = (Iden3D - outer_product(n0_q,n0_q))/fnorm_q;

		Tensor<1,3> dndxdotu_cos_i_q = scaledn0proj*df_cos_i_q;

		Tensor<2,2> Db2_cos_i_q;
		Db2_cos_i_q[0][0] = dndxdotu_cos_i_q*ddx0ddS_q;
		Db2_cos_i_q[0][1] = dndxdotu_cos_i_q*ddx0dSdtheta_q;
		Db2_cos_i_q[1][0] = Db2_cos_i_q[0][1];
		Db2_cos_i_q[1][1] = dndxdotu_cos_i_q*ddx0ddtheta_q;

		Tensor<1,3> df_sin_i_q = (cross_product(dx0dS_q,dUdtheta_sin_i_q) - cross_product(dx0dtheta_q,dUdS_sin_i_q));
		Tensor<1,3> dndxdotu_sin_i_q = scaledn0proj*df_sin_i_q;

		Tensor<2,2> Db2_sin_i_q;
		Db2_sin_i_q[0][0] = dndxdotu_sin_i_q*ddx0ddS_q;
		Db2_sin_i_q[0][1] = dndxdotu_sin_i_q*ddx0dSdtheta_q;
		Db2_sin_i_q[1][0] = Db2_sin_i_q[0][1];
		Db2_sin_i_q[1][1] = dndxdotu_sin_i_q*ddx0ddtheta_q;

		Tensor<2,2> Db_cos_i_q = Db1_cos_i_q + Db2_cos_i_q;
		Tensor<2,2> Db_sin_i_q = Db1_sin_i_q + Db2_sin_i_q;

		Tensor<2,2> DE_cos_i_q = 0.5*Pt*AC*Da_cos_i_q*AC*Pt;
		Tensor<2,2> DE_sin_i_q = 0.5*Pt*AC*Da_sin_i_q*AC*Pt;

		Tensor<2,2> DB_cos_i_q = Pt*(a0C*Db_cos_i_q*a0C + 2.0*Symmetric(a0C*Covariant2Form*DaC_cos_i_q))*Pt;
		Tensor<2,2> DB_sin_i_q = Pt*(a0C*Db_sin_i_q*a0C + 2.0*Symmetric(a0C*Covariant2Form*DaC_sin_i_q))*Pt;

		for (int j = 0; j < 6; j++){
			Tensor<1,3> Dj_v, Dj_w, Dj_vp, Dj_wp, Dj_vpp, Dj_wpp;

			if (dir2 == 0){
				if (j < 3){
					Dj_v[j] = 1.;
				} else {
					Dj_w[j%3] = 1.;
				}

			} else if (dir2 == 1){
				if (j < 3){
					Dj_vp[j] = 1.;
				} else {
					Dj_wp[j%3] = 1.;
				}

			} else if (dir2 == 2){
				if (j < 3){
					Dj_vpp[j] = 1.;
				} else {
					Dj_wpp[j%3] = 1.;
				}
			}


			const Tensor<1,3> dUdS_cos_j_q = Dj_vp;
			const Tensor<1,3> dUdS_sin_j_q = Dj_wp;
			const Tensor<1,3> dUdtheta_cos_j_q = cross_product(ez,Dj_v) + fmode*Dj_w;
			const Tensor<1,3> dUdtheta_sin_j_q = cross_product(ez,Dj_w) - fmode*Dj_v;

			const Tensor<1,3> ddUddS_cos_j_q = Dj_vpp;
			const Tensor<1,3> ddUddS_sin_j_q = Dj_wpp;
			const Tensor<1,3> ddUdSdtheta_cos_j_q = cross_product(ez,Dj_vp) + fmode*Dj_wp;
			const Tensor<1,3> ddUdSdtheta_sin_j_q = cross_product(ez,Dj_wp) - fmode*Dj_vp;
			const Tensor<1,3> ddUddtheta_cos_j_q = cross_product(ez,dUdtheta_cos_j_q) + fmode*cross_product(ez,Dj_w) - fmode*fmode*Dj_v;
			const Tensor<1,3> ddUddtheta_sin_j_q = cross_product(ez,dUdtheta_sin_j_q) - fmode*cross_product(ez,Dj_v) - fmode*fmode*Dj_w;




			Tensor<2,2> Da_cos_j_q;
			Da_cos_j_q[0][0] = 2.0 * dx0dS_q*dUdS_cos_j_q;
			Da_cos_j_q[0][1] = (dx0dS_q*dUdtheta_cos_j_q + dx0dtheta_q*dUdS_cos_j_q);
			Da_cos_j_q[1][0] = Da_cos_j_q[0][1];
			Da_cos_j_q[1][1] = 2.0*dx0dtheta_q*dUdtheta_cos_j_q;

			Tensor<2,2> DaC_cos_j_q = -a0C*Da_cos_j_q*a0C;

			Tensor<2,2> Da_sin_j_q;
			Da_sin_j_q[0][0] = 2.0 * dx0dS_q*dUdS_sin_j_q;
			Da_sin_j_q[0][1] = (dx0dS_q*dUdtheta_sin_j_q + dx0dtheta_q*dUdS_sin_j_q);
			Da_sin_j_q[1][0] = Da_sin_j_q[0][1];
			Da_sin_j_q[1][1] = 2.0*dx0dtheta_q*dUdtheta_sin_j_q;

			Tensor<2,2> DaC_sin_j_q = -a0C*Da_sin_j_q*a0C;

			Tensor<2,2> DDa_cos_q;
			DDa_cos_q[0][0] = 2.0*dUdS_cos_i_q*dUdS_cos_j_q;
			DDa_cos_q[0][1] = (dUdS_cos_i_q*dUdtheta_cos_j_q + dUdS_cos_j_q*dUdtheta_cos_i_q);
			DDa_cos_q[1][0] = DDa_cos_q[0][1];
			DDa_cos_q[1][1] = 2.0*(dUdtheta_cos_i_q*dUdtheta_cos_j_q);

			Tensor<2,2> DDa_sin_q;
			DDa_sin_q[0][0] = 2.0*dUdS_sin_i_q*dUdS_sin_j_q;
			DDa_sin_q[0][1] = (dUdS_sin_i_q*dUdtheta_sin_j_q + dUdS_sin_j_q*dUdtheta_sin_i_q);
			DDa_sin_q[1][0] = DDa_sin_q[0][1];
			DDa_sin_q[1][1] = 2.0*(dUdtheta_sin_i_q*dUdtheta_sin_j_q);

			Tensor<2,2> Db1_cos_j_q;
			Db1_cos_j_q[0][0] = n0_q*ddUddS_cos_j_q;
			Db1_cos_j_q[0][1] = n0_q*ddUdSdtheta_cos_j_q;
			Db1_cos_j_q[1][0] = Db1_cos_j_q[0][1];
			Db1_cos_j_q[1][1] = n0_q*ddUddtheta_cos_j_q;

			Tensor<2,2> Db1_sin_j_q;
			Db1_sin_j_q[0][0] = n0_q*ddUddS_sin_j_q;
			Db1_sin_j_q[0][1] = n0_q*ddUdSdtheta_sin_j_q;
			Db1_sin_j_q[1][0] = Db1_sin_j_q[0][1];
			Db1_sin_j_q[1][1] = n0_q*ddUddtheta_sin_j_q;

			Tensor<1,3> df_cos_j_q = (cross_product(dx0dS_q,dUdtheta_cos_j_q) - cross_product(dx0dtheta_q,dUdS_cos_j_q));

			Tensor<1,3> dndxdotu_cos_j_q = scaledn0proj*df_cos_j_q;

			Tensor<2,2> Db2_cos_j_q;
			Db2_cos_j_q[0][0] = dndxdotu_cos_j_q*ddx0ddS_q;
			Db2_cos_j_q[0][1] = dndxdotu_cos_j_q*ddx0dSdtheta_q;
			Db2_cos_j_q[1][0] = Db2_cos_j_q[0][1];
			Db2_cos_j_q[1][1] = dndxdotu_cos_j_q*ddx0ddtheta_q;

			Tensor<1,3> df_sin_j_q = (cross_product(dx0dS_q,dUdtheta_sin_j_q) - cross_product(dx0dtheta_q,dUdS_sin_j_q));
			Tensor<1,3> dndxdotu_sin_j_q = scaledn0proj*df_sin_j_q;


			Tensor<2,2> Db2_sin_j_q;
			Db2_sin_j_q[0][0] = dndxdotu_sin_j_q*ddx0ddS_q;
			Db2_sin_j_q[0][1] = dndxdotu_sin_j_q*ddx0dSdtheta_q;
			Db2_sin_j_q[1][0] = Db2_sin_j_q[0][1];
			Db2_sin_j_q[1][1] = dndxdotu_sin_j_q*ddx0ddtheta_q;

			Tensor<2,2> Db_cos_j_q = Db1_cos_j_q + Db2_cos_j_q;
			Tensor<2,2> Db_sin_j_q = Db1_sin_j_q + Db2_sin_j_q;

			Tensor<1,3> DDf_cos_q = cross_product(dUdS_cos_i_q,dUdtheta_cos_j_q) + cross_product(dUdS_cos_j_q,dUdtheta_cos_i_q);
			Tensor<1,3> DDn_cos_q =  -1.0*(n0_q*df_cos_j_q)*dndxdotu_cos_i_q/fnorm_q
					- 1.0*((outer_product(dndxdotu_cos_j_q,n0_q) + outer_product(n0_q,dndxdotu_cos_j_q))*df_cos_i_q)/fnorm_q
					+ scaledn0proj*DDf_cos_q;

			Tensor<1,3> DDf_sin_q = cross_product(dUdS_sin_i_q,dUdtheta_sin_j_q) + cross_product(dUdS_sin_j_q,dUdtheta_sin_i_q);
			Tensor<1,3> DDn_sin_q =  -1.0*(n0_q*df_sin_j_q)*dndxdotu_sin_i_q/fnorm_q
					- 1.0*((outer_product(dndxdotu_sin_j_q,n0_q) + outer_product(n0_q,dndxdotu_sin_j_q))*df_sin_i_q)/fnorm_q
					+ scaledn0proj*DDf_sin_q;

			Tensor<2,2> DDb_cos_q;
			DDb_cos_q[0][0] = DDn_cos_q*ddx0ddS_q + dndxdotu_cos_i_q*ddUddS_cos_j_q + dndxdotu_cos_j_q*ddUddS_cos_i_q;
			DDb_cos_q[0][1] = (DDn_cos_q*ddx0dSdtheta_q + dndxdotu_cos_i_q*ddUdSdtheta_cos_j_q + dndxdotu_cos_j_q*ddUdSdtheta_cos_i_q);
			DDb_cos_q[1][0] = DDb_cos_q[0][1];
			DDb_cos_q[1][1] = (DDn_cos_q*ddx0ddtheta_q + dndxdotu_cos_i_q*ddUddtheta_cos_j_q + dndxdotu_cos_j_q*ddUddtheta_cos_i_q);

			Tensor<2,2> DDb_sin_q;
			DDb_sin_q[0][0] = DDn_sin_q*ddx0ddS_q + dndxdotu_sin_i_q*ddUddS_sin_j_q + dndxdotu_sin_j_q*ddUddS_sin_i_q;
			DDb_sin_q[0][1] = (DDn_sin_q*ddx0dSdtheta_q + dndxdotu_sin_i_q*ddUdSdtheta_sin_j_q + dndxdotu_sin_j_q*ddUdSdtheta_sin_i_q);
			DDb_sin_q[1][0] = DDb_sin_q[0][1];
			DDb_sin_q[1][1] = (DDn_sin_q*ddx0ddtheta_q + dndxdotu_sin_i_q*ddUddtheta_sin_j_q + dndxdotu_sin_j_q*ddUddtheta_sin_i_q);


			//combining terms into the projected forms
			Tensor<2,2> DE_cos_j_q = 0.5*Pt*AC*Da_cos_j_q*AC*Pt;
			Tensor<2,2> DE_sin_j_q = 0.5*Pt*AC*Da_sin_j_q*AC*Pt;

			Tensor<2,2> DDE_cos_q = 0.5*Pt*AC*DDa_cos_q*AC*Pt;
			Tensor<2,2> DDE_sin_q = 0.5*Pt*AC*DDa_sin_q*AC*Pt;

			Tensor<2,2> DDaC_cos_q = -a0C*DDa_cos_q*a0C - 2.0*Symmetric(a0C*Da_cos_i_q*DaC_cos_j_q);
			Tensor<2,2> DDaC_sin_q = -a0C*DDa_sin_q*a0C - 2.0*Symmetric(a0C*Da_sin_i_q*DaC_sin_j_q);

			Tensor<2,2> DB_cos_j_q = Pt*(a0C*Db_cos_j_q*a0C + 2.0*Symmetric(a0C*Covariant2Form*DaC_cos_j_q))*Pt;
			Tensor<2,2> DB_sin_j_q = Pt*(a0C*Db_sin_j_q*a0C + 2.0*Symmetric(a0C*Covariant2Form*DaC_sin_j_q))*Pt;

			Tensor<2,2> DDB_cos_q = Pt*(a0C*DDb_cos_q*a0C +
					2.0*Symmetric(DaC_cos_j_q*Covariant2Form*DaC_cos_i_q + a0C*(Db_cos_i_q*DaC_cos_j_q + Db_cos_j_q*DaC_cos_i_q + Covariant2Form*DDaC_cos_q)))*Pt;


			Tensor<2,2> DDB_sin_q = Pt*(a0C*DDb_sin_q*a0C +
					2.0*Symmetric(DaC_sin_j_q*Covariant2Form*DaC_sin_i_q + a0C*(Db_sin_i_q*DaC_sin_j_q + Db_sin_j_q*DaC_sin_i_q + Covariant2Form*DDaC_sin_q)))*Pt;

			//
			double val = 0.;
			val += pi*Rq*(BilinearProduct(DE_cos_i_q,Material_Vector_InPlane[0].getddQ2ddF(),DE_cos_j_q));
			val += pi*Rq*(BilinearProduct(DE_sin_i_q,Material_Vector_InPlane[0].getddQ2ddF(),DE_sin_j_q));
			val += pi*Rq*(Tensor_Inner(Material_Vector_InPlane[0].getdQ2dF(),DDE_cos_q));
			val += pi*Rq*(Tensor_Inner(Material_Vector_InPlane[0].getdQ2dF(),DDE_sin_q));


			val += pi*Rq*hsc*(BilinearProduct(DB_cos_i_q,Material_Vector_Bending[0].getddQ2ddF(),DB_cos_j_q));
			val += pi*Rq*hsc*(BilinearProduct(DB_sin_i_q,Material_Vector_Bending[0].getddQ2ddF(),DB_sin_j_q));
			val += pi*Rq*hsc*(Tensor_Inner(Material_Vector_Bending[0].getdQ2dF(),DDB_cos_q));
			val += pi*Rq*hsc*(Tensor_Inner(Material_Vector_Bending[0].getdQ2dF(),DDB_sin_q));

			ddf[i][j] = val;
		}

	}


	return ddf;
}


void ElasticProblem::calculate_fd_fourier_stability(double fmode){
	extract_and_sort_values();


	const int Ndof = 6;
	std::vector<double> Si(rvals.size());
	std::vector<Tensor<1,2>> u(rvals.size());
	std::vector<Tensor<1,2>> up(rvals.size());
	std::vector<Tensor<1,2>> upp(rvals.size());

	std::vector<Tensor<2,Ndof>> f00(rvals.size());
	std::vector<Tensor<2,Ndof>> f01(rvals.size());
	std::vector<Tensor<2,Ndof>> f02(rvals.size());
	std::vector<Tensor<2,Ndof>> f11(rvals.size());
	std::vector<Tensor<2,Ndof>> f12(rvals.size());
	std::vector<Tensor<2,Ndof>> f22(rvals.size());


	for (int i = 0; i < rvals.size(); i++){
		Si[i] = rvals[i].first;
		u[i][0] = rvals[i].second;
		u[i][1] = zvals[i].second;
		up[i][0] = xirvals[i].second;
		up[i][1] = xizvals[i].second;


	}



	double dl = Si[1] - Si[0];
	//upp[0] = (up[1] - up[0])/dl;
	//upp[0] = -11.*up[0]/(6. * dl) + 3.*up[1]/dl - 1.5 * up[2]/dl + up[3]/(3.*dl);

	for (int i = 1; i < (rvals.size()-1); i++){
		upp[i] = (up[i+1] - up[i-1])/(2.*dl);

	}




	upp[0] = 2.5*upp[1] - 2.*upp[2] + 0.5*upp[3];
	upp[rvals.size()-1] = 2.5*upp[rvals.size()-2] - 2.*upp[rvals.size()-3] + 0.5*upp[rvals.size()-4];
	/*
	for (int i = 0; i < upp.size(); i++){
		std::cout << upp[i] << std::endl;
	}
	 */

	for (int i = 0; i < rvals.size(); i++){
		f00[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 0, i,fmode);
		f01[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 1, i,fmode);
		f02[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 2, i,fmode);
		f11[i] = calcDDf_Fourier(u[i], up[i], upp[i], 1, 1, i,fmode);
		f12[i] = calcDDf_Fourier(u[i], up[i], upp[i], 1, 2, i,fmode);
		f22[i] = calcDDf_Fourier(u[i], up[i], upp[i], 2, 2, i,fmode);

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;


	}
	Tensor<2,Ndof> Zero;
	int zeropad = 3;
	for (int i = 0; i < zeropad; i++){
		f00[i] = Zero;
		f01[i] = Zero;
		f02[i] = Zero;
		f11[i] = Zero;
		f12[i] = Zero;
		f22[i] = Zero;


		f00[rvals.size() - i - 1] = Zero;
		f01[rvals.size() - i - 1] = Zero;
		f02[rvals.size() - i - 1] = Zero;
		f11[rvals.size() - i - 1] = Zero;
		f12[rvals.size() - i - 1] = Zero;
		f22[rvals.size() - i - 1] = Zero;

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;
	}



	Tensor<2,Ndof> Iden;
	for (int i = 0; i < Ndof; i++){
		Iden[i][i] = 1.;
	}
	// We will also need f11 defined on the half nodes
	std::vector<Tensor<2,Ndof>> f11ext(Ndof*rvals.size()-1);

	for (int i = 0; i < f11ext.size(); i++){
		if (i%2==0){
			f11ext[i] = f11[i/2];
		} else {

			f11ext[i] = calcDDf_Fourier(0.5*(u[i/2]+u[i/2 + 1]),0.5*(up[i/2]+up[i/2 + 1]),0.5*(upp[i/2]+upp[i/2 + 1]),1,1,i/2,fmode);
		}


	}

	for (int i = 0; i < 2*zeropad; i++){
		f11ext[i] = Zero;
		f11ext[f11ext.size() - i - 1] = Zero;
	}




	//I now want to construct the local hessians. I will deal with boundary nodes later
	// the structure will be a vector of tensors at every node, where the vector is the local
	// constructed ones

	std::vector<std::vector<Tensor<2,Ndof>>> DDf22wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> DDf21wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> DDf20w(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df12wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df11wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> Df10w(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f02wpp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f01wp(rvals.size());
	std::vector<std::vector<Tensor<2,Ndof>>> f00w(rvals.size());

	std::vector<std::vector<Tensor<2,Ndof>>> localhessian(rvals.size());

	//initializing the stencil
	int stencilsize = 5;
	for (int i = 0; i < rvals.size(); i++){
		DDf22wpp[i].resize(stencilsize);
		DDf21wp[i].resize(stencilsize);
		DDf20w[i].resize(stencilsize);
		Df12wpp[i].resize(stencilsize);
		Df11wp[i].resize(stencilsize);
		Df10w[i].resize(stencilsize);
		f02wpp[i].resize(stencilsize);
		f01wp[i].resize(stencilsize);
		f00w[i].resize(stencilsize);

		localhessian[i].resize(stencilsize);
	}


	//setting the values for interior nodes
	for (int i = 1; i < (rvals.size() - 1); i++){

		DDf22wpp[i][0] = f22[i-1]/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][2] = (f22[i+1] + 4.*f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22[i+1] + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22[i+1]/pow(dl,4);


		DDf21wp[i][0] = -0.5*(transpose(f12[i-1])/pow(dl,3));
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f12[i-1] - f12[i+1];
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(transpose(f12[i+1]))/pow(dl,3);


		DDf20w[i][1] = transpose(f02[i-1])/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = transpose(f02[i+1])/pow(dl,2);


		Df12wpp[i][0] = -0.5*f12[i-1]/pow(dl,3);
		Df12wpp[i][1] = f12[i-1]/pow(dl,3);
		Df12wpp[i][2] = 0.5*(f12[i+1]-f12[i-1])/pow(dl,3);
		Df12wpp[i][3] = -f12[i+1]/pow(dl,3);
		Df12wpp[i][4] = 0.5*f12[i+1]/pow(dl,3);


		Df11wp[i][1] = f11ext[2*i-1]/pow(dl,2);
		Df11wp[i][2] = -(f11ext[2*i-1] + f11ext[2*i+1])/pow(dl,2);
		Df11wp[i][3] = f11ext[2*i+1]/pow(dl,2);




		Df10w[i][1] = -0.5*transpose(f01[i-1])/dl;
		Df10w[i][3] = 0.5*transpose(f01[i+1])/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];

	}


	//Dealing with boundary conditions. for now I will use both bcs from integration by parts
	//beginning with 0 boundary

	Tensor<2,Ndof> f20m1 = 0.*(2.*transpose(f02[0]) - transpose(f02[1]));
	Tensor<2,Ndof> f21m1 = 0.*(2.*transpose(f12[0]) - transpose(f12[1]));
	Tensor<2,Ndof> f11m1 = 0.*(2.*f11[0] - f11[1]);
	Tensor<2,Ndof> f21m1_2 = 0.*(0.5*(f21m1 + transpose(f12[0])));
	Tensor<2,Ndof> f21p1_2 = 0.*calcDDf_Fourier(0.5*(u[0]+u[1]),0.5*(up[0]+up[1]),0.5*(upp[0]+upp[1]),2,1,1,fmode);
	Tensor<2,Ndof> f22m1 = 0.*(2.*f22[0] - f22[1]);
	Tensor<2,Ndof> f01m1 = 0.*(2.*f01[0] - f01[1]);



	std::vector<Tensor<2,Ndof>> Df21wp0(3);
	Df21wp0[0] = f21m1_2/pow(dl,2);
	Df21wp0[1] = -1.0*(f21m1_2 + f21p1_2)/pow(dl,2);
	Df21wp0[2] = f21p1_2/pow(dl,2);

	std::vector<Tensor<2,Ndof>> Df20w0(3);
	Df20w0[0] = -.5*f20m1/dl;
	Df20w0[2] = .5*transpose(f02[1])/dl;

	std::vector<Tensor<2,Ndof>> Df22wpp0(5);
	Df22wpp0[0] = -0.5*f22m1/pow(dl,3);
	Df22wpp0[1] = f22m1/pow(dl,3);
	Df22wpp0[2] = 0.5*(f22[1] - f22m1)/pow(dl,3);
	Df22wpp0[3] = -f22[1]/pow(dl,3);
	Df22wpp0[4] = 0.5*f22[1]/pow(dl,3);






	{
		int i = 0;

		DDf22wpp[i][0] = f22m1/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22m1)/pow(dl,4);
		DDf22wpp[i][2] = (f22[i+1] + 4.*f22[i] + f22m1)/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22[i+1] + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22[i+1]/pow(dl,4);


		DDf21wp[i][0] = -0.5*f21m1/pow(dl,3);
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f21m1 - f12[i+1];
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(transpose(f12[i+1]))/pow(dl,3);


		DDf20w[i][1] = f20m1/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = transpose(f02[i+1])/pow(dl,2);


		Df12wpp[i][0] = -0.5*f21m1/pow(dl,3);
		Df12wpp[i][1] = f21m1/pow(dl,3);
		Df12wpp[i][2] = 0.5*(f12[i+1]-f21m1)/pow(dl,3);
		Df12wpp[i][3] = -f12[i+1]/pow(dl,3);
		Df12wpp[i][4] = 0.5*f12[i+1]/pow(dl,3);


		Df11wp[i][1] = f11m1/pow(dl,2);
		Df11wp[i][2] = -(f11m1 + f11ext[2*i+1])/pow(dl,2);
		Df11wp[i][3] = f11ext[2*i+1]/pow(dl,2);


		Df10w[i][1] = -0.5*transpose(f01m1)/dl;
		Df10w[i][3] = 0.5*transpose(f01[i+1])/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];
	}



	//Going to the L boundary boundary
	int Nend = rvals.size();


	Tensor<2,Ndof> f20Np1 = 0.*(2.*transpose(f02[Nend-1]) - transpose(f02[Nend-2]));
	Tensor<2,Ndof> f21Np1 = 0.*(2.*transpose(f12[Nend-1]) - transpose(f12[Nend-2]));
	Tensor<2,Ndof> f11Np1 = 0.*(2.*f11[Nend-1] - f11[Nend-2]);
	Tensor<2,Ndof> f21Np1_2 = 0.*(0.5*(f21Np1 + transpose(f12[Nend-1])));
	Tensor<2,Ndof> f21Nm1_2 = 0.*calcDDf_Fourier(0.5*(u[Nend-2]+u[Nend-1]),0.5*(up[Nend-2]+up[Nend-1]),0.5*(upp[Nend-2]+upp[Nend-1]),2,1,Nend-1,fmode);
	Tensor<2,Ndof> f22Np1 = 0.*(2.*f22[Nend-1] - f22[Nend-2]);
	Tensor<2,Ndof> f01Np1 = 0.*(2.*f01[Nend-1] - f01[Nend-2]);



	//finding el equations for the N boundaries
	{
		int i = Nend-1;

		DDf22wpp[i][0] = f22[i-1]/pow(dl,4);
		DDf22wpp[i][1] = -2.*(f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][2] = (f22Np1 + 4.*f22[i] + f22[i-1])/pow(dl,4);
		DDf22wpp[i][3] = -2.*(f22Np1 + f22[i])/pow(dl,4);
		DDf22wpp[i][4] = f22Np1/pow(dl,4);


		//this one has a problem
		DDf21wp[i][0] = -0.5*(transpose(f12[i-1])/pow(dl,3));
		DDf21wp[i][1] = transpose(f12[i])/pow(dl,3);
		Tensor<2,Ndof> f12difftemp = f12[i-1] - transpose(f21Np1);
		DDf21wp[i][2] = 0.5*(transpose(f12difftemp))/pow(dl,3);
		DDf21wp[i][3] = -transpose(f12[i])/pow(dl,3);
		DDf21wp[i][4] = 0.5*(f21Np1)/pow(dl,3);



		DDf20w[i][1] = transpose(f02[i-1])/pow(dl,2);
		DDf20w[i][2] = -2.*transpose(f02[i])/pow(dl,2);
		DDf20w[i][3] = f20Np1/pow(dl,2);


		Df12wpp[i][0] = -0.5*f12[i-1]/pow(dl,3);
		Df12wpp[i][1] = f12[i-1]/pow(dl,3);
		Df12wpp[i][2] = 0.5*(transpose(f21Np1)-f12[i-1])/pow(dl,3);
		Df12wpp[i][3] = -transpose(f21Np1)/pow(dl,3);
		Df12wpp[i][4] = 0.5*transpose(f21Np1)/pow(dl,3);

		Df11wp[i][1] = f11ext[2*i-1]/pow(dl,2);
		Tensor<2,Ndof> f11Np1_2 = 0.0*0.5*(f11[i] + f11Np1);
		Df11wp[i][2] = -(f11ext[2*i-1] + f11Np1_2)/pow(dl,2);
		Df11wp[i][3] = f11Np1_2/pow(dl,2);


		Df10w[i][1] = -0.5*transpose(f01[i-1])/dl;
		Df10w[i][3] = 0.5*transpose(f01Np1)/dl;


		f02wpp[i][1] = f02[i]/pow(dl,2);
		f02wpp[i][2] = -2.*f02[i]/pow(dl,2);
		f02wpp[i][3] = f02[i]/pow(dl,2);


		f01wp[i][1] = -0.5*f01[i]/dl;
		f01wp[i][3] = 0.5*f01[i]/dl;


		f00w[i][2] = f00[i];


	}




	for (int i = 0; i < localhessian.size(); i++){
		for (int j = 0; j < stencilsize; j++){
			localhessian[i][j] = DDf22wpp[i][j] + DDf21wp[i][j] + DDf20w[i][j] - Df12wpp[i][j] - Df11wp[i][j] - Df10w[i][j] + f02wpp[i][j] + f01wp[i][j] + f00w[i][j];

		}
	}




	LAPACKFullMatrix<double> FDHessian;


	FDHessian.reinit(Ndof*(rvals.size()+4),Ndof*(rvals.size()+4));
	FDHessian = 0.;



	for (int i = 0; i < (rvals.size()); i++){
		for (int j = 0; j < stencilsize; j++){
			for (int k = 0; k < Ndof; k++){
				for (int l = 0; l < Ndof; l++){
					FDHessian(Ndof*(i+2)+k,Ndof*(i+j)+l) = localhessian[i][j][k][l];
				}
			}
		}
	}

	save_matrix(FDHessian);
}
// @sect4{Step4::solve}


void ElasticProblem::calculate_fd_fourier_stability2(double fmode){


	extract_and_sort_values();


	const int Ndof = 6;
	std::vector<double> Si(rvals.size());
	std::vector<Tensor<1,2>> u(rvals.size());
	std::vector<Tensor<1,2>> up(rvals.size());
	std::vector<Tensor<1,2>> upp(rvals.size());

	std::vector<Tensor<2,Ndof>> f00(rvals.size());
	std::vector<Tensor<2,Ndof>> f01(rvals.size());
	std::vector<Tensor<2,Ndof>> f02(rvals.size());
	std::vector<Tensor<2,Ndof>> f11(rvals.size());
	std::vector<Tensor<2,Ndof>> f12(rvals.size());
	std::vector<Tensor<2,Ndof>> f22(rvals.size());


	for (int i = 0; i < rvals.size(); i++){
		Si[i] = rvals[i].first;
		u[i][0] = rvals[i].second;
		u[i][1] = zvals[i].second;
		up[i][0] = xirvals[i].second;
		up[i][1] = xizvals[i].second;


	}

	std::vector<Tensor<1,2>> utilde(rvals.size()-1);
	std::vector<Tensor<1,2>> utilde_p(rvals.size()-1);
	std::vector<Tensor<1,2>> utilde_pp(rvals.size()-1);

	std::vector<Tensor<2,Ndof>> f00tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f01tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f02tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f11tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f12tilde(rvals.size()-1);
	std::vector<Tensor<2,Ndof>> f22tilde(rvals.size()-1);

	double dl = Si[1] - Si[0];
	//upp[0] = (up[1] - up[0])/dl;
	//upp[0] = -11.*up[0]/(6. * dl) + 3.*up[1]/dl - 1.5 * up[2]/dl + up[3]/(3.*dl);

	for (int i = 1; i < (rvals.size()-1); i++){
		upp[i] = (up[i+1] - up[i-1])/(2.*dl);

	}




	upp[0] = 2.5*upp[1] - 2.*upp[2] + 0.5*upp[3];
	upp[rvals.size()-1] = 2.5*upp[rvals.size()-2] - 2.*upp[rvals.size()-3] + 0.5*upp[rvals.size()-4];
	/*
		for (int i = 0; i < upp.size(); i++){
			std::cout << upp[i] << std::endl;
		}
	 */

	for (int i = 0; i < rvals.size(); i++){
		f00[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 0, i,fmode);
		f01[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 1, i,fmode);
		f02[i] = calcDDf_Fourier(u[i], up[i], upp[i], 0, 2, i,fmode);
		f11[i] = calcDDf_Fourier(u[i], up[i], upp[i], 1, 1, i,fmode);
		f12[i] = calcDDf_Fourier(u[i], up[i], upp[i], 1, 2, i,fmode);
		f22[i] = calcDDf_Fourier(u[i], up[i], upp[i], 2, 2, i,fmode);

		//std::cout << f00[i][0][0] << "," << f00[i][0][1] << "," << f00[i][1][1] << "," << std::endl;


	}

	//setting half values

	for (int i = 0; i < utilde.size(); i++){
		utilde[i] = 0.5*(u[i] + u[i+1]);
		utilde_p[i] = 0.5*(up[i] + up[i+1]);
		utilde_pp[i] = 0.5*(upp[i] + upp[i+1]);

		f00tilde[i] = 0.5*(f00[i]+f00[i+1]);
		f01tilde[i] = 0.5*(f01[i]+f01[i+1]);
		f02tilde[i] = 0.5*(f02[i]+f02[i+1]);
		f11tilde[i] = 0.5*(f11[i]+f11[i+1]);
		f12tilde[i] = 0.5*(f12[i]+f12[i+1]);
		f22tilde[i] = 0.5*(f22[i]+f22[i+1]);

	}




	Tensor<2,Ndof> Zero;


	Tensor<2,Ndof> Iden;
	for (int i = 0; i < Ndof; i++){
		Iden[i][i] = 1.;
	}
	// We will also need f11 defined on the half nodes


	int Nend = rvals.size();


	std::vector<std::vector<Tensor<2,Ndof>>> Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> D2i(rvals.size()-1);

	int StencilSize2 = 4;
	for (int i = 0; i < (rvals.size()-1); i++){
		Gi[i].resize(StencilSize2);
		Di[i].resize(StencilSize2);
		D2i[i].resize(StencilSize2);
	}

	for (int i = 0; i < Gi.size(); i++){
		Gi[i][0] = Zero;
		Gi[i][1] = 1.0*Iden;
		Gi[i][2] = 0.0*Iden;
		Gi[i][3] = Zero;

		Di[i][0] = Zero;
		Di[i][1] = -Iden/dl;
		Di[i][2] = Iden/dl;
		Di[i][3] = Zero;

		if (i == 0) {
			D2i[i][0] = Zero;
			D2i[i][1] = Iden/pow(dl,2);
			D2i[i][2] = -2.*Iden/pow(dl,2);
			D2i[i][3] = Iden/pow(dl,2);

		} else if (i == (Gi.size()-1)) {

			D2i[i][0] = Iden/pow(dl,2);
			D2i[i][1] = -2.*Iden/pow(dl,2);
			D2i[i][2] = Iden/pow(dl,2);
			D2i[i][3] = Zero;

		} else {
			/*
				D2i[i][0] = Iden/(2.*pow(dl,2));
				D2i[i][1] = -Iden/(2.*pow(dl,2));
				D2i[i][2] = -Iden/(2.*pow(dl,2));
				D2i[i][3] = Iden/(2.*pow(dl,2));
			 */
			D2i[i][0] = Zero;
			D2i[i][1] = Iden/(pow(dl,2));
			D2i[i][2] = -2.*Iden/pow(dl,2);
			D2i[i][3] = Iden/pow(dl,2);


		}
	}


	std::vector<std::vector<Tensor<2,Ndof>>> f22D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f21Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f20Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f12D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f11Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f10Gi(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f02D2i(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f01Di(rvals.size()-1);
	std::vector<std::vector<Tensor<2,Ndof>>> f00Gi(rvals.size()-1);

	for (int i = 0; i < (rvals.size()-1); i++){
		f22D2i[i].resize(StencilSize2);
		f21Di[i].resize(StencilSize2);
		f20Gi[i].resize(StencilSize2);
		f12D2i[i].resize(StencilSize2);
		f11Di[i].resize(StencilSize2);
		f10Gi[i].resize(StencilSize2);
		f02D2i[i].resize(StencilSize2);
		f01Di[i].resize(StencilSize2);
		f00Gi[i].resize(StencilSize2);

		for (int j = 0; j < StencilSize2; j++){
			f22D2i[i][j] = f22tilde[i]*D2i[i][j];
			f21Di[i][j] = transpose(f12tilde[i])*Di[i][j];
			f20Gi[i][j] = transpose(f02tilde[i])*Gi[i][j];

			f12D2i[i][j] = f12tilde[i]*D2i[i][j];
			f11Di[i][j] = f11tilde[i]*Di[i][j];
			f10Gi[i][j] = transpose(f01tilde[i])*Gi[i][j];

			f02D2i[i][j]= f02tilde[i]*D2i[i][j];
			f01Di[i][j] = f01tilde[i]*Di[i][j];
			f00Gi[i][j] = f00tilde[i]*Gi[i][j];

		}
	}


	LAPACKFullMatrix<double> G; G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); G = 0.;
	LAPACKFullMatrix<double> D; D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); D = 0.;
	LAPACKFullMatrix<double> D2; D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); D2 = 0.;
	LAPACKFullMatrix<double> f22D2; f22D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f22D2 = 0.;
	LAPACKFullMatrix<double> f21D; f21D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f21D = 0.;
	LAPACKFullMatrix<double> f20G; f20G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f20G = 0.;
	LAPACKFullMatrix<double> f12D2; f12D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f12D2 = 0.;
	LAPACKFullMatrix<double> f11D; f11D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f11D = 0.;
	LAPACKFullMatrix<double> f10G; f10G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f10G = 0.;
	LAPACKFullMatrix<double> f02D2; f02D2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f02D2 = 0.;
	LAPACKFullMatrix<double> f01D; f01D.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f01D = 0.;
	LAPACKFullMatrix<double> f00G; f00G.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); f00G = 0.;


	for (int i = 0; i < Di.size(); i++){
		if (i==0){
			for (int j = 1; j < StencilSize2; j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];

					}
				}
			}

		} else if (i == (Di.size()-1)){
			for (int j = 0; j < StencilSize2-1; j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];
					}
				}
			}
		} else {
			for (int j = 0; j < (StencilSize2); j++){
				for (int k = 0; k < Ndof; k++){
					for (int l = 0; l < Ndof; l++){
						G(Ndof * i + k, Ndof * (i + j - 1) + l) += Gi[i][j][k][l];
						D(Ndof * i + k, Ndof * (i + j - 1) + l) += Di[i][j][k][l];
						D2(Ndof * i + k, Ndof * (i + j - 1) + l) += D2i[i][j][k][l];

						f22D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f22D2i[i][j][k][l];
						f21D(Ndof * i + k, Ndof * (i + j - 1) + l) += f21Di[i][j][k][l];
						f20G(Ndof * i + k, Ndof * (i + j - 1) + l) += f20Gi[i][j][k][l];

						f12D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f12D2i[i][j][k][l];
						f11D(Ndof * i + k, Ndof * (i + j - 1) + l) += f11Di[i][j][k][l];
						f10G(Ndof * i + k, Ndof * (i + j - 1) + l) += f10Gi[i][j][k][l];

						f02D2(Ndof * i + k, Ndof * (i + j - 1) + l) += f02D2i[i][j][k][l];
						f01D(Ndof * i + k, Ndof * (i + j - 1) + l) += f01Di[i][j][k][l];
						f00G(Ndof * i + k, Ndof * (i + j - 1) + l) += f00Gi[i][j][k][l];
					}
				}
			}

		}
	}

	LAPACKFullMatrix<double> Hout; Hout.reinit(Ndof*rvals.size(),Ndof*rvals.size()); Hout = 0.;
	LAPACKFullMatrix<double> F2; F2.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F2 = 0.;
	LAPACKFullMatrix<double> F1; F1.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F1 = 0.;
	LAPACKFullMatrix<double> F0; F0.reinit(Ndof*(rvals.size()-1),Ndof*rvals.size()); F0 = 0.;

	for (int i = 0; i < f22D2.m(); i++){
		for (int j = 0; j < f22D2.n(); j++){
			F2(i,j) = f22D2(i,j) + f21D(i,j) + f20G(i,j);
			F1(i,j) = f12D2(i,j) + f11D(i,j) + f10G(i,j);
			F0(i,j) = f02D2(i,j) + f01D(i,j) + f00G(i,j);
		}
	}


	D2.Tmmult(Hout,F2,true);
	D.Tmmult(Hout,F1,true);
	G.Tmmult(Hout,F0,true);


	save_matrix(Hout);

}

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

void ElasticProblem::extract_and_sort_values(){
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

	rvals.clear();
	zvals.clear();
	xirvals.clear();
	xizvals.clear();

	for (unsigned int i = 0; i < ndofs; i++) {

		if (is_r_comp[i]){
			rvals.push_back(std::make_pair(support_points[i][0],solution[i]));
		} else if (is_z_comp[i]) {
			zvals.push_back(std::make_pair(support_points[i][0],solution[i]));
		} else if (is_xi_r_comp[i]) {
			xirvals.push_back(std::make_pair(support_points[i][0],solution[i]));
		} else if (is_xi_z_comp[i]){
			xizvals.push_back(std::make_pair(support_points[i][0],solution[i]));
		}
	}

	sort(rvals.begin(),rvals.end());
	sort(zvals.begin(),zvals.end());
	sort(xirvals.begin(),xirvals.end());
	sort(xizvals.begin(),xizvals.end());

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


Tensor<2,2> ElasticProblem::transpose2D(Tensor<2,2> & min){
	Tensor<2,2> mout;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			mout[i][j] = min[j][i];
		}
	}
	return mout;
}

Tensor<2,6> ElasticProblem::transpose6D(Tensor<2,6> & min){
	Tensor<2,6> mout;
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			mout[i][j] = min[j][i];
		}
	}
	return mout;
}
Tensor<2,2> ElasticProblem::inverse2D(Tensor<2,2> & min){
	Tensor<2,2> mout;
	double detmin = min[0][0]*min[1][1] - min[0][1]*min[1][0];
	mout[0][0] = min[1][1]/detmin;
	mout[0][1] = -1.0*min[0][1]/detmin;
	mout[1][0] = -1.0*min[1][0]/detmin;
	mout[1][1] = min[0][0]/detmin;
	return mout;
}
Tensor<2,6> ElasticProblem::inverse6D(Tensor<2,6> & min){
	Tensor<2,6> mout;
	std::cout << "Trying the 6x6 matrix invert" << std::endl;
	LAPACKFullMatrix<double> mla;
	mla.reinit(6, 6);
	mla = 0.;

	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			mla(i,j) = min[i][j];
			std::cout << min[i][j] << ",";
		}
		std::cout << std::endl;
	}

	mla.invert();
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			mout[i][j] = mla(i,j);
		}
	}

	std::cout << "Did it" << std::endl;
	return mout;

}
template<int d>
Tensor<2,d> ElasticProblem::Symmetric(const Tensor<2,d> & A){
	Tensor<2,d> Aout;
	for (unsigned int i = 0; i < d; i++){
		for (unsigned int j = 0; j < d; j++){
			Aout[i][j] = 0.5*(A[i][j] + A[j][i]);
		}
	}
	return Aout;
}

Tensor<1,3> ElasticProblem::cross_product(const Tensor<1,3> & v1, const Tensor<1,3> & v2){
	Tensor<1,3> vout;
	vout[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vout[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vout[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return vout;
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
void ElasticProblem::save_matrix(const LAPACKFullMatrix<double> & min){
	std::ofstream rfile;
	std::string rfilename = "stabmatrices/stabmatrix" + std::to_string(solveiteration) + ".csv";
	rfile.open(rfilename);

	for (unsigned int i = 0; i < min.n_rows(); ++i){
		for (unsigned int j = 0; j < min.n_cols(); j++){
			rfile << std::setprecision(10) << min(i,j) << ',';
		}
		rfile << '\n';

	}


	rfile.close();
}
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
	std::cout << "Set refinement level: " ;
	std::cin >> refinelevel;
	std::cout << std::endl;

	std::cout << "Choose Fourier Mode (int): ";
	std::cin >> fmodein;
	std::cout << std::endl;
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

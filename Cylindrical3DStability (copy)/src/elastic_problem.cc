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
: fe(FESystem<DIM>(FE_Q<DIM>(2), 3),4),
  x0fe(FE_Q<DIM>(2), 6),
  x0dof_handler(x0triangulation)
, dof_handler(triangulation){}


void ElasticProblem::solve_path(){


	h = dat.h;
	int Ntsteps = dat.Nsteps;
	nu = dat.nu;
	double bendingdefmag = dat.bendingdefmag;
	std::cout << bendingdefmag << std::endl;
	double inplanedefmag = dat.inplanedefmag;
	std::cout << inplanedefmag << std::endl;
	std::cout << "--------------------------------" << std::endl;
	int Nmax = 5;
	Tensor<2,2> bending_dir;
	bending_dir[0][0] = dat.bending00;
	bending_dir[0][1] = dat.bending01;
	bending_dir[1][0] = dat.bending01;
	bending_dir[1][1] = dat.bending11;
	Tensor<2,2> inplane_dir;
	inplane_dir[0][0] = dat.inplane00;
	inplane_dir[0][1] = dat.inplane01;
	inplane_dir[1][0] = dat.inplane01;
	inplane_dir[1][1] = dat.inplane11;

	std::vector<double> nvec = linspace(1.0,(double)Nmax,Nmax);
	std::vector<double> bendingdefmagvec = linspace(0.0,bendingdefmag,Ntsteps);
	std::vector<double> inplanedefmagvec = linspace(0.0,inplanedefmag,Ntsteps);

	initialize_reference_config();
	assemble_constraint_system();



	Vector<double> ev;
	Evals.resize(Nmax, Ntsteps); Evals.setZero();

	std::vector<std::vector<Vector<double>>> Evalsall;
	Evalsall.resize(Ntsteps);
	for (unsigned int i = 0; i < Ntsteps; i++){
		Evalsall[i].resize(Nmax);
		for (unsigned int j = 0; j < Nmax; j++){
			Evalsall[i][j].reinit(10);
		}

	}

	for (unsigned int ts = 0; ts < Ntsteps; ts++){

		defstep = bendingdefmagvec[ts];
		Evalsall[ts].resize(Nmax);
		std::cout << "TimeStep: " << ts << std::endl;
		timestep = ts;
		load_state(timestep);

		//current_bending_a = defstep * bending_dir;
		current_bending_a = bendingdefmagvec[ts]*bending_dir;
		current_inplane_a = inplanedefmagvec[ts]*inplane_dir;
		for (unsigned int i = 0; i < Nmax; i++){
			nval = nvec[i];
			update_internal_metrics();
			std::cout << "-----------------------" << std::endl;
			assemble_system();

			construct_reduced_mappings();
			Evalsall[ts][i] = calculate_stability(ts);
		}
		save_stability(Evalsall);
	}

}

void ElasticProblem::make_grid()
{


	GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
	triangulation.refine_global(dat.refinelevel);
	//triangulation.refine_global(7);
	std::cout << "   Number of active cells: " << triangulation.n_active_cells()
																																																																																					<< std::endl << "   Total number of cells: "
																																																																																					<< triangulation.n_cells() << std::endl;
}



void ElasticProblem::setup_system()
{


	dof_handler.distribute_dofs(fe);
	solution.reinit(dof_handler.n_dofs());
	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            																																																																																																												<< std::endl;


	DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler,
			dsp,
			constraints,
			/*keep_constrained_dofs = */ false);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);
	constraint_matrix.reinit(sparsity_pattern);

	Material_Vector_InPlane.resize(triangulation.n_active_cells());
	Material_Vector_Bending.resize(triangulation.n_active_cells());
	Reference_Configuration_Vec.resize(triangulation.n_active_cells());

	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);
	inplane_a.resize(triangulation.n_active_cells(), quadrature_formula.size());
	bending_a.resize(triangulation.n_active_cells(), quadrature_formula.size());



	std::ofstream out("sparsity_pattern.svg");
	sparsity_pattern.print_svg(out);




}

void ElasticProblem::setup_constraints(){
	constraints.clear();

	const int ndofs = dof_handler.n_dofs();

	//Creating masks for degrees of freedom
	////////////////////
	std::vector<bool> vr_components = create_bool_vector(12,0);
	ComponentMask vr_mask(vr_components);

	std::vector<bool> is_vr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vr_mask, is_vr_comp);
	////////////////////
	std::vector<bool> vtheta_components = create_bool_vector(12,1);
	ComponentMask vtheta_mask(vtheta_components);

	std::vector<bool> is_vtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vtheta_mask, is_vtheta_comp);
	/////////////////////////
	std::vector<bool> vz_components = create_bool_vector(12,2);
	ComponentMask vz_mask(vz_components);

	std::vector<bool> is_vz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vz_mask, is_vz_comp);

	////////////////
	std::vector<bool> wr_components = create_bool_vector(12,3);
	ComponentMask wr_mask(wr_components);

	std::vector<bool> is_wr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wr_mask, is_wr_comp);
	///////////////////
	std::vector<bool> wtheta_components = create_bool_vector(12,4);
	ComponentMask wtheta_mask(wtheta_components);

	std::vector<bool> is_wtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wtheta_mask, is_wtheta_comp);
	//////////////////
	std::vector<bool> wz_components = create_bool_vector(12,5);
	ComponentMask wz_mask(wz_components);

	std::vector<bool> is_wz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wz_mask, is_wz_comp);
	//////////////////////
	std::vector<bool> xir_components = create_bool_vector(12,6);
	ComponentMask xir_mask(xir_components);

	std::vector<bool> is_xir_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xir_mask, is_xir_comp);
	/////////////////////
	std::vector<bool> xitheta_components = create_bool_vector(12,7);
	ComponentMask xitheta_mask(xitheta_components);

	std::vector<bool> is_xitheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xitheta_mask, is_xitheta_comp);
	//////////////////////
	std::vector<bool> xiz_components = create_bool_vector(12,8);
	ComponentMask xiz_mask(xiz_components);

	std::vector<bool> is_xiz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xiz_mask, is_xiz_comp);
	//////////////////////
	std::vector<bool> etar_components = create_bool_vector(12,9);
	ComponentMask etar_mask(etar_components);

	std::vector<bool> is_etar_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etar_mask, is_etar_comp);
	//////////////////////
	std::vector<bool> etatheta_components = create_bool_vector(12,10);
	ComponentMask etatheta_mask(etatheta_components);

	std::vector<bool> is_etatheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etatheta_mask, is_etatheta_comp);
	//////////////////////
	std::vector<bool> etaz_components = create_bool_vector(12,11);
	ComponentMask etaz_mask(etaz_components);

	std::vector<bool> is_etaz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etaz_mask, is_etaz_comp);



	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	for (unsigned int i = 0; i < ndofs; i++){

		if (fabs(support_points[i][0]) < pow(10.0,-7.0) && is_vtheta_comp[i]){
			//constraints.add_line(i);
		}

	}

	constraints.close();

}



void ElasticProblem::assemble_system()
{
	system_matrix = 0;

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
	FEValues<DIM> fe_values_x0(x0fe,quadrature_formula,
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


	std::vector<double> vr_n_q(n_q_points);
	std::vector<double> vtheta_n_q(n_q_points);
	std::vector<double> vz_n_q(n_q_points);

	std::vector<double> wr_n_q(n_q_points);
	std::vector<double> wtheta_n_q(n_q_points);
	std::vector<double> wz_n_q(n_q_points);

	std::vector<Tensor<1,1>> dvr_n_q(n_q_points);
	std::vector<Tensor<1,1>> dvtheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> dvz_n_q(n_q_points);

	std::vector<Tensor<1,1>> dwr_n_q(n_q_points);
	std::vector<Tensor<1,1>> dwtheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> dwz_n_q(n_q_points);

	std::vector<double> xir_n_q(n_q_points);
	std::vector<double> xitheta_n_q(n_q_points);
	std::vector<double> xiz_n_q(n_q_points);

	std::vector<double> etar_n_q(n_q_points);
	std::vector<double> etatheta_n_q(n_q_points);
	std::vector<double> etaz_n_q(n_q_points);

	std::vector<Tensor<1,1>> dxir_n_q(n_q_points);
	std::vector<Tensor<1,1>> dxitheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> dxiz_n_q(n_q_points);

	std::vector<Tensor<1,1>> detar_n_q(n_q_points);
	std::vector<Tensor<1,1>> detatheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> detaz_n_q(n_q_points);

	std::vector<double> r_q(n_q_points);
	std::vector<double> z_q(n_q_points);
	std::vector<double> R_q(n_q_points);
	std::vector<double> Phi_q(n_q_points);
	std::vector<double> dPhi_q(n_q_points);
	std::vector<double> dr_q(n_q_points);
	std::vector<double> dz_q(n_q_points);
	std::vector<Tensor<1,1>> ddr_q(n_q_points);
	std::vector<Tensor<1,1>> ddz_q(n_q_points);

	Tensor<1,3> ez;
	ez[2] = 1.0;
	Tensor<2,3> Iden3D;
	Iden3D[0][0] = 1.0;
	Iden3D[1][1] = 1.0;
	Iden3D[2][2] = 1.0;


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
		//const auto * x0cell = x0dof_handler.active_cell_iterators();
		//fe_values_x0.reinit(x0cell);
		cell_matrix = 0;
		cell_rhs    = 0;





		//
		const FEValuesExtractors::Scalar vr(0);
		const FEValuesExtractors::Scalar vtheta(1);
		const FEValuesExtractors::Scalar vz(2);
		const FEValuesExtractors::Scalar wr(3);
		const FEValuesExtractors::Scalar wtheta(4);
		const FEValuesExtractors::Scalar wz(5);

		const FEValuesExtractors::Scalar xir(6);
		const FEValuesExtractors::Scalar xitheta(7);
		const FEValuesExtractors::Scalar xiz(8);
		const FEValuesExtractors::Scalar etar(9);
		const FEValuesExtractors::Scalar etatheta(10);
		const FEValuesExtractors::Scalar etaz(11);

		const FEValuesExtractors::Scalar r(0);
		const FEValuesExtractors::Scalar z(1);
		const FEValuesExtractors::Scalar lr(2);
		const FEValuesExtractors::Scalar dr(4);
		const FEValuesExtractors::Scalar dz(5);
		//

		//
		fe_values[vr].get_function_values(solution,vr_n_q);
		fe_values[vtheta].get_function_values(solution,vtheta_n_q);
		fe_values[vz].get_function_values(solution,vz_n_q);

		fe_values[wr].get_function_values(solution,wr_n_q);
		fe_values[wtheta].get_function_values(solution,wtheta_n_q);
		fe_values[wz].get_function_values(solution,wz_n_q);

		fe_values[xir].get_function_values(solution,xir_n_q);
		fe_values[xitheta].get_function_values(solution,xitheta_n_q);
		fe_values[xiz].get_function_values(solution,xiz_n_q);

		fe_values[etar].get_function_values(solution,etar_n_q);
		fe_values[etatheta].get_function_values(solution,etatheta_n_q);
		fe_values[etaz].get_function_values(solution,etaz_n_q);


		fe_values[vr].get_function_gradients(solution,dvr_n_q);
		fe_values[vtheta].get_function_gradients(solution,dvtheta_n_q);
		fe_values[vz].get_function_gradients(solution,dvz_n_q);

		fe_values[wr].get_function_gradients(solution,dwr_n_q);
		fe_values[wtheta].get_function_gradients(solution,dwtheta_n_q);
		fe_values[wz].get_function_gradients(solution,dwz_n_q);

		fe_values[xir].get_function_gradients(solution,dxir_n_q);
		fe_values[xitheta].get_function_gradients(solution,dxitheta_n_q);
		fe_values[xiz].get_function_gradients(solution,dxiz_n_q);

		fe_values[etar].get_function_gradients(solution,detar_n_q);
		fe_values[etatheta].get_function_gradients(solution,detatheta_n_q);
		fe_values[etaz].get_function_gradients(solution,detaz_n_q);


		const auto &x_q0 = fe_values.quadrature_point(0);

		const auto &x0cell_iterator = GridTools::find_active_cell_around_point(x0dof_handler,x_q0);
		fe_values_x0.reinit(x0cell_iterator);
		unsigned int x0cell_index = x0cell_iterator->active_cell_index();


		fe_values_x0[r].get_function_values(x0solution,r_q);
		fe_values_x0[z].get_function_values(x0solution,z_q);
		fe_values_x0[dr].get_function_values(x0solution,dr_q);
		fe_values_x0[dz].get_function_values(x0solution,dz_q);
		fe_values_x0[dr].get_function_gradients(x0solution,ddr_q);
		fe_values_x0[dz].get_function_gradients(x0solution,ddz_q);
		fe_values_x0[r].get_function_values(R0solution,R_q);
		fe_values_x0[z].get_function_values(R0solution,Phi_q);
		fe_values_x0[lr].get_function_values(R0solution,dPhi_q);

		//




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
			const auto &x_q = fe_values.quadrature_point(q_index);





			// Deformed configuration quantities
			Tensor<1,3> x0_q;
			x0_q[0] = r_q[q_index]; x0_q[2] = z_q[q_index];
			Tensor<1,3> dx0dS_q;
			dx0dS_q[0] = dr_q[q_index]; dx0dS_q[2] = dz_q[q_index];
			Tensor<1,3> dx0dtheta_q;
			dx0dtheta_q[1] = r_q[q_index];

			Tensor<1,3> n0_q = normalize(cross_product(dx0dS_q,dx0dtheta_q));


			Tensor<1,3> ddx0ddS_q;
			ddx0ddS_q[0] = ddr_q[q_index][0]; ddx0ddS_q[2] = ddz_q[q_index][0];
			Tensor<1,3> ddx0dSdtheta_q;
			ddx0dSdtheta_q[1] = dr_q[q_index];
			Tensor<1,3> ddx0ddtheta_q;
			ddx0ddtheta_q[0] = -r_q[q_index];

			double Rref = R_q[q_index];

			//
			Tensor<1,3> v_n_q;
			v_n_q[0] = vr_n_q[q_index]; v_n_q[1] = vtheta_n_q[q_index]; v_n_q[2] = vz_n_q[q_index];
			Tensor<1,3> w_n_q;
			w_n_q[0] = wr_n_q[q_index]; w_n_q[1] = wtheta_n_q[q_index]; w_n_q[2] = wz_n_q[q_index];
			Tensor<1,3> xi_n_q;
			xi_n_q[0] = xir_n_q[q_index]; xi_n_q[1] = xitheta_n_q[q_index]; xi_n_q[2] = xiz_n_q[q_index];
			Tensor<1,3> eta_n_q;
			eta_n_q[0] = etar_n_q[q_index]; eta_n_q[1] = etatheta_n_q[q_index]; eta_n_q[2] = etaz_n_q[q_index];


			Tensor<1,3> dv_n_q;
			dv_n_q[0] = dvr_n_q[q_index][0]; dv_n_q[1] = dvtheta_n_q[q_index][0]; dv_n_q[2] = dvz_n_q[q_index][0];
			Tensor<1,3> dw_n_q;
			dw_n_q[0] = dwr_n_q[q_index][0]; dw_n_q[1] = dwtheta_n_q[q_index][0]; dw_n_q[2] = dwz_n_q[q_index][0];
			Tensor<1,3> dxi_n_q;
			dxi_n_q[0] = dxir_n_q[q_index][0]; dxi_n_q[1] = dxitheta_n_q[q_index][0]; dxi_n_q[2] = dxiz_n_q[q_index][0];
			Tensor<1,3> deta_n_q;
			deta_n_q[0] = detar_n_q[q_index][0]; deta_n_q[1] = detatheta_n_q[q_index][0]; deta_n_q[2] = detaz_n_q[q_index][0];
			//


			const Tensor<1,3> dudS_cos_q = xi_n_q;
			const Tensor<1,3> dudS_sin_q = eta_n_q;
			const Tensor<1,3> dudtheta_cos_q = cross_product(ez,v_n_q) + nval*w_n_q;
			const Tensor<1,3> dudtheta_sin_q = cross_product(ez,w_n_q) - nval*v_n_q;

			const Tensor<1,3> dduddS_cos_q = dxi_n_q;
			const Tensor<1,3> dduddS_sin_q = deta_n_q;
			const Tensor<1,3> ddudSdtheta_cos_q = cross_product(ez,xi_n_q) + nval*eta_n_q;
			const Tensor<1,3> ddudSdtheta_sin_q = cross_product(ez,eta_n_q) - nval*xi_n_q;
			const Tensor<1,3> dduddtheta_cos_q = cross_product(ez,dudtheta_cos_q) + nval*cross_product(ez,w_n_q) - nval*nval*v_n_q;
			const Tensor<1,3> dduddtheta_sin_q = cross_product(ez,dudtheta_sin_q) - nval*cross_product(ez,v_n_q) - nval*nval*w_n_q;



			Tensor<2,2> a0c;
			a0c[0][0] = dx0dS_q*dx0dS_q;
			a0c[0][1] = dx0dtheta_q*dx0dS_q;
			a0c[1][0] = a0c[0][1];
			a0c[1][1] = dx0dtheta_q*dx0dtheta_q;

			Tensor<2,2> a0C = inverse2x2(a0c);

			Tensor<2,2> Ac;
			Ac[0][0] = 1.0;
			Ac[1][1] = Rref*Rref;

			Tensor<2,2> AC;
			AC[0][0] = 1.0;
			AC[1][1] = 1.0/(Rref*Rref);

			Tensor<2,2> Pt;
			Pt[0][0] = 1.0;
			Pt[1][1] = Rref;

			double hsc = 12.0*pow(h,2.0)/(1.0 - nu*nu);


			Tensor<2,2> CovariantMetric;

			CovariantMetric[0][0] = dx0dS_q*dx0dS_q;
			CovariantMetric[0][1] = dx0dS_q*dx0dtheta_q;
			CovariantMetric[1][0] = CovariantMetric[0][1];
			CovariantMetric[1][1] = dx0dtheta_q*dx0dtheta_q;


			Tensor<2,2> ProjectedMetric = 0.5*Pt*(AC*CovariantMetric*AC - AC)*Pt;
			const double stretch_q = sqrt(pow(dr_q[q_index],2.0) + pow(dz_q[q_index],2.0));


			Tensor<2,2> Covariant2Form;

			Covariant2Form[0][0] = n0_q*ddx0ddS_q;
			Covariant2Form[0][1] = n0_q*ddx0dSdtheta_q;
			Covariant2Form[1][0] = Covariant2Form[0][1];
			Covariant2Form[1][1] = n0_q*ddx0ddtheta_q;


			Tensor<2,2> Projected2Form;
			Projected2Form = Pt*a0C*Covariant2Form*a0C*Pt;

			Tensor<2,2> Reference2Form;
			Reference2Form[0][0] = -dPhi_q[q_index];
			Reference2Form[1][1] = cos(Phi_q[q_index])/Rref;



			Tensor<2,2> InPlane = ProjectedMetric - inplane_a.get_value(cell_index, q_index);
			Tensor<2,2> Bending = Projected2Form - Reference2Form - bending_a.get_value(cell_index, q_index);


			Material_Vector_InPlane[cell_index].set_Params(Emodv,nu, InPlane);
			Material_Vector_Bending[cell_index].set_Params(Emodv, nu, Bending);



			for (const unsigned int i : fe_values.dof_indices())
			{
				//
				const double Vr_n_i_q = fe_values[vr].value(i,q_index);
				const double Vtheta_n_i_q = fe_values[vtheta].value(i,q_index);
				const double Vz_n_i_q = fe_values[vz].value(i,q_index);

				const double Wr_n_i_q = fe_values[wr].value(i,q_index);
				const double Wtheta_n_i_q = fe_values[wtheta].value(i,q_index);
				const double Wz_n_i_q = fe_values[wz].value(i,q_index);

				const double Xir_n_i_q = fe_values[xir].value(i,q_index);
				const double Xitheta_n_i_q = fe_values[xitheta].value(i,q_index);
				const double Xiz_n_i_q = fe_values[xiz].value(i,q_index);

				const double Etar_n_i_q = fe_values[etar].value(i,q_index);
				const double Etatheta_n_i_q = fe_values[etatheta].value(i,q_index);
				const double Etaz_n_i_q = fe_values[etaz].value(i,q_index);

				Tensor<1,3> V_n_i_q;
				V_n_i_q[0] = Vr_n_i_q; V_n_i_q[1] = Vtheta_n_i_q; V_n_i_q[2] = Vz_n_i_q;

				Tensor<1,3> W_n_i_q;
				W_n_i_q[0] = Wr_n_i_q; W_n_i_q[1] = Wtheta_n_i_q; W_n_i_q[2] = Wz_n_i_q;

				Tensor<1,3> Xi_n_i_q;
				Xi_n_i_q[0] = Xir_n_i_q; Xi_n_i_q[1] = Xitheta_n_i_q; Xi_n_i_q[2] = Xiz_n_i_q;

				Tensor<1,3> Eta_n_i_q;
				Eta_n_i_q[0] = Etar_n_i_q; Eta_n_i_q[1] = Etatheta_n_i_q; Eta_n_i_q[2] = Etaz_n_i_q;


				const Tensor<1,1> dVr_n_i_q = fe_values[vr].gradient(i,q_index);
				const Tensor<1,1> dVtheta_n_i_q = fe_values[vtheta].gradient(i,q_index);
				const Tensor<1,1> dVz_n_i_q = fe_values[vz].gradient(i,q_index);

				const Tensor<1,1> dWr_n_i_q = fe_values[wr].gradient(i,q_index);
				const Tensor<1,1> dWtheta_n_i_q = fe_values[wtheta].gradient(i,q_index);
				const Tensor<1,1> dWz_n_i_q = fe_values[wz].gradient(i,q_index);

				const Tensor<1,1> dXir_n_i_q = fe_values[xir].gradient(i,q_index);
				const Tensor<1,1> dXitheta_n_i_q = fe_values[xitheta].gradient(i,q_index);
				const Tensor<1,1> dXiz_n_i_q = fe_values[xiz].gradient(i,q_index);

				const Tensor<1,1> dEtar_n_i_q = fe_values[etar].gradient(i,q_index);
				const Tensor<1,1> dEtatheta_n_i_q = fe_values[etatheta].gradient(i,q_index);
				const Tensor<1,1> dEtaz_n_i_q = fe_values[etaz].gradient(i,q_index);

				Tensor<1,3> dV_n_i_q; ///adding in zeros to deal with the tensor to double conversion
				dV_n_i_q[0] = dVr_n_i_q[0]; dV_n_i_q[1] = dVtheta_n_i_q[0]; dV_n_i_q[2] = dVz_n_i_q[0];

				Tensor<1,3> dW_n_i_q;
				dW_n_i_q[0] = dWr_n_i_q[0]; dW_n_i_q[1] = dWtheta_n_i_q[0]; dW_n_i_q[2] = dWz_n_i_q[0];

				Tensor<1,3> dXi_n_i_q;
				dXi_n_i_q[0] = dXir_n_i_q[0]; dXi_n_i_q[1] = dXitheta_n_i_q[0]; dXi_n_i_q[2] = dXiz_n_i_q[0];

				Tensor<1,3> dEta_n_i_q;
				dEta_n_i_q[0] = dEtar_n_i_q[0]; dEta_n_i_q[1] = dEtatheta_n_i_q[0]; dEta_n_i_q[2] = dEtaz_n_i_q[0];


				const Tensor<1,3> dUdS_cos_i_q = Xi_n_i_q;
				const Tensor<1,3> dUdS_sin_i_q = Eta_n_i_q;
				const Tensor<1,3> dUdtheta_cos_i_q = cross_product(ez,V_n_i_q) + nval*W_n_i_q;
				const Tensor<1,3> dUdtheta_sin_i_q = cross_product(ez,W_n_i_q) - nval*V_n_i_q;

				const Tensor<1,3> ddUddS_cos_i_q = dXi_n_i_q;
				const Tensor<1,3> ddUddS_sin_i_q = dEta_n_i_q;
				const Tensor<1,3> ddUdSdtheta_cos_i_q = cross_product(ez,Xi_n_i_q) + nval*Eta_n_i_q;
				const Tensor<1,3> ddUdSdtheta_sin_i_q = cross_product(ez,Eta_n_i_q) - nval*Xi_n_i_q;
				const Tensor<1,3> ddUddtheta_cos_i_q = cross_product(ez,dUdtheta_cos_i_q) + nval*cross_product(ez,W_n_i_q) - nval*nval*V_n_i_q;
				const Tensor<1,3> ddUddtheta_sin_i_q = cross_product(ez,dUdtheta_sin_i_q) - nval*cross_product(ez,V_n_i_q) - nval*nval*W_n_i_q;


				//

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
				//

				for (const unsigned int j : fe_values.dof_indices()) {


					//
					const double Vr_n_j_q = fe_values[vr].value(j,q_index);
					const double Vtheta_n_j_q = fe_values[vtheta].value(j,q_index);
					const double Vz_n_j_q = fe_values[vz].value(j,q_index);

					const double Wr_n_j_q = fe_values[wr].value(j,q_index);
					const double Wtheta_n_j_q = fe_values[wtheta].value(j,q_index);
					const double Wz_n_j_q = fe_values[wz].value(j,q_index);

					const double Xir_n_j_q = fe_values[xir].value(j,q_index);
					const double Xitheta_n_j_q = fe_values[xitheta].value(j,q_index);
					const double Xiz_n_j_q = fe_values[xiz].value(j,q_index);

					const double Etar_n_j_q = fe_values[etar].value(j,q_index);
					const double Etatheta_n_j_q = fe_values[etatheta].value(j,q_index);
					const double Etaz_n_j_q = fe_values[etaz].value(j,q_index);

					Tensor<1,3> V_n_j_q;
					V_n_j_q[0] = Vr_n_j_q; V_n_j_q[1] = Vtheta_n_j_q; V_n_j_q[2] = Vz_n_j_q;

					Tensor<1,3> W_n_j_q;
					W_n_j_q[0] = Wr_n_j_q; W_n_j_q[1] = Wtheta_n_j_q; W_n_j_q[2] = Wz_n_j_q;

					Tensor<1,3> Xi_n_j_q;
					Xi_n_j_q[0] = Xir_n_j_q; Xi_n_j_q[1] = Xitheta_n_j_q; Xi_n_j_q[2] = Xiz_n_j_q;

					Tensor<1,3> Eta_n_j_q;
					Eta_n_j_q[0] = Etar_n_j_q; Eta_n_j_q[1] = Etatheta_n_j_q; Eta_n_j_q[2] = Etaz_n_j_q;


					const Tensor<1,1> dVr_n_j_q = fe_values[vr].gradient(j,q_index);
					const Tensor<1,1> dVtheta_n_j_q = fe_values[vtheta].gradient(j,q_index);
					const Tensor<1,1> dVz_n_j_q = fe_values[vz].gradient(j,q_index);

					const Tensor<1,1> dWr_n_j_q = fe_values[wr].gradient(j,q_index);
					const Tensor<1,1> dWtheta_n_j_q = fe_values[wtheta].gradient(j,q_index);
					const Tensor<1,1> dWz_n_j_q = fe_values[wz].gradient(j,q_index);

					const Tensor<1,1> dXir_n_j_q = fe_values[xir].gradient(j,q_index);
					const Tensor<1,1> dXitheta_n_j_q = fe_values[xitheta].gradient(j,q_index);
					const Tensor<1,1> dXiz_n_j_q = fe_values[xiz].gradient(j,q_index);

					const Tensor<1,1> dEtar_n_j_q = fe_values[etar].gradient(j,q_index);
					const Tensor<1,1> dEtatheta_n_j_q = fe_values[etatheta].gradient(j,q_index);
					const Tensor<1,1> dEtaz_n_j_q = fe_values[etaz].gradient(j,q_index);

					Tensor<1,3> dV_n_j_q; ///adding in zeros to deal with the tensor to double conversion
					dV_n_j_q[0] = dVr_n_j_q[0]; dV_n_j_q[1] = dVtheta_n_j_q[0]; dV_n_j_q[2] = dVz_n_j_q[0];

					Tensor<1,3> dW_n_j_q;
					dW_n_j_q[0] = dWr_n_j_q[0]; dW_n_j_q[1] = dWtheta_n_j_q[0]; dW_n_j_q[2] = dWz_n_j_q[0];

					Tensor<1,3> dXi_n_j_q;
					dXi_n_j_q[0] = dXir_n_j_q[0]; dXi_n_j_q[1] = dXitheta_n_j_q[0]; dXi_n_j_q[2] = dXiz_n_j_q[0];

					Tensor<1,3> dEta_n_j_q;
					dEta_n_j_q[0] = dEtar_n_j_q[0]; dEta_n_j_q[1] = dEtatheta_n_j_q[0]; dEta_n_j_q[2] = dEtaz_n_j_q[0];


					const Tensor<1,3> dUdS_cos_j_q = Xi_n_j_q;
					const Tensor<1,3> dUdS_sin_j_q = Eta_n_j_q;
					const Tensor<1,3> dUdtheta_cos_j_q = cross_product(ez,V_n_j_q) + nval*W_n_j_q;
					const Tensor<1,3> dUdtheta_sin_j_q = cross_product(ez,W_n_j_q) - nval*V_n_j_q;

					const Tensor<1,3> ddUddS_cos_j_q = dXi_n_j_q;
					const Tensor<1,3> ddUddS_sin_j_q = dEta_n_j_q;
					const Tensor<1,3> ddUdSdtheta_cos_j_q = cross_product(ez,Xi_n_j_q) + nval*Eta_n_j_q;
					const Tensor<1,3> ddUdSdtheta_sin_j_q = cross_product(ez,Eta_n_j_q) - nval*Xi_n_j_q;
					const Tensor<1,3> ddUddtheta_cos_j_q = cross_product(ez,dUdtheta_cos_j_q) + nval*cross_product(ez,W_n_j_q) - nval*nval*V_n_j_q;
					const Tensor<1,3> ddUddtheta_sin_j_q = cross_product(ez,dUdtheta_sin_j_q) - nval*cross_product(ez,V_n_j_q) - nval*nval*W_n_j_q;


					//

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

					cell_matrix(i,j) += pi*Rref*(BilinearProduct(DE_cos_i_q,Material_Vector_InPlane[cell_index].getddQ2ddF(),DE_cos_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*(BilinearProduct(DE_sin_i_q,Material_Vector_InPlane[cell_index].getddQ2ddF(),DE_sin_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*(Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),DDE_cos_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*(Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),DDE_sin_q))*fe_values.JxW(q_index);


					cell_matrix(i,j) += pi*Rref*hsc*(BilinearProduct(DB_cos_i_q,Material_Vector_Bending[cell_index].getddQ2ddF(),DB_cos_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*hsc*(BilinearProduct(DB_sin_i_q,Material_Vector_Bending[cell_index].getddQ2ddF(),DB_sin_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*hsc*(Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),DDB_cos_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += pi*Rref*hsc*(Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),DDB_sin_q))*fe_values.JxW(q_index);


				}


				/*
				cell_rhs(i) += ((dr_q[q_index][0]*dR_i_q[0] + dz_q[q_index][0]*dZ_i_q[0]) * // phi_i(x_q)
						fe_values.JxW(q_index));            // dx

				cell_rhs(i) += 100.0*((r_q[q_index]-0.5)*R_i_q + (z_q[q_index] - 0.75)*Z_i_q)*fe_values.JxW(q_index);
				 */






			}
		}

		cell->get_dof_indices(local_dof_indices);

		//constraints.distribute_local_to_global( cell_matrix, cell_rhs, local_dof_indices,system_matrix,system_rhs);

		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < local_dof_indices.size(); i++) {


			for (unsigned int j = 0; j < local_dof_indices.size(); j++) {
				system_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[i][j]);
			}
		}
	}

	constraints.condense(system_matrix);

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



	std::vector<Tensor<1,1>> dvr_n_q(n_q_points);
	std::vector<Tensor<1,1>> dvtheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> dvz_n_q(n_q_points);

	std::vector<Tensor<1,1>> dwr_n_q(n_q_points);
	std::vector<Tensor<1,1>> dwtheta_n_q(n_q_points);
	std::vector<Tensor<1,1>> dwz_n_q(n_q_points);

	std::vector<double> xir_n_q(n_q_points);
	std::vector<double> xitheta_n_q(n_q_points);
	std::vector<double> xiz_n_q(n_q_points);

	std::vector<double> etar_n_q(n_q_points);
	std::vector<double> etatheta_n_q(n_q_points);
	std::vector<double> etaz_n_q(n_q_points);







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



		//
		const FEValuesExtractors::Scalar vr(0);
		const FEValuesExtractors::Scalar vtheta(1);
		const FEValuesExtractors::Scalar vz(2);
		const FEValuesExtractors::Scalar wr(3);
		const FEValuesExtractors::Scalar wtheta(4);
		const FEValuesExtractors::Scalar wz(5);

		const FEValuesExtractors::Scalar xir(6);
		const FEValuesExtractors::Scalar xitheta(7);
		const FEValuesExtractors::Scalar xiz(8);
		const FEValuesExtractors::Scalar etar(9);
		const FEValuesExtractors::Scalar etatheta(10);
		const FEValuesExtractors::Scalar etaz(11);

		//

		//

		fe_values[xir].get_function_values(solution,xir_n_q);
		fe_values[xitheta].get_function_values(solution,xitheta_n_q);
		fe_values[xiz].get_function_values(solution,xiz_n_q);

		fe_values[etar].get_function_values(solution,etar_n_q);
		fe_values[etatheta].get_function_values(solution,etatheta_n_q);
		fe_values[etaz].get_function_values(solution,etaz_n_q);


		fe_values[vr].get_function_gradients(solution,dvr_n_q);
		fe_values[vtheta].get_function_gradients(solution,dvtheta_n_q);
		fe_values[vz].get_function_gradients(solution,dvz_n_q);

		fe_values[wr].get_function_gradients(solution,dwr_n_q);
		fe_values[wtheta].get_function_gradients(solution,dwtheta_n_q);
		fe_values[wz].get_function_gradients(solution,dwz_n_q);


		//



		for (const unsigned int q_index : fe_values.quadrature_point_indices()){



			const auto &x_q = fe_values.quadrature_point(q_index);



			for (const unsigned int i : fe_values.dof_indices())
			{
				//
				const double Vr_n_i_q = fe_values[vr].value(i,q_index);
				const double Vtheta_n_i_q = fe_values[vtheta].value(i,q_index);
				const double Vz_n_i_q = fe_values[vz].value(i,q_index);

				const double Wr_n_i_q = fe_values[wr].value(i,q_index);
				const double Wtheta_n_i_q = fe_values[wtheta].value(i,q_index);
				const double Wz_n_i_q = fe_values[wz].value(i,q_index);

				const double Xir_n_i_q = fe_values[xir].value(i,q_index);
				const double Xitheta_n_i_q = fe_values[xitheta].value(i,q_index);
				const double Xiz_n_i_q = fe_values[xiz].value(i,q_index);

				const double Etar_n_i_q = fe_values[etar].value(i,q_index);
				const double Etatheta_n_i_q = fe_values[etatheta].value(i,q_index);
				const double Etaz_n_i_q = fe_values[etaz].value(i,q_index);

				Tensor<1,3> V_n_i_q;
				V_n_i_q[0] = Vr_n_i_q; V_n_i_q[1] = Vtheta_n_i_q; V_n_i_q[2] = Vz_n_i_q;

				Tensor<1,3> W_n_i_q;
				W_n_i_q[0] = Wr_n_i_q; W_n_i_q[1] = Wtheta_n_i_q; W_n_i_q[2] = Wz_n_i_q;

				Tensor<1,3> Xi_n_i_q;
				Xi_n_i_q[0] = Xir_n_i_q; Xi_n_i_q[1] = Xitheta_n_i_q; Xi_n_i_q[2] = Xiz_n_i_q;

				Tensor<1,3> Eta_n_i_q;
				Eta_n_i_q[0] = Etar_n_i_q; Eta_n_i_q[1] = Etatheta_n_i_q; Eta_n_i_q[2] = Etaz_n_i_q;


				const Tensor<1,1> dVr_n_i_q = fe_values[vr].gradient(i,q_index);
				const Tensor<1,1> dVtheta_n_i_q = fe_values[vtheta].gradient(i,q_index);
				const Tensor<1,1> dVz_n_i_q = fe_values[vz].gradient(i,q_index);

				const Tensor<1,1> dWr_n_i_q = fe_values[wr].gradient(i,q_index);
				const Tensor<1,1> dWtheta_n_i_q = fe_values[wtheta].gradient(i,q_index);
				const Tensor<1,1> dWz_n_i_q = fe_values[wz].gradient(i,q_index);

				const Tensor<1,1> dXir_n_i_q = fe_values[xir].gradient(i,q_index);
				const Tensor<1,1> dXitheta_n_i_q = fe_values[xitheta].gradient(i,q_index);
				const Tensor<1,1> dXiz_n_i_q = fe_values[xiz].gradient(i,q_index);

				const Tensor<1,1> dEtar_n_i_q = fe_values[etar].gradient(i,q_index);
				const Tensor<1,1> dEtatheta_n_i_q = fe_values[etatheta].gradient(i,q_index);
				const Tensor<1,1> dEtaz_n_i_q = fe_values[etaz].gradient(i,q_index);

				Tensor<1,3> dV_n_i_q; ///adding in zeros to deal with the tensor to double conversion
				dV_n_i_q[0] = dVr_n_i_q[0]; dV_n_i_q[1] = dVtheta_n_i_q[0]; dV_n_i_q[2] = dVz_n_i_q[0];

				Tensor<1,3> dW_n_i_q;
				dW_n_i_q[0] = dWr_n_i_q[0]; dW_n_i_q[1] = dWtheta_n_i_q[0]; dW_n_i_q[2] = dWz_n_i_q[0];

				Tensor<1,3> dXi_n_i_q;
				dXi_n_i_q[0] = dXir_n_i_q[0]; dXi_n_i_q[1] = dXitheta_n_i_q[0]; dXi_n_i_q[2] = dXiz_n_i_q[0];

				Tensor<1,3> dEta_n_i_q;
				dEta_n_i_q[0] = dEtar_n_i_q[0]; dEta_n_i_q[1] = dEtatheta_n_i_q[0]; dEta_n_i_q[2] = dEtaz_n_i_q[0];





				for (const unsigned int j : fe_values.dof_indices()) {


					//
					const double Vr_n_j_q = fe_values[vr].value(j,q_index);
					const double Vtheta_n_j_q = fe_values[vtheta].value(j,q_index);
					const double Vz_n_j_q = fe_values[vz].value(j,q_index);

					const double Wr_n_j_q = fe_values[wr].value(j,q_index);
					const double Wtheta_n_j_q = fe_values[wtheta].value(j,q_index);
					const double Wz_n_j_q = fe_values[wz].value(j,q_index);

					const double Xir_n_j_q = fe_values[xir].value(j,q_index);
					const double Xitheta_n_j_q = fe_values[xitheta].value(j,q_index);
					const double Xiz_n_j_q = fe_values[xiz].value(j,q_index);

					const double Etar_n_j_q = fe_values[etar].value(j,q_index);
					const double Etatheta_n_j_q = fe_values[etatheta].value(j,q_index);
					const double Etaz_n_j_q = fe_values[etaz].value(j,q_index);

					Tensor<1,3> V_n_j_q;
					V_n_j_q[0] = Vr_n_j_q; V_n_j_q[1] = Vtheta_n_j_q; V_n_j_q[2] = Vz_n_j_q;

					Tensor<1,3> W_n_j_q;
					W_n_j_q[0] = Wr_n_j_q; W_n_j_q[1] = Wtheta_n_j_q; W_n_j_q[2] = Wz_n_j_q;

					Tensor<1,3> Xi_n_j_q;
					Xi_n_j_q[0] = Xir_n_j_q; Xi_n_j_q[1] = Xitheta_n_j_q; Xi_n_j_q[2] = Xiz_n_j_q;

					Tensor<1,3> Eta_n_j_q;
					Eta_n_j_q[0] = Etar_n_j_q; Eta_n_j_q[1] = Etatheta_n_j_q; Eta_n_j_q[2] = Etaz_n_j_q;


					const Tensor<1,1> dVr_n_j_q = fe_values[vr].gradient(j,q_index);
					const Tensor<1,1> dVtheta_n_j_q = fe_values[vtheta].gradient(j,q_index);
					const Tensor<1,1> dVz_n_j_q = fe_values[vz].gradient(j,q_index);

					const Tensor<1,1> dWr_n_j_q = fe_values[wr].gradient(j,q_index);
					const Tensor<1,1> dWtheta_n_j_q = fe_values[wtheta].gradient(j,q_index);
					const Tensor<1,1> dWz_n_j_q = fe_values[wz].gradient(j,q_index);

					const Tensor<1,1> dXir_n_j_q = fe_values[xir].gradient(j,q_index);
					const Tensor<1,1> dXitheta_n_j_q = fe_values[xitheta].gradient(j,q_index);
					const Tensor<1,1> dXiz_n_j_q = fe_values[xiz].gradient(j,q_index);

					const Tensor<1,1> dEtar_n_j_q = fe_values[etar].gradient(j,q_index);
					const Tensor<1,1> dEtatheta_n_j_q = fe_values[etatheta].gradient(j,q_index);
					const Tensor<1,1> dEtaz_n_j_q = fe_values[etaz].gradient(j,q_index);

					Tensor<1,3> dV_n_j_q; ///adding in zeros to deal with the tensor to double conversion
					dV_n_j_q[0] = dVr_n_j_q[0]; dV_n_j_q[1] = dVtheta_n_j_q[0]; dV_n_j_q[2] = dVz_n_j_q[0];

					Tensor<1,3> dW_n_j_q;
					dW_n_j_q[0] = dWr_n_j_q[0]; dW_n_j_q[1] = dWtheta_n_j_q[0]; dW_n_j_q[2] = dWz_n_j_q[0];

					Tensor<1,3> dXi_n_j_q;
					dXi_n_j_q[0] = dXir_n_j_q[0]; dXi_n_j_q[1] = dXitheta_n_j_q[0]; dXi_n_j_q[2] = dXiz_n_j_q[0];

					Tensor<1,3> dEta_n_j_q;
					dEta_n_j_q[0] = dEtar_n_j_q[0]; dEta_n_j_q[1] = dEtatheta_n_j_q[0]; dEta_n_j_q[2] = dEtaz_n_j_q[0];


					cell_matrix(i,j) += 10000.*(dV_n_i_q - Xi_n_i_q)*Xi_n_j_q*fe_values.JxW(q_index);
					cell_matrix(i,j) += 10000.*(dW_n_i_q - Eta_n_i_q)*Eta_n_j_q*fe_values.JxW(q_index);




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


void ElasticProblem::construct_reduced_mappings(){


	x_global_to_reduced.clear();
	x_reduced_to_global.clear();
	xi_global_to_reduced.clear();
	xi_reduced_to_global.clear();


	const int ndofs = dof_handler.n_dofs();

	x_global_to_reduced.resize(ndofs);
	xi_global_to_reduced.resize(ndofs);
	// Constraint stuff
	//Creating masks for degrees of freedom
	////////////////////
	std::vector<bool> vr_components = create_bool_vector(12,0);
	ComponentMask vr_mask(vr_components);

	std::vector<bool> is_vr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vr_mask, is_vr_comp);
	////////////////////
	std::vector<bool> vtheta_components = create_bool_vector(12,1);
	ComponentMask vtheta_mask(vtheta_components);

	std::vector<bool> is_vtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vtheta_mask, is_vtheta_comp);
	/////////////////////////
	std::vector<bool> vz_components = create_bool_vector(12,2);
	ComponentMask vz_mask(vz_components);

	std::vector<bool> is_vz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vz_mask, is_vz_comp);

	////////////////
	std::vector<bool> wr_components = create_bool_vector(12,3);
	ComponentMask wr_mask(wr_components);

	std::vector<bool> is_wr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wr_mask, is_wr_comp);
	///////////////////
	std::vector<bool> wtheta_components = create_bool_vector(12,4);
	ComponentMask wtheta_mask(wtheta_components);

	std::vector<bool> is_wtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wtheta_mask, is_wtheta_comp);
	//////////////////
	std::vector<bool> wz_components = create_bool_vector(12,5);
	ComponentMask wz_mask(wz_components);

	std::vector<bool> is_wz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wz_mask, is_wz_comp);
	//////////////////////
	std::vector<bool> xir_components = create_bool_vector(12,6);
	ComponentMask xir_mask(xir_components);

	std::vector<bool> is_xir_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xir_mask, is_xir_comp);
	/////////////////////
	std::vector<bool> xitheta_components = create_bool_vector(12,7);
	ComponentMask xitheta_mask(xitheta_components);

	std::vector<bool> is_xitheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xitheta_mask, is_xitheta_comp);
	//////////////////////
	std::vector<bool> xiz_components = create_bool_vector(12,8);
	ComponentMask xiz_mask(xiz_components);

	std::vector<bool> is_xiz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, xiz_mask, is_xiz_comp);
	//////////////////////
	std::vector<bool> etar_components = create_bool_vector(12,9);
	ComponentMask etar_mask(etar_components);

	std::vector<bool> is_etar_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etar_mask, is_etar_comp);
	//////////////////////
	std::vector<bool> etatheta_components = create_bool_vector(12,10);
	ComponentMask etatheta_mask(etatheta_components);

	std::vector<bool> is_etatheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etatheta_mask, is_etatheta_comp);
	//////////////////////
	std::vector<bool> etaz_components = create_bool_vector(12,11);
	ComponentMask etaz_mask(etaz_components);

	std::vector<bool> is_etaz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, etaz_mask, is_etaz_comp);



	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	unsigned int xcntr = 0;
	unsigned int xicntr = 0;
	nxdofs = 0;
	nxidofs = 0;



	for (unsigned int i = 0; i < ndofs; i++) {

		if (is_vr_comp[i] || is_vtheta_comp[i] || is_vz_comp[i] || is_wr_comp[i] || is_wtheta_comp[i] || is_wz_comp[i] ){
			x_reduced_to_global.push_back(i);
			x_global_to_reduced[i] = xcntr;
			xcntr++;
			nxdofs++;
		} else {
			x_global_to_reduced[i] = -1;
		}

		if (is_xir_comp[i] || is_xitheta_comp[i] || is_xiz_comp[i] || is_etar_comp[i] || is_etatheta_comp[i] || is_etaz_comp[i] ){
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
	Eigen::MatrixXd Cxix(nxidofs,nxdofs); Cxix.setZero();
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
			if(fabs(entry->value()) > 1e-12 ){
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
				if(fabs(entry->value()) > 1e-12 ){
					//       system_matrix_petsc.set(row, entry->column(),entry->value());
					if (x_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){

						Cxix(xi_global_to_reduced[global_column] ,x_global_to_reduced[global_row]) += entry->value();
					} else if (xi_global_to_reduced[global_row] > -1 && xi_global_to_reduced[global_column] > -1){
						Cxixi(xi_global_to_reduced[global_row] ,xi_global_to_reduced[global_column]) += entry->value();
					}
				}

			}
		}



		XiProj = -1.0*Cxixi.inverse()*Cxix;

		for (unsigned int i = 0; i < XiProj.rows(); i++){
			for (unsigned int j = 0; j < XiProj.cols(); j++){
				if (fabs(XiProj(i,j)) < 1e-10){
					//XiProj(i,j) = 0.0;
				}

			}
		}
		XiProjtrans = XiProj.transpose();

	}
	std::cout << "Calculating Ktilde" << std::endl;

	//std::thread th1(&calc_Kx1,this,Kxxi,XiProj);
	//calc_Kx1(&Kx1,&Kxxi,&XiProj);


	/*
	{


		//calc_Kx1(&Kx1,&Kxxi,&XiProj);
		std::future<void> resultFromT1 = std::async(std::launch::async,ElasticProblem::calc_MMult,&Kx1,&Kxxi,&XiProj);

		std::cout << "1: Calculating Kxxi*XiR" << std::endl;
		//double dv = 1.2;
		//std::future<double> resultFromT1 = std::async(std::launch::async, [](double * m1){return 2.0*(*m1);},&dv);
		Eigen::MatrixXd K2Temp,K3Temp,K4Temp,K5Temp;
		std::future<void> resultFromT2 = std::async(std::launch::async,ElasticProblem::calc_MMult,&K2Temp,&XiProjtrans,&Kxixi);
		std::cout << "2: Calculating Xir*Kxixi" << std::endl;
		resultFromT1.get();
		resultFromT2.get();

		std::future<void> resultFromT3 = std::async(std::launch::async,ElasticProblem::calc_MAdd,&K3Temp, &Kxx, &Kx1, 1.0, 2.0);

		std::cout << "3: Calculating Kxx + 2*Kxxi*XiR" << std::endl;
		std::future<void> resultFromT4 = std::async(std::launch::async,ElasticProblem::calc_MMult,&K4Temp,&K2Temp,&XiProj);
		std::cout << "4: Calculating XiR*Kxixi*XiR" << std::endl;
		//calc_Kx2(Kxixi);
		resultFromT3.get();
		resultFromT4.get();
		K5Temp = K3Temp + K4Temp;
		Ktilde = 0.5*(K5Temp + K5Temp.transpose());


	}
	 */
	Eigen::MatrixXd Ktildeasym = Kxx + Kxxi*XiProj + XiProjtrans*Kxix + XiProjtrans*Kxixi*XiProj;


	//Ktilde = 0.5*(Ktildeasym + Ktildeasym.transpose());
	Ktilde.setZero();
	Ktilde = Ktildeasym;


	for (unsigned int i = 0; i < Ktilde.rows(); i++) {
		for (unsigned int j = 0; j < Ktilde.cols(); j++){
			/*
			if (fabs(Ktilde(i,j)) < 1e-10){
				Ktilde(i,j) = 0.0;
			} else {
				KtildeLA(i,j) = Ktilde(i,j);
			}
			*/

			KtildeLA(i,j) = Ktilde(i,j);
		}
	}

	//std::string outname = "matrices/ktildeout" + std::to_string((int)nval) + "_" + std::to_string(timestep) + ".csv";
	//writeToCSVfile(outname,Ktilde);
	//Ktilde = Kxx + Kx1 + Kx2;

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
			inplane_a.set_value(cell_index, q_index, current_inplane_a);
			bending_a.set_value(cell_index,q_index,current_bending_a);
		}


	}
}


void ElasticProblem::calc_MMult(Eigen::MatrixXd * Kxout, const Eigen::MatrixXd * Kxxi, const Eigen::MatrixXd * XiR){
	*Kxout = (*Kxxi)*(*XiR);
}
void ElasticProblem::calc_MAdd(Eigen::MatrixXd * Kxout, const Eigen::MatrixXd * Kxxi, const Eigen::MatrixXd * XiR,const double a, const double b){
	*Kxout = a*(*Kxxi) + b*(*XiR);
}
void ElasticProblem::calc_Kx2(const Eigen::MatrixXd & Kxixi){
	Kx2 = XiProjtrans*Kxixi*XiProj;
}
Vector<double> ElasticProblem::calculate_stability(int ts){
	std::cout << "Fourier Mode " << (int) nval << std::endl;
	Vector<double> eigenvalues;
	FullMatrix<double> eigenvectors;
	Vector<double> evalsout(10);
	double tol = 1.0e-13;

	KtildeLA.compute_eigenvalues_symmetric(-Emodv*1.0, Emodv*1000.0*h, tol, eigenvalues, eigenvectors);
	std::cout << "Em rows: " << eigenvectors.n_rows() << "Em cols: " <<eigenvectors.n_cols() << std::endl;
	for (unsigned int i = 0; i < eigenvalues.size(); i++){
		//std::cout << eigenvalues[i] << std::endl;
	}

	for (unsigned int i = 0; i < evalsout.size(); i++){

		evalsout[i] = eigenvalues[i]/Emodv;
		std::cout << "Eigval: " << evalsout[i] << std::endl;
		if (evalsout[i] < -1.e-7 ){
			Vector<double> evout(eigenvectors.n_rows());
			for (unsigned int j = 0; j < eigenvectors.n_rows(); j++){
				evout[j] = eigenvectors[j][i]/Emodv;
			}

			for (unsigned int kk = 0; kk < evout.size(); kk++){
				evout[kk] = evout[kk]/sqrt(evout.norm_sqr());
			}

			savestate(evout,ts);
		}
	}

	return evalsout;


}


void ElasticProblem::savestate(Vector<double> & vcomp, int index){
	const int ndofs = dof_handler.n_dofs();


	//Creating masks for degrees of freedom
	////////////////////
	std::vector<bool> vr_components = create_bool_vector(12,0);
	ComponentMask vr_mask(vr_components);

	std::vector<bool> is_vr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vr_mask, is_vr_comp);
	////////////////////
	std::vector<bool> vtheta_components = create_bool_vector(12,1);
	ComponentMask vtheta_mask(vtheta_components);

	std::vector<bool> is_vtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vtheta_mask, is_vtheta_comp);
	/////////////////////////
	std::vector<bool> vz_components = create_bool_vector(12,2);
	ComponentMask vz_mask(vz_components);

	std::vector<bool> is_vz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, vz_mask, is_vz_comp);

	////////////////
	std::vector<bool> wr_components = create_bool_vector(12,3);
	ComponentMask wr_mask(wr_components);

	std::vector<bool> is_wr_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wr_mask, is_wr_comp);
	///////////////////
	std::vector<bool> wtheta_components = create_bool_vector(12,4);
	ComponentMask wtheta_mask(wtheta_components);

	std::vector<bool> is_wtheta_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wtheta_mask, is_wtheta_comp);
	//////////////////
	std::vector<bool> wz_components = create_bool_vector(12,5);
	ComponentMask wz_mask(wz_components);

	std::vector<bool> is_wz_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, wz_mask, is_wz_comp);





	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);


	std::string strout = "_" + std::to_string(index) + ".csv";
	std::ofstream vrfile;
	vrfile.open("unstablemodes/vr" + strout);
	vrfile << "S,val\n";

	std::ofstream vthetafile;
	vthetafile.open("unstablemodes/vtheta" + strout);
	vthetafile << "S,val\n";
	std::ofstream vzfile;
	vzfile.open("unstablemodes/vz" + strout);
	vzfile << "S,val\n";
	std::ofstream wrfile;
	wrfile.open("unstablemodes/wr" + strout);
	wrfile << "S,val\n";
	std::ofstream wthetafile;
	wthetafile.open("unstablemodes/wtheta" + strout);
	wthetafile << "S,val\n";
	std::ofstream wzfile;
	wzfile.open("unstablemodes/wz" + strout);
	wzfile << "S,val\n";
	for (int i = 0; i < vcomp.size(); i++){
		int iglobal = x_reduced_to_global[i];
		if (is_vr_comp[iglobal]) {
			vrfile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		} else if (is_vtheta_comp[iglobal]){
			vthetafile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		} else if (is_vz_comp[iglobal]){
			vzfile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		} else if (is_wr_comp[iglobal]){
			wrfile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		} else if (is_wtheta_comp[iglobal]){
			wthetafile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		} else if (is_wz_comp[iglobal]){
			wzfile << support_points[iglobal][0] << "," << vcomp[i] << '\n';
		}

	}

	vrfile.close();
	vthetafile.close();
	vzfile.close();
	wrfile.close();
	wthetafile.close();
	wzfile.close();






}


void ElasticProblem::save_stability(std::vector<std::vector<Vector<double>>> eall){
	int nsize = eall[0].size();
	int tsize = eall.size();

	for (unsigned int in = 0; in < nsize; in++){
		std::string strout = "evals" + std::to_string(in+1) + ".csv";
		std::ofstream rfile;
		rfile.open(strout);
		//std::cout << "in: " <<in << std::endl;
		for (unsigned int it = 0; it < tsize; it++){
			//std::cout << "it: " << it << std::endl;
			for (unsigned int nv = 0; nv < eall[it][in].size(); nv++){
				//std::cout << "nv: " << nv << std::endl;
				rfile << std::setprecision(17) << eall[it][in][nv];
				if (nv < (eall[it][in].size() - 1)){
					rfile << ',';
				}
			}
			rfile << '\n';

		}

		rfile.close();
	}
}

std::vector<bool> ElasticProblem::create_bool_vector(int nsize,int ntrue){
	std::vector<bool> vout(nsize);

	for (unsigned int i = 0; i < nsize; i++){
		if (i==ntrue){
			vout[i] = true;
		} else {
			vout[i] = false;
		}
	}
	return vout;
}

std::vector<double> ElasticProblem::linspace(double start_in, double end_in, int num_in){


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
template<typename T>
T ElasticProblem::max_val(std::vector<T> & vin){
	T maxval = vin[0];

	for (unsigned int i = 0; i < vin.size(); i++){
		if (vin[i] > maxval){
			maxval = vin[i];
		}
	}
	return maxval;
}

Tensor<2,2> ElasticProblem::inverse2x2(const Tensor<2,2> & input){
	double detinput = input[0][0]*input[1][1] - input[0][1]*input[1][0];
	Tensor<2,2> output;
	output[0][0] = input[1][1]/detinput;
	output[0][1] = -1.0*input[0][1]/detinput;
	output[1][0] = -1.0*input[1][0]/detinput;
	output[1][1] = input[0][0]/detinput;
	return output;
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

template<int T>
Tensor<2,T> ElasticProblem::outer_product(const Tensor<1,T> & v1, const Tensor<1,T> & v2){
	Tensor<2,T> Mout;

	for (unsigned int i = 0; i < T; i++) {
		for (unsigned int j = 0; j < T; j++) {
			Mout[i][j] = v1[i]*v2[j];
		}
	}
	return Mout;
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

Tensor<1,3> ElasticProblem::normalize(const Tensor<1,3> & v) {
	return v/v.norm();
}

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

	double r0 = 0.3;

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		unsigned int cell_index = cell->active_cell_index();


		Reference_Configuration_Vec[cell_index].resize(n_q_points);

		for (const unsigned int q_index : fe_values.quadrature_point_indices()){
			const auto &x_q = fe_values.quadrature_point(q_index);
			Reference_Configuration_Vec[cell_index][q_index].set_deformation_param(1.0);
			Reference_Configuration_Vec[cell_index][q_index].set_R0(r0);
			Reference_Configuration_Vec[cell_index][q_index].set_point(x_q[0]);

		}


	}
}

// Commented out because this is causing the compiler error with fesystem
void ElasticProblem::load_state(unsigned int indx){
	if (indx==0){


		std::string triang_file_name = "../../CylindricalSystem/build/de/triang.tridat";
		std::ifstream triang_in_u(triang_file_name.c_str());
		boost::archive::text_iarchive triang_ar_u(triang_in_u);
		x0triangulation.load(triang_ar_u,1);

		x0dof_handler.distribute_dofs(x0fe);

		std::string dof_file_name = "../../CylindricalSystem/build/de/dof_handler.dofdat";
		std::ifstream dof_in_u(dof_file_name.c_str());
		boost::archive::text_iarchive dof_ar_u(dof_in_u);
		x0dof_handler.load(dof_ar_u,1);




		x0solution.reinit(x0dof_handler.n_dofs());
		R0solution.reinit(x0dof_handler.n_dofs());

		std::string R0_file_name = "../../CylindricalSystem/build/de/R0.soldat";
		std::ifstream R0_in_u(R0_file_name.c_str());
		boost::archive::text_iarchive R0_ar_u(R0_in_u);
		R0solution.load(R0_ar_u,1);

	}

	std::string sol_file_name = "../../CylindricalSystem/build/de/solution_" + std::to_string(indx) + ".soldat";
	std::ifstream sol_in_u(sol_file_name.c_str());
	boost::archive::text_iarchive sol_ar_u(sol_in_u);
	x0solution.load(sol_ar_u,1);


}

void ElasticProblem::writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
	std::ofstream file(name.c_str());

	for(int  i = 0; i < matrix.rows(); i++){
		for(int j = 0; j < matrix.cols(); j++){
			std::string str = boost::lexical_cast<std::string>(matrix(i,j));
			if(j+1 == matrix.cols()){
				file<<str;
			}else{
				file<<str<<',';
			}
		}
		file<<'\n';
	}
}

void ElasticProblem::load_external_params(){
	std::string str = "../../CylindricalSystem/build/params.dat";
	char *cstr = new char[str.length() + 1];
	strcpy(cstr, str.c_str());
	rf.readInputFile(cstr, dat);
}

void ElasticProblem::run()
{
	std::cout << "Solving problem in " << DIM << " space dimensions."
			<< std::endl;

	load_external_params();
	make_grid();
	setup_system();
	setup_constraints();
	solve_path();
}


}

#endif // ELASTIC_PROBLEM_CC_

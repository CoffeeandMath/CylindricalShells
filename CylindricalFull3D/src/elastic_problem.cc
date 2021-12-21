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
: fe(FESystem<DIM>(FE_Q<DIM>(shapefunctiondeg), 3),1,
		FESystem<DIM>(FE_Q<DIM>(shapefunctiondeg),2),3,
		FESystem<DIM>(FE_Q<DIM>(shapefunctiondeg),2),3)
, dof_handler(triangulation){}


void ElasticProblem::solve_path(){


	int Nsteps = 100;
	double defmagmax = 1.0;

	std::vector<double> defmagvec = linspace(0.0,defmagmax,Nsteps);
	Tensor<2,2> Btemp;
	for (unsigned int i = 0; i < Nsteps; i++){
		Btemp[1][1] = 0.2*defmagvec[i];
		Btemp[0][0] = defmagvec[i];
		Applied_Bending.set_all_values(Btemp);
		newton_raphson();
		output_results(i);
	}






}
void ElasticProblem::assemble_system(){
	system_matrix = 0;
	system_rhs = 0;

	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);

	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);


	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);

	std::vector<double> x1_n_q(n_q_points);
	std::vector<double> x2_n_q(n_q_points);
	std::vector<double> x3_n_q(n_q_points);

	std::vector<double> x1old_n_q(n_q_points);
	std::vector<double> x2old_n_q(n_q_points);
	std::vector<double> x3old_n_q(n_q_points);

	std::vector<Tensor<1,2>> Dx1_n_q(n_q_points);
	std::vector<Tensor<1,2>> Dx2_n_q(n_q_points);
	std::vector<Tensor<1,2>> Dx3_n_q(n_q_points);

	std::vector<Tensor<1,2>> f1_n_q(n_q_points);
	std::vector<Tensor<1,2>> f2_n_q(n_q_points);
	std::vector<Tensor<1,2>> f3_n_q(n_q_points);

	std::vector<Tensor<2,2>> Df1_n_q(n_q_points);
	std::vector<Tensor<2,2>> Df2_n_q(n_q_points);
	std::vector<Tensor<2,2>> Df3_n_q(n_q_points);

	std::vector<Tensor<1,2>> lambda1_n_q(n_q_points);
	std::vector<Tensor<1,2>> lambda2_n_q(n_q_points);
	std::vector<Tensor<1,2>> lambda3_n_q(n_q_points);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	Tensor<2,2> Acov_q;
	Tensor<2,2> Acontra_q;
	Tensor<2,2> FR_q;
	Tensor<2,2> FRinv_q;
	Tensor<2,2> P_q;

	Tensor<2,2> E_q;
	Tensor<2,2> B_q;
	Tensor<2,3> Identity3D; Identity3D[0][0] = 1.0; Identity3D[1][1] = 1.0; Identity3D[2][2] = 1.0;


	for (const auto &cell : dof_handler.active_cell_iterators()){
		fe_values.reinit(cell);
		cell_matrix = 0;
		cell_rhs    = 0;
		unsigned int cell_index = cell->active_cell_index();

		const FEValuesExtractors::Scalar x1(0);
		const FEValuesExtractors::Scalar x2(1);
		const FEValuesExtractors::Scalar x3(2);
		const FEValuesExtractors::Vector f1(3);
		const FEValuesExtractors::Vector f2(5);
		const FEValuesExtractors::Vector f3(7);
		const FEValuesExtractors::Vector lambda1(9);
		const FEValuesExtractors::Vector lambda2(11);
		const FEValuesExtractors::Vector lambda3(13);



		fe_values[x1].get_function_values(solution,x1_n_q);
		fe_values[x2].get_function_values(solution,x2_n_q);
		fe_values[x3].get_function_values(solution,x3_n_q);
		fe_values[x1].get_function_values(prev_solution,x1old_n_q);
		fe_values[x2].get_function_values(prev_solution,x2old_n_q);
		fe_values[x3].get_function_values(prev_solution,x3old_n_q);
		fe_values[x1].get_function_gradients(solution,Dx1_n_q);
		fe_values[x2].get_function_gradients(solution,Dx2_n_q);
		fe_values[x3].get_function_gradients(solution,Dx3_n_q);
		fe_values[f1].get_function_values(solution,f1_n_q);
		fe_values[f2].get_function_values(solution,f2_n_q);
		fe_values[f3].get_function_values(solution,f3_n_q);
		fe_values[f1].get_function_gradients(solution,Df1_n_q);
		fe_values[f2].get_function_gradients(solution,Df2_n_q);
		fe_values[f3].get_function_gradients(solution,Df3_n_q);
		fe_values[lambda1].get_function_values(solution,lambda1_n_q);
		fe_values[lambda2].get_function_values(solution,lambda2_n_q);
		fe_values[lambda3].get_function_values(solution,lambda3_n_q);

		////////////// Asserting data structures /////////////////
		Tensor<1,3> x_q;
		Tensor<1,3> a1_q;
		Tensor<1,3> a2_q;
		Tensor<1,3> n_q;
		Tensor<2,3> scalednproj_q;
		Tensor<2,2> acov_q;
		Tensor<2,2> acontra_q;
		Tensor<2,2> bcov_q;
		/////////i components//////
		Tensor<1,3> dx_i_q;
		Tensor<1,2> df1_i_q;
		Tensor<1,2> df2_i_q;
		Tensor<1,2> df3_i_q;
		Tensor<2,2> Ddf1_i_q;
		Tensor<2,2> Ddf2_i_q;
		Tensor<2,2> Ddf3_i_q;

		Tensor<1,2> dlambda1_i_q;
		Tensor<1,2> dlambda2_i_q;
		Tensor<1,2> dlambda3_i_q;

		Tensor<1,3> da1_i_q;
		Tensor<1,3> da2_i_q;
		Tensor<2,2> dacov_i_q;
		Tensor<2,2> dacontra_i_q;
		Tensor<1,3> dN_i_q;
		Tensor<1,3> dn_i_q;
		Tensor<2,2> dbcov_i_q;
		Tensor<2,2> dE_i_q;
		Tensor<2,2> dB_i_q;

		/////////////// j components ///////////
		Tensor<1,3> dx_j_q;
		Tensor<1,2> df1_j_q;
		Tensor<1,2> df2_j_q;
		Tensor<1,2> df3_j_q;

		Tensor<2,2> Ddf1_j_q;
		Tensor<2,2> Ddf2_j_q;
		Tensor<2,2> Ddf3_j_q;

		Tensor<1,2> dlambda1_j_q;
		Tensor<1,2> dlambda2_j_q;
		Tensor<1,2> dlambda3_j_q;

		Tensor<1,3> da1_j_q;
		Tensor<1,3> da2_j_q;

		Tensor<2,2> dacov_j_q;
		Tensor<2,2> dacontra_j_q;

		Tensor<1,3> dN_j_q;
		Tensor<1,3> dn_j_q;

		Tensor<2,2> dbcov_j_q;

		Tensor<2,2> dE_j_q;
		Tensor<2,2> dB_j_q;
		////////// ij components ///////////
		Tensor<2,2> ddacov_ij_q;
		Tensor<1,3> ddN_ij_q;
		Tensor<1,3> ddn_ij_q;
		Tensor<2,2> ddbcov_ij_q;

		Tensor<2,2> ddE_ij_q;
		Tensor<2,2> ddB_ij_q;

		/////////// loop temporaries ///////////
		Tensor<1,3> DDdx;
		Tensor<1,3> DDx;
		Tensor<1,3> DDxi;
		Tensor<1,3> DDxj;
		for (const unsigned int q_index : fe_values.quadrature_point_indices()){
			const auto &xi_q = fe_values.quadrature_point(q_index);

			x_q[0] = x1_n_q[q_index]; x_q[1] = x2_n_q[q_index]; x_q[2] = x3_n_q[q_index];


			Acov_q[0][0] = 1.0; Acov_q[0][1] = 0.0; Acov_q[1][0] = 0.0; Acov_q[1][1] = 1.0;
			FR_q = Acov_q;
			Acontra_q = inverse2x2(Acov_q);
			FRinv_q = Acontra_q;

			P_q[0][0] = sqrt(Acov_q[0][0]);
			P_q[0][1] = Acov_q[0][1] / sqrt(Acov_q[0][0]);
			P_q[1][0] = 0.0;
			P_q[1][1] = sqrt(Acov_q[1][1] - Acov_q[0][1]*Acov_q[0][1]/(Acov_q[0][0]));

			a1_q[0] = f1_n_q[q_index][0]; a1_q[1] = f2_n_q[q_index][0]; a1_q[2] = f3_n_q[q_index][0];
			a2_q[0] = f1_n_q[q_index][1]; a2_q[1] = f2_n_q[q_index][1]; a2_q[2] = f3_n_q[q_index][1];

			n_q = Normalize(Cross(a1_q,a2_q));
			scalednproj_q = (Identity3D - Outer_Product(n_q,n_q))/Cross(a1_q,a2_q).norm();


			acov_q[0][0] = a1_q*a1_q; acov_q[0][1] = a1_q*a2_q; acov_q[1][0] = acov_q[0][1]; acov_q[1][1] = a2_q*a2_q;

			acontra_q = inverse2x2(acov_q);



			for (unsigned int i = 0; i < 2; i++){
				for (unsigned int j = 0; j < 2; j++){
					Tensor<1,3> DDx;
					DDx[0] = Df1_n_q[q_index][i][j]; DDx[1] = Df2_n_q[q_index][i][j]; DDx[2] = Df3_n_q[q_index][i][j];
					bcov_q[i][j] = n_q*DDx;
				}
			}


			E_q = 0.5*P_q*(Acontra_q*acov_q*Acontra_q - Acontra_q)*Transpose(P_q);
			B_q = P_q*acontra_q*bcov_q*acontra_q*Transpose(P_q);


			Material_Vector_InPlane[cell_index].set_Params(Emodv,nu, E_q - Applied_Strain.get_value(cell_index,q_index));
			Material_Vector_Bending[cell_index].set_Params(Emodv, nu, B_q - Applied_Bending.get_value(cell_index, q_index));

			for (const unsigned int i : fe_values.dof_indices()){


				dx_i_q[0] = fe_values[x1].value(i,q_index); dx_i_q[1] = fe_values[x2].value(i,q_index); dx_i_q[2] = fe_values[x3].value(i,q_index);

				df1_i_q = fe_values[f1].value(i,q_index);
				df2_i_q = fe_values[f2].value(i,q_index);
				df3_i_q = fe_values[f3].value(i,q_index);

				Ddf1_i_q = fe_values[f1].gradient(i,q_index);
				Ddf2_i_q = fe_values[f2].gradient(i,q_index);
				Ddf3_i_q = fe_values[f3].gradient(i,q_index);

				dlambda1_i_q = fe_values[lambda1].value(i,q_index);
				dlambda2_i_q = fe_values[lambda2].value(i,q_index);
				dlambda3_i_q = fe_values[lambda3].value(i,q_index);


				da1_i_q[0] = df1_i_q[0]; da1_i_q[1] = df2_i_q[0]; da1_i_q[2] = df3_i_q[0];

				da2_i_q[0] = df1_i_q[1]; da2_i_q[1] = df2_i_q[1]; da2_i_q[2] = df3_i_q[1];


				dacov_i_q[0][0] = 2.0*a1_q*da1_i_q;
				dacov_i_q[0][1] = a1_q*da2_i_q + da1_i_q*a2_q;
				dacov_i_q[1][0] = dacov_i_q[0][1];
				dacov_i_q[1][1] = 2.0*a2_q*da2_i_q;

				dacontra_i_q = -1.0*acontra_q*dacov_i_q*acontra_q;

				dN_i_q = (Cross(da1_i_q,a2_q) + Cross(a1_q,da2_i_q));
				dn_i_q = scalednproj_q*dN_i_q;

				for (unsigned int ii = 0; ii < 2; ii++){
					for (unsigned int ij = 0; ij < 2; ij++){

						DDx[0] = Df1_n_q[q_index][ii][ij]; DDx[1] = Df2_n_q[q_index][ii][ij]; DDx[2] = Df3_n_q[q_index][ii][ij];
						DDdx[0] = Ddf1_i_q[ii][ij]; DDdx[1] = Ddf2_i_q[ii][ij]; DDdx[2] = Ddf3_i_q[ii][ij];
						dbcov_i_q[ii][ij] = n_q*DDdx + dn_i_q*DDx;
					}
				}

				dE_i_q = 0.5*P_q*(Acontra_q*dacov_i_q*Acontra_q)*Transpose(P_q);
				dB_i_q = P_q*(acontra_q * dbcov_i_q * acontra_q + 2.0*Symmetric(acontra_q*bcov_q*dacontra_i_q))*Transpose(P_q);


				for (const unsigned int j : fe_values.dof_indices()){

					dx_j_q[0] = fe_values[x1].value(j,q_index); dx_j_q[1] = fe_values[x2].value(j,q_index); dx_j_q[2] = fe_values[x3].value(j,q_index);

					df1_j_q = fe_values[f1].value(j,q_index);
					df2_j_q = fe_values[f2].value(j,q_index);
					df3_j_q = fe_values[f3].value(j,q_index);

					Ddf1_j_q = fe_values[f1].gradient(j,q_index);
					Ddf2_j_q = fe_values[f2].gradient(j,q_index);
					Ddf3_j_q = fe_values[f3].gradient(j,q_index);

					dlambda1_j_q = fe_values[lambda1].value(j,q_index);
					dlambda2_j_q = fe_values[lambda2].value(j,q_index);
					dlambda3_j_q = fe_values[lambda3].value(j,q_index);


					da1_j_q[0] = df1_j_q[0]; da1_j_q[1] = df2_j_q[0]; da1_j_q[2] = df3_j_q[0];

					da2_j_q[0] = df1_j_q[1]; da2_j_q[1] = df2_j_q[1]; da2_j_q[2] = df3_j_q[1];


					dacov_j_q[0][0] = 2.0*a1_q*da1_j_q;
					dacov_j_q[0][1] = a1_q*da2_j_q + da1_j_q*a2_q;
					dacov_j_q[1][0] = dacov_j_q[0][1];
					dacov_j_q[1][1] = 2.0*a2_q*da2_j_q;

					dacontra_j_q = -1.0*acontra_q*dacov_j_q*acontra_q;

					dN_j_q = (Cross(da1_j_q,a2_q) + Cross(a1_q,da2_j_q));
					dn_j_q = scalednproj_q*dN_j_q;

					for (unsigned int ii = 0; ii < 2; ii++){
						for (unsigned int ij = 0; ij < 2; ij++){

							DDx[0] = Df1_n_q[q_index][ii][ij]; DDx[1] = Df2_n_q[q_index][ii][ij]; DDx[2] = Df3_n_q[q_index][ii][ij];
							DDdx[0] = Ddf1_j_q[ii][ij]; DDdx[1] = Ddf2_j_q[ii][ij]; DDdx[2] = Ddf3_j_q[ii][ij];
							dbcov_j_q[ii][ij] = n_q*DDdx + dn_j_q*DDx;
						}
					}

					dE_j_q = 0.5*P_q*(Acontra_q*dacov_j_q*Acontra_q)*Transpose(P_q);
					dB_j_q = P_q*(acontra_q * dbcov_j_q * acontra_q + 2.0*Symmetric(acontra_q*bcov_q*dacontra_j_q))*Transpose(P_q);


					ddacov_ij_q[0][0] = da1_i_q*da1_j_q;
					ddacov_ij_q[0][1] = da1_i_q*da2_j_q + da1_j_q*da2_i_q;
					ddacov_ij_q[1][0] = ddacov_ij_q[0][1];
					ddacov_ij_q[1][1] = da2_i_q*da2_j_q;

					ddN_ij_q = Cross(da1_j_q,da2_i_q) + Cross(da1_i_q,da2_j_q);
					ddn_ij_q = (-1.0*(n_q*dN_j_q)*dn_i_q  - (Symmetric(Outer_Product(dn_j_q,n_q)))*dN_i_q)/Cross(a1_q,a2_q).norm() + scalednproj_q*ddN_ij_q;



					for (unsigned int ii = 0; ii < 2; ii++){
						for (unsigned int jj = 0; jj < 2; jj++){

							DDxi[0] = Ddf1_i_q[ii][jj]; DDxi[1] = Ddf2_i_q[ii][jj]; DDxi[2] = Ddf3_i_q[ii][jj];

							DDxj[0] = Ddf1_j_q[ii][jj]; DDxj[1] = Ddf2_j_q[ii][jj]; DDxj[2] = Ddf3_j_q[ii][jj];

							DDx[0] = Df1_n_q[q_index][ii][jj]; DDx[1] = Df2_n_q[q_index][ii][jj]; DDx[2] = Df3_n_q[q_index][ii][jj];

							ddbcov_ij_q[ii][jj] = ddn_ij_q*DDx + dn_i_q*DDxj + dn_j_q*DDxi;

						}
					}
					ddE_ij_q = 0.5*P_q*(Acontra_q*ddacov_ij_q*Acontra_q)*Transpose(P_q);
					ddB_ij_q = P_q*(acontra_q*ddbcov_ij_q*acontra_q + 2.0*Symmetric(dacontra_j_q*bcov_q*dacontra_i_q + acontra_q*(dbcov_i_q*dacontra_j_q + dbcov_j_q*dacontra_i_q + bcov_q*ddacov_ij_q)))*Transpose(P_q);



					cell_matrix(i,j) += BilinearProduct(dE_i_q,Material_Vector_InPlane[cell_index].getddQ2ddF(),dE_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),ddE_ij_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += hsc*BilinearProduct(dB_i_q,Material_Vector_Bending[cell_index].getddQ2ddF(),dB_j_q)*fe_values.JxW(q_index);
					cell_matrix(i,j) += hsc*Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),ddB_ij_q)*fe_values.JxW(q_index);

					cell_matrix(i,j) -= (dlambda1_i_q*(fe_values[x1].gradient(j,q_index) - df1_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) -= (dlambda1_j_q*(fe_values[x1].gradient(i,q_index) - df1_i_q))*fe_values.JxW(q_index);

					cell_matrix(i,j) -= (dlambda2_i_q*(fe_values[x2].gradient(j,q_index) - df2_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) -= (dlambda2_j_q*(fe_values[x2].gradient(i,q_index) - df2_i_q))*fe_values.JxW(q_index);

					cell_matrix(i,j) -= (dlambda3_i_q*(fe_values[x3].gradient(j,q_index) - df3_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) -= (dlambda3_j_q*(fe_values[x3].gradient(i,q_index) - df3_i_q))*fe_values.JxW(q_index);

					cell_matrix(i,j) += mu*((fe_values[x1].gradient(i,q_index) - df1_i_q)*(fe_values[x1].gradient(j,q_index) - df1_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += mu*((fe_values[x2].gradient(i,q_index) - df2_i_q)*(fe_values[x2].gradient(j,q_index) - df2_j_q))*fe_values.JxW(q_index);
					cell_matrix(i,j) += mu*((fe_values[x3].gradient(i,q_index) - df3_i_q)*(fe_values[x3].gradient(j,q_index) - df3_j_q))*fe_values.JxW(q_index);

					cell_matrix(i,j) += xdrag*(fe_values[x1].value(i,q_index)*fe_values[x1].value(j,q_index) + fe_values[x2].value(i,q_index)*fe_values[x2].value(j,q_index) + fe_values[x3].value(i,q_index)*fe_values[x3].value(j,q_index))*fe_values.JxW(q_index);
				}

				cell_rhs(i) += Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),dE_i_q)*fe_values.JxW(q_index);
				cell_rhs(i) += hsc*Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),dB_i_q)*fe_values.JxW(q_index);

				cell_rhs(i) -= lambda1_n_q[q_index]*((fe_values[x1].gradient(i,q_index) - df1_i_q))*fe_values.JxW(q_index);
				cell_rhs(i) -= dlambda1_i_q*(Dx1_n_q[q_index] - f1_n_q[q_index])*fe_values.JxW(q_index);

				cell_rhs(i) -= lambda2_n_q[q_index]*((fe_values[x2].gradient(i,q_index) - df2_i_q))*fe_values.JxW(q_index);
				cell_rhs(i) -= dlambda2_i_q*(Dx2_n_q[q_index] - f2_n_q[q_index])*fe_values.JxW(q_index);

				cell_rhs(i) -= lambda3_n_q[q_index]*((fe_values[x3].gradient(i,q_index) - df3_i_q))*fe_values.JxW(q_index);
				cell_rhs(i) -= dlambda3_i_q*(Dx3_n_q[q_index] - f3_n_q[q_index])*fe_values.JxW(q_index);

				cell_rhs(i) += mu*(Dx1_n_q[q_index] - f1_n_q[q_index])*((fe_values[x1].gradient(i,q_index) - df1_i_q))*fe_values.JxW(q_index);
				cell_rhs(i) += mu*(Dx2_n_q[q_index] - f2_n_q[q_index])*((fe_values[x2].gradient(i,q_index) - df2_i_q))*fe_values.JxW(q_index);
				cell_rhs(i) += mu*(Dx3_n_q[q_index] - f3_n_q[q_index])*((fe_values[x3].gradient(i,q_index) - df3_i_q))*fe_values.JxW(q_index);


				cell_rhs(i) += xdrag*((x1_n_q[q_index] - x1old_n_q[q_index])*fe_values[x1].value(i,q_index)
						+(x2_n_q[q_index] - x2old_n_q[q_index])*fe_values[x2].value(i,q_index)
						+(x3_n_q[q_index] - x3old_n_q[q_index])*fe_values[x3].value(i,q_index))*fe_values.JxW(q_index);
			}

		}
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

void ElasticProblem::initialize_configuration(){
	const int ndofs = dof_handler.n_dofs();

	//////////////////////////////////
	std::vector<bool> x1_components = create_bool_vector(15,0);
	ComponentMask x1_mask(x1_components);
	std::vector<bool> is_x1_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, x1_mask, is_x1_comp);
	//////////////////////////////////
	std::vector<bool> x2_components = create_bool_vector(15,1);
	ComponentMask x2_mask(x2_components);
	std::vector<bool> is_x2_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, x2_mask, is_x2_comp);
	//////////////////////////////////
	std::vector<bool> x3_components = create_bool_vector(15,2);
	ComponentMask x3_mask(x3_components);
	std::vector<bool> is_x3_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, x3_mask, is_x3_comp);
	//////////////////////////////////
	std::vector<bool> f11_components = create_bool_vector(15,3);
	ComponentMask f11_mask(f11_components);
	std::vector<bool> is_f11_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f11_mask, is_f11_comp);
	//////////////////////////////////
	std::vector<bool> f12_components = create_bool_vector(15,4);
	ComponentMask f12_mask(f12_components);
	std::vector<bool> is_f12_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f12_mask, is_f12_comp);
	//////////////////////////////////
	std::vector<bool> f21_components = create_bool_vector(15,5);
	ComponentMask f21_mask(f21_components);
	std::vector<bool> is_f21_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f21_mask, is_f21_comp);
	//////////////////////////////////
	std::vector<bool> f22_components = create_bool_vector(15,6);
	ComponentMask f22_mask(f22_components);
	std::vector<bool> is_f22_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f22_mask, is_f22_comp);
	//////////////////////////////////
	std::vector<bool> f31_components = create_bool_vector(15,7);
	ComponentMask f31_mask(f31_components);
	std::vector<bool> is_f31_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f31_mask, is_f31_comp);
	//////////////////////////////////
	std::vector<bool> f32_components = create_bool_vector(15,8);
	ComponentMask f32_mask(f32_components);
	std::vector<bool> is_f32_comp(ndofs, false);
	DoFTools::extract_dofs(dof_handler, f32_mask, is_f32_comp);

	std::vector<Point<DIM>> support_points(ndofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);


	for (unsigned int i = 0; i < ndofs; i++){
		if (is_x1_comp[i]){
			solution[i] = support_points[i][0];
		} else if (is_x2_comp[i]){
			solution[i] = support_points[i][1];

		} else if (is_x3_comp[i]){
			solution[i] = 0.0;

		} else if (is_f11_comp[i]){
			solution[i] = 1.0;
		} else if (is_f12_comp[i]){
			solution[i] = 0.0;
		} else if (is_f21_comp[i]){
			solution[i] = 0.0;
		} else if (is_f22_comp[i]){
			solution[i] = 1.0;
		}
	}

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

Tensor<1,3> ElasticProblem::Normalize(const Tensor<1,3> & v){
	return v/v.norm();
}

Tensor<1,3> ElasticProblem::Cross(const Tensor<1,3> & v1, const Tensor<1,3> & v2){
	Tensor<1,3> vout;

	vout[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vout[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vout[2] = v1[0]*v2[1] - v1[1]*v2[0];

	return vout;

}

Tensor<2,2> ElasticProblem::Transpose(const Tensor<2,2> & A){
	Tensor<2,2> Aout;
	for (unsigned int i = 0; i < 2; i++){
		for (unsigned int j = 0; j < 2; j++){
			Aout[i][j] = A[j][i];
		}
	}
	return Aout;
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

template<int d>
Tensor<2,d> ElasticProblem::Outer_Product(const Tensor<1,d> & v1, const Tensor<1,d>& v2){
	Tensor<2,d> Tout;

	for (unsigned int i = 0; i < d; i++){
		for (unsigned int j = 0; j < d; j++){
			Tout[i][j] = v1[i]*v2[j];
		}
	}
	return Tout;

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
		std::cout << "Assembling System" << std::endl;
		assemble_system();
		std::cout << "Solving Linear System" << std::endl;

		solve();


		stepsize = sqrt(linearsolve.norm_sqr())/linearsolve.size();
		//residual = step_direction;

		std::cout << "Iteration: " << cntr << std::endl;
		std::cout << "Step Size: " << stepsize<< std::endl;
	}
	std::cout << "Number of Steps: " << cntr << std::endl;
}



void ElasticProblem::make_grid()
{


	GridGenerator::hyper_cube(triangulation, -1.0, 1.0);
	triangulation.refine_global(refinelevel);
	std::cout << "   Number of active cells: " << triangulation.n_active_cells()
																																																																																																																																																																	<< std::endl << "   Total number of cells: "
																																																																																																																																																																	<< triangulation.n_cells() << std::endl;
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



/*
void ElasticProblem::output_results(){
	std::vector<std::string> solution_names(3,"x");
	solution_names.push_back("f");
	solution_names.push_back("f");
	solution_names.push_back("f");
	solution_names.push_back("f");
	solution_names.push_back("f");
	solution_names.push_back("f");
	solution_names.push_back("lambda");
	solution_names.push_back("lambda");
	solution_names.push_back("lambda");
	solution_names.push_back("lambda");
	solution_names.push_back("lambda");
	solution_names.push_back("lambda");

	std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
		      3, DataComponentInterpretation::component_is_part_of_vector);

	for (unsigned int i = 0; i < 6; i++){
		data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
	}
	for (unsigned int i = 0; i < 6; i++){
			data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
		}

	 DataOut<2> data_out;
	  data_out.attach_dof_handler(dof_handler);
	  data_out.add_data_vector(solution,
	                           solution_names,
	                           DataOut<2>::type_dof_data,
	                           data_component_interpretation);
	  data_out.build_patches();
	  std::ofstream output("output.vtk");
	  data_out.write_vtk(output);


}
 */

void ElasticProblem::output_results(const unsigned int cycle) const
{
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);

	std::vector<std::string> solution_names;

	solution_names.emplace_back("x_displacement");
	solution_names.emplace_back("y_displacement");
	solution_names.emplace_back("z_displacement");
	solution_names.emplace_back("w_1_gradient");
	solution_names.emplace_back("w_2_gradient");
	solution_names.emplace_back("w_3_gradient");
	solution_names.emplace_back("w_4_gradient");
	solution_names.emplace_back("w_5_gradient");
	solution_names.emplace_back("w_6_gradient");
	solution_names.emplace_back("l_1");
	solution_names.emplace_back("l_2");
	solution_names.emplace_back("l_3");
	solution_names.emplace_back("l_4");
	solution_names.emplace_back("l_5");
	solution_names.emplace_back("l_6");


	// After setting up the names for the different components of the
	// solution vector, we can add the solution vector to the list of
	// data vectors scheduled for output. Note that the following
	// function takes a vector of strings as second argument, whereas
	// the one which we have used in all previous examples accepted a
	// string there. (In fact, the function we had used before would
	// convert the single string into a vector with only one element
	// and forwards that to the other function.)
	data_out.add_data_vector(solution, solution_names);
	data_out.build_patches();

	std::ofstream output("solutions/solution-" + std::to_string(cycle) + ".vtk");
	data_out.write_vtk(output);
}


void ElasticProblem::setup_system()
{


	dof_handler.distribute_dofs(fe);
	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
	linearsolve.reinit(dof_handler.n_dofs());
	prev_solution.reinit(dof_handler.n_dofs());
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


	QGauss<DIM> quadrature_formula(fe.degree + quadegadd);




	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);


	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	Material_Vector_InPlane.resize(triangulation.n_active_cells());
	Material_Vector_Bending.resize(triangulation.n_active_cells());
	Applied_Strain.resize(triangulation.n_active_cells(),n_q_points);
	Applied_Bending.resize(triangulation.n_active_cells(),n_q_points);




}

void ElasticProblem::run()
{
	std::cout << "Solving problem in " << DIM << " space dimensions."
			<< std::endl;
	make_grid();
	setup_system();
	initialize_configuration();
	solve_path();


}


}

#endif // ELASTIC_PROBLEM_CC_

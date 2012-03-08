/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//Defines the post-processing steps in the class "Ctr_pd_ml" for a PHASOR-DOMAIN near-field-to-far-field transformer in a general N-layered medium

#include "headers.h"

#include "Ctr_pd_ml.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//definition of Cpointsources needed
#include "pointsources/Cpointsources.h"

extern double dx,dt;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif

extern int number_of_layers;

extern Array<int,1> LayerLowerZIndices;
extern Array<int,1> LayerThicknesses;

extern int rank;

namespace{
	int i,j,k,pln;
};

void Ctr_pd_ml::ConstructFarField()
{
	//variables needed for displaying the progress of post-processing (in % completed)
	int d_counter = 0;	//counter is needed because angles may be done out of order

	int l,d1,d2;
	for (d2=0; d2<D2; d2++)
	{
		for (d1=0; d1<D1; d1++)
		{
			for (l=0; l<L; l++)
			{
				if (observation_angle_exists(l,d1,d2))
				{//compute the far-field arrays only if there is a valid observation angle for the direction parameter (dir1,dir2)
					if (abs(CosT(l,d1,d2))>abs(threshold))		// execute the loop only if angle is in the LOWER or UPPER half space
					{
						//Calculate transmission-line parameters for the given l (frequency) and d1,d2 (observation-angle parameters)
						Calculate_TL_GreenFunctions(l,d1,d2);

						ConstructPotential_A(l,d1,d2);	//compute A_theta,A_phi
						ConstructPotential_F(l,d1,d2);	//compute F_theta,F_phi
						//construct the electric far field from the electric and magnetic potentials
						E_theta(l,d1,d2) = dt*pow2(dx)*1/(4*M_PI)*(ii)*(ww(l)*(-mu_r_o*mu_0*A_theta(l,d1,d2)) - kk_o_g*F_phi(l,d1,d2));
						E_phi(l,d1,d2) = dt*pow2(dx)*1/(4*M_PI)*(ii)*(ww(l)*(-mu_r_o*mu_0*A_phi(l,d1,d2)) + kk_o_g*F_theta(l,d1,d2));
						//apply the correction due to the distance of the actual far-field origin from the uppermost interface
						E_theta(l,d1,d2) *= exp(-ii*kk_o_g*(dx*(FarFieldOriginZ-TemporaryFarFieldOriginZ))*cosT);
						E_phi(l,d1,d2) *= exp(-ii*kk_o_g*(dx*(FarFieldOriginZ-TemporaryFarFieldOriginZ))*cosT);

					}
					d_counter++;
					//Wait until all nodes are here, for the far-field processing progress counter to be accurate
					//(a more efficient way can be found later)
#ifndef MPI_DISABLE
					MPI_Barrier(MPI_CartSubComm);
#endif
                }
                else
                {
                    d_counter++;	//count the non-existent angles here
                }
                if (rank==0)
                {
                    cout << "Far-field post-processing : " << ceil(100.0*d_counter/(L*D1*D2)) << "% completed.\r" ;
                }
            }
		}
	}
	if (rank==0)
	{
		cout << endl;
	}
	//*****************************************************************************************//
	//*****************************************************************************************//

	// Gather the far-field arrays and write into file
	GatherAndWriteFarField();
}

void Ctr_pd_ml::ConstructPotential_A(const int& l, const int& d1, const int& d2)
{
	//Back face of the virtual surface
	// J_y_back is stored at the nodes containing H-field components in the outer H-shell
	if (Back_Hz_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hz_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hz_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_y_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_y_plus(pln));
					A_phi(l,d1,d2) += J_y_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_y_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_y_plus(pln));
				}
			}
		}
	}
	// J_z_back is stored at the nodes containing H-field components in the outer H-shell
	if (Back_Hy_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hy_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hy_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_z_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_z_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_z_plus(pln));
				}
			}
		}
	}

	//Front face of the virtual surface
	// J_y_front is stored at the nodes containing H-field components in the outer H-shell
	if (Front_Hz_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hz_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hz_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_y_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_y_plus(pln));
					A_phi(l,d1,d2) += J_y_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_y_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_y_plus(pln));
				}
			}
		}
	}
	// J_z_front is stored at the nodes containing H-field components in the outer H-shell
	if (Front_Hy_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hy_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hy_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_z_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_z_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_z_plus(pln));
				}
			}
		}
	}

	//Left face of the virtual surface
	// J_x_left is stored at the nodes containing H-field components in the outer H-shell
	if (Left_Hz_PassesThroughNode)
	{
		j=Left_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hz_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hz_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
				{
					A_theta(l,d1,d2) += J_x_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_x_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_x_plus(pln));
					A_phi(l,d1,d2) += J_x_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_x_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_x_plus(pln));
				}
			}
		}
	}
	// J_z_left is stored at the nodes containing H-field components in the outer H-shell
	if (Left_Hx_PassesThroughNode)
	{
		j=Left_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hx_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hx_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
				{
					A_theta(l,d1,d2) += J_z_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_z_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_z_plus(pln));
				}
			}
		}
	}

	//Right face of the virtual surface
	// J_x_right is stored at the nodes containing H-field components in the outer H-shell
	if (Right_Hz_PassesThroughNode)
	{
		j=Right_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hz_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hz_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
				{
					A_theta(l,d1,d2) += J_x_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_x_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_x_plus(pln));
					A_phi(l,d1,d2) += J_x_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_x_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_x_plus(pln));
				}
			}
		}
	}
	// J_z_right is stored at the nodes containing H-field components in the outer H-shell
	if (Right_Hx_PassesThroughNode)
	{
		j=Right_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Hx_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Hx_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
				{
					A_theta(l,d1,d2) += J_z_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_z_minus(pln)
							 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_z_plus(pln));
				}
			}
		}
	}

	//Bottom face of the virtual surface
	// J_x_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hy_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
				{
					for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
					{
						A_theta(l,d1,d2) += J_x_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_x_plus(pln));
						A_phi(l,d1,d2) += J_x_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_x_plus(pln));
					}
				}
			}
		}
	}
	// J_y_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hx_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
				{
					for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
					{
						A_theta(l,d1,d2) += J_y_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_y_plus(pln));
						A_phi(l,d1,d2) += J_y_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_y_plus(pln));
					}
				}
			}
		}
	}

	//Upper face of the virtual surface
	// J_x_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hy_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
				{
					for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
					{
						A_theta(l,d1,d2) += J_x_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_x_plus(pln));
						A_phi(l,d1,d2) += J_x_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_x_plus(pln));
					}
				}
			}
		}
	}
	// J_y_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hx_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
				{
					for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
					{
						A_theta(l,d1,d2) += J_y_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_theta_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_theta_J_y_plus(pln));
						A_phi(l,d1,d2) += J_y_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_A_phi_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_A_phi_J_y_plus(pln));
					}
				}
			}
		}
	}
}

void Ctr_pd_ml::ConstructPotential_F(const int& l, const int& d1, const int& d2)
{
	//Back face of the virtual surface
	// M_y_back is stored at the nodes containing E-field components in the outer E-shell
	if (Back_Ez_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ez_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ez_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_y_plus(pln));
					F_phi(l,d1,d2) += M_y_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_y_plus(pln));
				}
			}
		}
	}
	// M_z_back is stored at the nodes containing E-field components in the outer E-shell
	if (Back_Ey_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ey_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ey_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_z_back(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_z_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_z_plus(pln));
				}
			}
		}
	}

	//Front face of the virtual surface
	// M_y_front is stored at the nodes containing E-field components in the outer E-shell
	if (Front_Ez_PassesThroughNode)
	{
		i=Front_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ez_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ez_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_y_plus(pln));
					F_phi(l,d1,d2) += M_y_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_y_plus(pln));
				}
			}
		}
	}
	// M_z_front is stored at the nodes containing E-field components in the outer E-shell
	if (Front_Ey_PassesThroughNode)
	{
		i=Front_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ey_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ey_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_z_front(l,j,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_z_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_z_plus(pln));
				}
			}
		}
	}

	//Left face of the virtual surface
	// M_x_left is stored at the nodes containing E-field components in the outer E-shell
	if (Left_Ez_PassesThroughNode)
	{
		j=Left_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ez_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ez_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
				{
					F_theta(l,d1,d2) += M_x_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_x_plus(pln));
					F_phi(l,d1,d2) += M_x_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_x_plus(pln));
				}
			}
		}
	}
	// M_z_left is stored at the nodes containing E-field components in the outer E-shell
	if (Left_Ex_PassesThroughNode)
	{
		j=Left_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ex_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ex_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
				{
					F_theta(l,d1,d2) += M_z_left(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_z_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_z_plus(pln));
				}
			}
		}
	}

	//Right face of the virtual surface
	// M_x_right is stored at the nodes containing E-field components in the outer E-shell
	if (Right_Ez_PassesThroughNode)
	{
		j=Right_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ez_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ez_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
				{
					F_theta(l,d1,d2) += M_x_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_x_plus(pln));
					F_phi(l,d1,d2) += M_x_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_x_plus(pln));
				}
			}
		}
	}
	// M_z_right is stored at the nodes containing E-field components in the outer E-shell
	if (Right_Ex_PassesThroughNode)
	{
		j=Right_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			for (k=max(LayerLowerZIndices(pln),Ex_lowerlimit_in_node);k<=min(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1,Ex_upperlimit_in_node);k++)
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln);
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln));
				for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
				{
					F_theta(l,d1,d2) += M_z_right(l,i,k)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_z_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_z_plus(pln));
				}
			}
		}
	}

	//Bottom face of the virtual surface
	// M_x_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ey_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
				{
					for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
					{
						F_theta(l,d1,d2) += M_x_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_x_plus(pln));
						F_phi(l,d1,d2) += M_x_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_x_plus(pln));
					}
				}
			}
		}
	}
	// M_y_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ex_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
				{
					for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
					{
						F_theta(l,d1,d2) += M_y_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_y_plus(pln));
						F_phi(l,d1,d2) += M_y_lower(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_y_plus(pln));
					}
				}
			}
		}
	}

	//Upper face of the virtual surface
	// M_x_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ey_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
				{
					for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
					{
						F_theta(l,d1,d2) += M_x_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_x_plus(pln));
						F_phi(l,d1,d2) += M_x_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_x_plus(pln));
					}
				}
			}
		}
	}
	// M_y_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ex_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (pln=0; pln<number_of_layers; pln++)
		{
			if ((k>=LayerLowerZIndices(pln))&&(k<=(LayerLowerZIndices(pln)+LayerThicknesses(pln)-1)))
			{
				z_pr_minus_z_n = k - LayerLowerZIndices(pln) + 0.5;
				z_pr_minus_z_n_1 = k - (LayerLowerZIndices(pln)+LayerThicknesses(pln)) + 0.5;
				for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
				{
					for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
					{
						F_theta(l,d1,d2) += M_y_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_theta_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_theta_M_y_plus(pln));
						F_phi(l,d1,d2) += M_y_upper(l,i,j)
								*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
								*(exp(ii*kt_dx(pln)*z_pr_minus_z_n_1)*Green_F_phi_M_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*z_pr_minus_z_n)*Green_F_phi_M_y_plus(pln));
					}
				}
			}
		}
	}
}

complex<double> Ctr_pd_ml::TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2)
//Theoretical theta-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	/** FIXME: Could be made more efficient **/
	Calculate_TL_GreenFunctions(l,d1,d2); //recalculate the TL parameters
	/** FIXME: Could be made more efficient **/

	total_theta_component = 0;
	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		int srcindex;
		for (srcindex=0; srcindex<Data.PointSourcesPtr->NumberOfElectricDipoles(); srcindex++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(srcindex)->Orientation;
			x = Data.PointSourcesPtr->ElectricDipole(srcindex)->x;
			y = Data.PointSourcesPtr->ElectricDipole(srcindex)->y;
			z = Data.PointSourcesPtr->ElectricDipole(srcindex)->z;
			if (Orientation=="x_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX+0.5;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ;
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY+0.5;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ;
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			theta_component = Data.PointSourcesPtr->ElectricDipole(srcindex)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			for (pln=0; pln<number_of_layers; pln++)
			{
				if ((z>=LayerLowerZIndices(pln))&&(z<=LayerLowerZIndices(pln)+LayerThicknesses(pln)-1))
				{
					if (Orientation=="x_directed")
					{
						SourceZ_minus_z_n = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-LayerLowerZIndices(pln);
						SourceZ_minus_z_n_1 = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-(LayerLowerZIndices(pln)+LayerThicknesses(pln));
						theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
								*(exp(ii*kt_dx(pln)*SourceZ_minus_z_n_1)*Green_A_theta_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*SourceZ_minus_z_n)*Green_A_theta_J_x_plus(pln));
					}
					else if (Orientation=="y_directed")
					{
						SourceZ_minus_z_n = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-LayerLowerZIndices(pln);
						SourceZ_minus_z_n_1 = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-(LayerLowerZIndices(pln)+LayerThicknesses(pln));
						theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
								*(exp(ii*kt_dx(pln)*SourceZ_minus_z_n_1)*Green_A_theta_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*SourceZ_minus_z_n)*Green_A_theta_J_y_plus(pln));
					}
					else if (Orientation=="z_directed")
					{
						SourceZ_minus_z_n = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-LayerLowerZIndices(pln)+0.5;
						SourceZ_minus_z_n_1 = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-(LayerLowerZIndices(pln)+LayerThicknesses(pln))+0.5;
						theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
								*(exp(ii*kt_dx(pln)*SourceZ_minus_z_n_1)*Green_A_theta_J_z_minus(pln)
								 +exp(-ii*kt_dx(pln)*SourceZ_minus_z_n)*Green_A_theta_J_z_plus(pln));
					}
				}
			}

			//apply the final prefactors
			theta_component *= 1/(4*M_PI)*(ii)*ww(l)*(-mu_r_o*mu_0);

			theta_component *= exp(-ii*kk_o_g*(dx*(FarFieldOriginZ-TemporaryFarFieldOriginZ))*cosT);

			//add contribution of this dipole
			total_theta_component += theta_component;
		}
	}

	return total_theta_component;
// 	return 0;
}

complex<double> Ctr_pd_ml::TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2)
//Theoretical phi-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	/** FIXME: Could be made more efficient **/
	Calculate_TL_GreenFunctions(l,d1,d2); //recalculate the TL parameters
	/** FIXME: Could be made more efficient **/

	total_phi_component = 0;
	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		int srcindex;
		for (srcindex=0; srcindex<Data.PointSourcesPtr->NumberOfElectricDipoles(); srcindex++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(srcindex)->Orientation;
			x = Data.PointSourcesPtr->ElectricDipole(srcindex)->x;
			y = Data.PointSourcesPtr->ElectricDipole(srcindex)->y;
			z = Data.PointSourcesPtr->ElectricDipole(srcindex)->z;
			if (Orientation=="x_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX+0.5;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ;
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY+0.5;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ;
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(srcindex)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(srcindex)->y-FarFieldOriginY;
//				SourceZ = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-TemporaryFarFieldOriginZ+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			phi_component = Data.PointSourcesPtr->ElectricDipole(srcindex)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			for (pln=0; pln<number_of_layers; pln++)
			{
				if ((z>=LayerLowerZIndices(pln))&&(z<=LayerLowerZIndices(pln)+LayerThicknesses(pln)-1))
				{
					if (Orientation=="x_directed")
					{
						SourceZ_minus_z_n = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-LayerLowerZIndices(pln);
						SourceZ_minus_z_n_1 = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-(LayerLowerZIndices(pln)+LayerThicknesses(pln));
						phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
								*(exp(ii*kt_dx(pln)*SourceZ_minus_z_n_1)*Green_A_phi_J_x_minus(pln)
								 +exp(-ii*kt_dx(pln)*SourceZ_minus_z_n)*Green_A_phi_J_x_plus(pln));
					}
					else if (Orientation=="y_directed")
					{
						SourceZ_minus_z_n = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-LayerLowerZIndices(pln);
						SourceZ_minus_z_n_1 = Data.PointSourcesPtr->ElectricDipole(srcindex)->z-(LayerLowerZIndices(pln)+LayerThicknesses(pln));
						phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
								*(exp(ii*kt_dx(pln)*SourceZ_minus_z_n_1)*Green_A_phi_J_y_minus(pln)
								 +exp(-ii*kt_dx(pln)*SourceZ_minus_z_n)*Green_A_phi_J_y_plus(pln));
					}
					else if (Orientation=="z_directed")
					{
						phi_component *= 0;	//z-directed dipole does not contribute to the phi-component of the radiated E-field
					}
				}
			}

			//apply the final prefactors
			phi_component *= 1/(4*M_PI)*(ii)*ww(l)*(-mu_r_o*mu_0);

			//apply the correction due to the distance of the actual far-field origin from the interface
			phi_component *= exp(-ii*kk_o_g*(dx*(FarFieldOriginZ-TemporaryFarFieldOriginZ))*cosT);

			//add contribution of this dipole
			total_phi_component += phi_component;
		}
	}

	return total_phi_component;
// 	return 0;
}

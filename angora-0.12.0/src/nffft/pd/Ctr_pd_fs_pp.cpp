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

//Defines the post-processing steps in the PHASOR-DOMAIN near-field-to-far-field transformer object "Ctr_pd_fs" in free space

#include "headers.h"

#include "Ctr_pd_fs.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//definition of Cpointsources needed
#include "pointsources/Cpointsources.h"

extern double dx,dt;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif

extern int rank;

extern int i,j,k;


void Ctr_pd_fs::ConstructFarField()
{
	//counter variable for displaying the progress of post-processing (in % completed)
	int d_counter = 0;

	for (int d2=0; d2<D2; d2++)
	{
		for (int d1=0; d1<D1; d1++)
		{
			for (int l=0; l<L; l++)
			{
				if (observation_angle_exists(l,d1,d2))
				{//compute the far-field arrays only if there is a valid observation angle for the direction parameter (dir1,dir2)
					ConstructPotential_A(l,d1,d2);	//compute the cartesian components of the electric potential A_x,A_y,A_z
					ConstructPotential_F(l,d1,d2);	//compute the cartesian components of the magnetic potential F_x,F_y,F_z
					//construct the theta and phi components from the cartesian components
					A_theta(l,d1,d2) = CosTCosP(l,d1,d2)*A_x(l,d1,d2) + CosTSinP(l,d1,d2)*A_y(l,d1,d2) - SinT(l,d1,d2)*A_z(l,d1,d2);
					A_phi(l,d1,d2) = -SinP(l,d1,d2)*A_x(l,d1,d2) + CosP(l,d1,d2)*A_y(l,d1,d2);
					F_theta(l,d1,d2) = CosTCosP(l,d1,d2)*F_x(l,d1,d2) + CosTSinP(l,d1,d2)*F_y(l,d1,d2) - SinT(l,d1,d2)*F_z(l,d1,d2);
					F_phi(l,d1,d2) = -SinP(l,d1,d2)*F_x(l,d1,d2) + CosP(l,d1,d2)*F_y(l,d1,d2);
					//construct the electric far field from the electric and magnetic potentials
					E_theta(l,d1,d2) = dt*pow2(dx)*K*(ii)*kk_g(l,d1,d2)*(-Z_0*A_theta(l,d1,d2) - F_phi(l,d1,d2));
					E_phi(l,d1,d2) = dt*pow2(dx)*K*(ii)*kk_g(l,d1,d2)*(-Z_0*A_phi(l,d1,d2) + F_theta(l,d1,d2));

					d_counter++;
                }
                else
                {// all arrays are left 0
                    d_counter++;	//count the non-existent angles here
                }
                if (rank==0)
                {
                    cout << "Far-field post-processing : " << ceil(100.0*d_counter/(L*D1*D2)) << "% completed.\r" ;
                }
			}
		}

		//Wait until all nodes are here, for the far-field processing progress counter to be accurate
		//(a more efficient way can be found later)
#ifndef MPI_DISABLE
		MPI_Barrier(MPI_CartSubComm);
#endif
	}
	if (rank==0)
	{
		cout << endl;
	}

	// Gather the far-field arrays and write into file
	GatherAndWriteFarField();
}

void Ctr_pd_fs::ConstructPotential_A(const int& l, const int& d1, const int& d2)
{
	//constructs the cartesian components of the magnetic potential A_x,A_y,A_z from the magnetic equivalent currents J_x,J_y,J_z

    //define the temporary variables for fast post-processing
    sinTcosP = SinTCosP(l,d1,d2);
    sinTsinP = SinTSinP(l,d1,d2);
    cosT = CosT(l,d1,d2);
    k_dx = kk_g(l,d1,d2)*dx;

	//Back face of the virtual surface
	// J_y_back is stored at the nodes containing E-field components in the outer H-shell
	if (Back_Hz_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				A_y(l,d1,d2) += J_y_back(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_z_back is stored at the nodes containing H-field components in the outer H-shell
	if (Back_Hy_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=Hy_lowerlimit_in_node;k<=Hy_upperlimit_in_node;k++)
			{
				A_z(l,d1,d2) += J_z_back(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}

	//Front face of the virtual surface
	// J_y_front is stored at the nodes containing E-field components in the outer H-shell
	if (Front_Hz_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				A_y(l,d1,d2) += J_y_front(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_z_front is stored at the nodes containing H-field components in the outer H-shell
	if (Front_Hy_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=Hy_lowerlimit_in_node;k<=Hy_upperlimit_in_node;k++)
			{
				A_z(l,d1,d2) += J_z_front(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}

	//Left face of the virtual surface
	// J_x_left is stored at the nodes containing H-field components in the outer H-shell
	if (Left_Hz_PassesThroughNode)
	{
		j=Left_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				A_x(l,d1,d2) += J_x_left(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_z_left is stored at the nodes containing H-field components in the outer H-shell
	if (Left_Hx_PassesThroughNode)
	{
		j=Left_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=Hx_lowerlimit_in_node;k<=Hx_upperlimit_in_node;k++)
			{
				A_z(l,d1,d2) += J_z_left(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}

	//Right face of the virtual surface
	// J_x_right is stored at the nodes containing H-field components in the outer H-shell
	if (Right_Hz_PassesThroughNode)
	{
		j=Right_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				A_x(l,d1,d2) += J_x_right(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_z_right is stored at the nodes containing H-field components in the outer H-shell
	if (Right_Hx_PassesThroughNode)
	{
		j=Right_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=Hx_lowerlimit_in_node;k<=Hx_upperlimit_in_node;k++)
			{
				A_z(l,d1,d2) += J_z_right(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}

	//Bottom face of the virtual surface
	// J_x_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hy_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
		{
			for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
			{
				A_x(l,d1,d2) += J_x_lower(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_y_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hx_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
			{
				A_y(l,d1,d2) += J_y_lower(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}

	//Upper face of the virtual surface
	// J_x_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hy_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
		{
			for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
			{
				A_x(l,d1,d2) += J_x_upper(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
	// J_y_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hx_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
			{
				A_y(l,d1,d2) += J_y_upper(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}
}

void Ctr_pd_fs::ConstructPotential_F(const int& l, const int& d1, const int& d2)
{
	//constructs the cartesian components of the magnetic potential F_x,F_y,F_z from the magnetic equivalent currents M_x,M_y,M_z

    //define the temporary variables for fast post-processing
    sinTcosP = SinTCosP(l,d1,d2);
    sinTsinP = SinTSinP(l,d1,d2);
    cosT = CosT(l,d1,d2);
    k_dx = kk_g(l,d1,d2)*dx;

	//Back face of the virtual surface
	// M_y_back is stored at the nodes containing E-field components in the outer E-shell
	if (Back_Ez_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				F_y(l,d1,d2) += M_y_back(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_z_back is stored at the nodes containing H-field components in the outer E-shell
	if (Back_Ey_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=Ey_lowerlimit_in_node;k<=Ey_upperlimit_in_node;k++)
			{
				F_z(l,d1,d2) += M_z_back(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}

	//Front face of the virtual surface
	// M_y_front is stored at the nodes containing E-field components in the outer E-shell
	if (Front_Ez_PassesThroughNode)
	{
		i=Front_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				F_y(l,d1,d2) += M_y_front(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_z_front is stored at the nodes containing H-field components in the outer E-shell
	if (Front_Ey_PassesThroughNode)
	{
		i=Front_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=Ey_lowerlimit_in_node;k<=Ey_upperlimit_in_node;k++)
			{
				F_z(l,d1,d2) += M_z_front(l,j,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}

	//Left face of the virtual surface
	// M_x_left is stored at the nodes containing E-field components in the outer E-shell
	if (Left_Ez_PassesThroughNode)
	{
		j=Left_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				F_x(l,d1,d2) += M_x_left(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_z_left is stored at the nodes containing E-field components in the outer E-shell
	if (Left_Ex_PassesThroughNode)
	{
		j=Left_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=Ex_lowerlimit_in_node;k<=Ex_upperlimit_in_node;k++)
			{
				F_z(l,d1,d2) += M_z_left(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}

	//Right face of the virtual surface
	// M_x_right is stored at the nodes containing E-field components in the outer E-shell
	if (Right_Ez_PassesThroughNode)
	{
		j=Right_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				F_x(l,d1,d2) += M_x_right(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_z_right is stored at the nodes containing E-field components in the outer E-shell
	if (Right_Ex_PassesThroughNode)
	{
		j=Right_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=Ex_lowerlimit_in_node;k<=Ex_upperlimit_in_node;k++)
			{
				F_z(l,d1,d2) += M_z_right(l,i,k)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ)*cosT));
			}
		}
	}

	//Bottom face of the virtual surface
	// M_x_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ey_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
		{
			for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
			{
				F_x(l,d1,d2) += M_x_lower(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_y_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ex_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
			{
				F_y(l,d1,d2) += M_y_lower(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}

	//Upper face of the virtual surface
	// M_x_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ey_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
		{
			for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
			{
				F_x(l,d1,d2) += M_x_upper(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
	// M_y_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ex_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
			{
				F_y(l,d1,d2) += M_y_upper(l,i,j)*exp(ii*k_dx*
						((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP+(k-FarFieldOriginZ+0.5)*cosT));
			}
		}
	}
}

complex<double> Ctr_pd_fs::TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2)
//Theoretical theta-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	total_theta_component = 0;
	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		for (int i=0; i<Data.PointSourcesPtr->NumberOfElectricDipoles(); i++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(i)->Orientation;
			x = Data.PointSourcesPtr->ElectricDipole(i)->x;
			y = Data.PointSourcesPtr->ElectricDipole(i)->y;
			z = Data.PointSourcesPtr->ElectricDipole(i)->z;
			if (Orientation=="x_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX+0.5;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ;
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY+0.5;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ;
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			theta_component = Data.PointSourcesPtr->ElectricDipole(i)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			theta_component *= exp(ii*kk_g(l,d1,d2)*(SourceX*dx*SinTCosP(l,d1,d2)+SourceY*dx*SinTSinP(l,d1,d2)+SourceZ*dx*CosT(l,d1,d2)));

			//apply the final prefactors
			if (Orientation=="x_directed")
			{
				theta_component *= K*(ii)*kk_g(l,d1,d2)*(-Z_0)*CosTCosP(l,d1,d2);
			}
			else if (Orientation=="y_directed")
			{
				theta_component *= K*(ii)*kk_g(l,d1,d2)*(-Z_0)*CosTSinP(l,d1,d2);
			}
			else if (Orientation=="z_directed")
			{
				theta_component *= K*(ii)*kk_g(l,d1,d2)*(-Z_0)*(-SinT(l,d1,d2));
			}

			//add contribution of this dipole
			total_theta_component += theta_component;
		}
	}

	return total_theta_component;
// 	return 0;
}

complex<double> Ctr_pd_fs::TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2)
//Theoretical phi-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	total_phi_component = 0;
	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		for (int i=0; i<Data.PointSourcesPtr->NumberOfElectricDipoles(); i++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(i)->Orientation;
			x = Data.PointSourcesPtr->ElectricDipole(i)->x;
			y = Data.PointSourcesPtr->ElectricDipole(i)->y;
			z = Data.PointSourcesPtr->ElectricDipole(i)->z;
			if (Orientation=="x_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX+0.5;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ;
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY+0.5;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ;
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-FarFieldOriginZ+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			phi_component = Data.PointSourcesPtr->ElectricDipole(i)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			phi_component *= exp(ii*kk_g(l,d1,d2)*(SourceX*dx*SinTCosP(l,d1,d2)+SourceY*dx*SinTSinP(l,d1,d2)+SourceZ*dx*CosT(l,d1,d2)));

			if (Orientation=="x_directed")
			{
				//apply the final prefactors
				phi_component *= K*(ii)*kk_g(l,d1,d2)*(-Z_0)*(-SinP(l,d1,d2));
			}
			else if (Orientation=="y_directed")
			{
				//apply the final prefactors
				phi_component *= K*(ii)*kk_g(l,d1,d2)*(-Z_0)*(CosP(l,d1,d2));
			}
			else if (Orientation=="z_directed")
			{
				//z component does not contribute
				phi_component = 0;
			}

			//add contribution of this dipole
			total_phi_component += phi_component;
		}
	}

	return total_phi_component;
}

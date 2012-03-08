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

//Defines the post-processing steps in the PHASOR-DOMAIN near-field-to-far-field transformer object "Ctr_pd_2l" in a 2-layered medium

#include "headers.h"

#include "Ctr_pd_2l.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//definition of Cpointsources needed
#include "pointsources/Cpointsources.h"

extern double dx;

extern int i,j,k;


void Ctr_pd_2l::ConstructPotential_A(const int& l, const int& d1, const int& d2)
{
	//constructs the theta and phi components of the magnetic potential A from the magnetic equivalent currents J_x,J_y,J_z
	//uses the formulas (21)-(24) in Capoglu thesis pg. 17

	//Back face of the virtual surface
	// J_y_back is stored at the nodes containing H-field components in the outer H-shell
	if (Back_Hz_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hz_lowerlimit_in_node);k<=Hz_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_y_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_y_downward_0);
				A_phi(l,d1,d2) += J_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_y_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_y_downward_0);
			}
			for (k=Hz_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hz_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_y_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_y_downward_1);
				A_phi(l,d1,d2) += J_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_y_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_y_downward_1);
			}
		}
	}
	// J_z_back is stored at the nodes containing H-field components in the outer H-shell
	if (Back_Hy_PassesThroughNode)
	{
		i=Back_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hy_lowerlimit_in_node);k<=Hy_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_z_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_z_downward_0);
			}
			for (k=Hy_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hy_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_z_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_z_downward_1);
			}
		}
	}

	//Front face of the virtual surface
	// J_y_front is stored at the nodes containing H-field components in the outer H-shell
	if (Front_Hz_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hz_lowerlimit_in_node);k<=Hz_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_y_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_y_downward_0);
				A_phi(l,d1,d2) += J_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_y_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_y_downward_0);
			}
			for (k=Hz_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hz_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_y_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_y_downward_1);
				A_phi(l,d1,d2) += J_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_y_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_y_downward_1);
			}
		}
	}
	// J_z_front is stored at the nodes containing H-field components in the outer H-shell
	if (Front_Hy_PassesThroughNode)
	{
		i=Front_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hy_lowerlimit_in_node);k<=Hy_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_z_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_z_downward_0);
			}
			for (k=Hy_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hy_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_z_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_z_downward_1);
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
			for (k=max(LowerHalfSpaceUpper+1,Hz_lowerlimit_in_node);k<=Hz_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_x_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_x_downward_0);
				A_phi(l,d1,d2) += J_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_x_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_x_downward_0);
			}
			for (k=Hz_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hz_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_x_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_x_downward_1);
				A_phi(l,d1,d2) += J_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_x_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_x_downward_1);
			}
		}
	}
	// J_z_left is stored at the nodes containing H-field components in the outer H-shell
	if (Left_Hx_PassesThroughNode)
	{
		j=Left_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hx_lowerlimit_in_node);k<=Hx_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_z_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_z_downward_0);
			}
			for (k=Hx_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hx_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_z_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_z_downward_1);
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
			for (k=max(LowerHalfSpaceUpper+1,Hz_lowerlimit_in_node);k<=Hz_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_x_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_x_downward_0);
				A_phi(l,d1,d2) += J_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_x_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_x_downward_0);
			}
			for (k=Hz_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hz_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				A_theta(l,d1,d2) += J_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_x_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_x_downward_1);
				A_phi(l,d1,d2) += J_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_x_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_x_downward_1);
			}
		}
	}
	// J_z_right is stored at the nodes containing H-field components in the outer H-shell
	if (Right_Hx_PassesThroughNode)
	{
		j=Right_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Hx_lowerlimit_in_node);k<=Hx_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_z_upward_0
						 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_z_downward_0);
			}
			for (k=Hx_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Hx_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				A_theta(l,d1,d2) += J_z_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
						*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_z_upward_1
						 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_z_downward_1);
			}
		}
	}

	//Bottom face of the virtual surface
	// J_x_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hy_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1);
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
			{
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_x_downward_0);
					A_phi(l,d1,d2) += J_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_x_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
			{
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_x_downward_1);
					A_phi(l,d1,d2) += J_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_x_downward_1);
				}
			}
		}
	}
	// J_y_lower is stored at the nodes containing H-field components in the outer H-shell
	if (Lower_Hx_PassesThroughNode)
	{
		k=Lower_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1);
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
			{
				for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_y_downward_0);
					A_phi(l,d1,d2) += J_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_y_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
			{
				for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_y_downward_1);
					A_phi(l,d1,d2) += J_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_y_downward_1);
				}
			}
		}
	}

	//Upper face of the virtual surface
	// J_x_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hy_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1);
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
			{
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_x_downward_0);
					A_phi(l,d1,d2) += J_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_x_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
			{
				for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_x_downward_1);
					A_phi(l,d1,d2) += J_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_x_downward_1);
				}
			}
		}
	}
	// J_y_upper is stored at the nodes containing H-field components in the outer H-shell
	if (Upper_Hx_PassesThroughNode)
	{
		k=Upper_E_shell_index;  //The equivalent electric current J coincides with the E field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1);
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
			{
				for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_theta_J_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_theta_J_y_downward_0);
					A_phi(l,d1,d2) += J_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_A_phi_J_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_A_phi_J_y_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
			{
				for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
				{
					A_theta(l,d1,d2) += J_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_theta_J_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_theta_J_y_downward_1);
					A_phi(l,d1,d2) += J_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_A_phi_J_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_A_phi_J_y_downward_1);
				}
			}
		}
	}
}

void Ctr_pd_2l::ConstructPotential_F(const int& l, const int& d1, const int& d2)
{
	//constructs the theta and phi components of the magnetic potential F from the magnetic equivalent currents M_x,M_y,M_z
	//uses the formulas (21)-(24) in Capoglu thesis pg. 17

	//Back face of the virtual surface
	// M_y_back is stored at the nodes containing E-field components in the outer E-shell
	if (Back_Ez_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Ez_lowerlimit_in_node);k<=Ez_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_y_downward_0);
				F_phi(l,d1,d2) += M_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_y_downward_0);
			}
			for (k=Ez_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ez_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_y_downward_1);
				F_phi(l,d1,d2) += M_y_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_y_downward_1);
			}
		}
	}
	// M_z_back is stored at the nodes containing E-field components in the outer E-shell
	if (Back_Ey_PassesThroughNode)
	{
		i=Back_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Ey_lowerlimit_in_node);k<=Ey_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_z_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_z_downward_0);
			}
			for (k=Ey_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ey_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_back(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_z_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_z_downward_1);
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
			for (k=max(LowerHalfSpaceUpper+1,Ez_lowerlimit_in_node);k<=Ez_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_y_downward_0);
				F_phi(l,d1,d2) += M_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_y_downward_0);
			}
			for (k=Ez_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ez_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_y_downward_1);
				F_phi(l,d1,d2) += M_y_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_y_downward_1);
			}
		}
	}
	// M_z_front is stored at the nodes containing E-field components in the outer E-shell
	if (Front_Ey_PassesThroughNode)
	{
		i=Front_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Ey_lowerlimit_in_node);k<=Ey_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_z_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_z_downward_0);
			}
			for (k=Ey_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ey_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_front(l,j,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_z_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_z_downward_1);
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
			for (k=max(LowerHalfSpaceUpper+1,Ez_lowerlimit_in_node);k<=Ez_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_x_downward_0);
				F_phi(l,d1,d2) += M_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_x_downward_0);
			}
			for (k=Ez_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ez_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_x_downward_1);
				F_phi(l,d1,d2) += M_x_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_x_downward_1);
			}
		}
	}
	// M_z_left is stored at the nodes containing E-field components in the outer E-shell
	if (Left_Ex_PassesThroughNode)
	{
		j=Left_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Ex_lowerlimit_in_node);k<=Ex_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_z_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_z_downward_0);
			}
			for (k=Ex_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ex_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_left(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_z_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_z_downward_1);
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
			for (k=max(LowerHalfSpaceUpper+1,Ez_lowerlimit_in_node);k<=Ez_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_x_downward_0);
				F_phi(l,d1,d2) += M_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_x_downward_0);
			}
			for (k=Ez_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ez_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
				F_theta(l,d1,d2) += M_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_x_downward_1);
				F_phi(l,d1,d2) += M_x_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_x_downward_1);
			}
		}
	}
	// M_z_right is stored at the nodes containing E-field components in the outer E-shell
	if (Right_Ex_PassesThroughNode)
	{
		j=Right_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=max(LowerHalfSpaceUpper+1,Ex_lowerlimit_in_node);k<=Ex_upperlimit_in_node;k++)
			{//upper half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_z_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_z_downward_0);
			}
			for (k=Ex_lowerlimit_in_node;k<=min(LowerHalfSpaceUpper,Ex_upperlimit_in_node);k++)
			{//lower half space
				z_pr = k-(LowerHalfSpaceUpper+1);
				F_theta(l,d1,d2) += M_z_right(l,i,k)
						*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_z_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_z_downward_1);
			}
		}
	}

	//Bottom face of the virtual surface
	// M_x_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ey_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
			{
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_x_downward_0);
					F_phi(l,d1,d2) += M_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_x_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
			{
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_x_downward_1);
					F_phi(l,d1,d2) += M_x_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_x_downward_1);
				}
			}
		}
	}
	// M_y_lower is stored at the nodes containing E-field components in the outer E-shell
	if (Lower_Ex_PassesThroughNode)
	{
		k=Lower_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
			{
				for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_y_downward_0);
					F_phi(l,d1,d2) += M_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_y_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
			{
				for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_y_downward_1);
					F_phi(l,d1,d2) += M_y_lower(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_y_downward_1);
				}
			}
		}
	}

	//Upper face of the virtual surface
	// M_x_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ey_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
			{
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_x_downward_0);
					F_phi(l,d1,d2) += M_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_x_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_x_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
			{
				for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_x_downward_1);
					F_phi(l,d1,d2) += M_x_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX)*sinTcosP+(j-FarFieldOriginY+0.5)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_x_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_x_downward_1);
				}
			}
		}
	}
	// M_y_upper is stored at the nodes containing E-field components in the outer E-shell
	if (Upper_Ex_PassesThroughNode)
	{
		k=Upper_H_shell_index;  //The equivalent magnetic current M coincides with the H field on the NFFFT surface
		z_pr = k-(LowerHalfSpaceUpper+1)+0.5;
		if (k>LowerHalfSpaceUpper)
		{//upper half space
			for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
			{
				for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_theta_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_theta_M_y_downward_0);
					F_phi(l,d1,d2) += M_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt0_dx*z_pr)*Green_F_phi_M_y_upward_0
						 	 +exp(-ii*kt0_dx*z_pr)*Green_F_phi_M_y_downward_0);
				}
			}
		}
		else
		{//lower half space
			for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
			{
				for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
				{
					F_theta(l,d1,d2) += M_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_theta_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_theta_M_y_downward_1);
					F_phi(l,d1,d2) += M_y_upper(l,i,j)
							*exp(ii*ko_dx*((i-FarFieldOriginX+0.5)*sinTcosP+(j-FarFieldOriginY)*sinTsinP))
							*(exp(ii*kt1_dx*z_pr)*Green_F_phi_M_y_upward_1
						 	 +exp(-ii*kt1_dx*z_pr)*Green_F_phi_M_y_downward_1);
				}
			}
		}
	}
}

complex<double> Ctr_pd_2l::TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2)
//Theoretical theta-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	/** FIXME: Could be made more efficient **/
	Calculate_TL_GreenFunctions(l,d1,d2); //recalculate the TL parameters
	/** FIXME: Could be made more efficient **/

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
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1);
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY+0.5;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1);
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1)+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			theta_component = Data.PointSourcesPtr->ElectricDipole(i)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			if (SourceZ>=0)
			{//if dipole is in upper half space
				if (Orientation=="x_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt0_dx*SourceZ)*Green_A_theta_J_x_upward_0
						 	 +exp(-ii*kt0_dx*SourceZ)*Green_A_theta_J_x_downward_0);
				}
				else if (Orientation=="y_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt0_dx*SourceZ)*Green_A_theta_J_y_upward_0
						 	 +exp(-ii*kt0_dx*SourceZ)*Green_A_theta_J_y_downward_0);
				}
				else if (Orientation=="z_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt0_dx*SourceZ)*Green_A_theta_J_z_upward_0
						 	 +exp(-ii*kt0_dx*SourceZ)*Green_A_theta_J_z_downward_0);
				}
			}
			else
			{//if dipole is in lower half space
				if (Orientation=="x_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt1_dx*SourceZ)*Green_A_theta_J_x_upward_1
						 	 +exp(-ii*kt1_dx*SourceZ)*Green_A_theta_J_x_downward_1);
				}
				else if (Orientation=="y_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt1_dx*SourceZ)*Green_A_theta_J_y_upward_1
						 	 +exp(-ii*kt1_dx*SourceZ)*Green_A_theta_J_y_downward_1);
				}
				else if (Orientation=="z_directed")
				{
					theta_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt1_dx*SourceZ)*Green_A_theta_J_z_upward_1
						 	 +exp(-ii*kt1_dx*SourceZ)*Green_A_theta_J_z_downward_1);
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

complex<double> Ctr_pd_2l::TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2)
//Theoretical phi-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr) (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
{
	/** FIXME: Could be made more efficient **/
	Calculate_TL_GreenFunctions(l,d1,d2); //recalculate the TL parameters
	/** FIXME: Could be made more efficient **/

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
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1);
			}
			else if (Orientation=="y_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY+0.5;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1);
			}
			else if (Orientation=="z_directed")
			{
				SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
				SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
				SourceZ = Data.PointSourcesPtr->ElectricDipole(i)->z-(LowerHalfSpaceUpper+1)+0.5;
			}

			//start with the temporal phasor (Fourier transform)
			phi_component = Data.PointSourcesPtr->ElectricDipole(i)->current_moment_Fourier_component(ww(l));

			//apply the spatial phasor
			if (SourceZ>=0)
			{//if dipole is in upper half space
				if (Orientation=="x_directed")
				{
					phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt0_dx*SourceZ)*Green_A_phi_J_x_upward_0
						 	 +exp(-ii*kt0_dx*SourceZ)*Green_A_phi_J_x_downward_0);
				}
				else if (Orientation=="y_directed")
				{
					phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt0_dx*SourceZ)*Green_A_phi_J_y_upward_0
						 	 +exp(-ii*kt0_dx*SourceZ)*Green_A_phi_J_y_downward_0);
				}
				else if (Orientation=="z_directed")
				{
					phi_component *= 0;	//z-directed dipole does not contribute to the phi-component of the radiated E-field
				}
			}
			else
			{//if dipole is in lower half space
				if (Orientation=="x_directed")
				{
					phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt1_dx*SourceZ)*Green_A_phi_J_x_upward_1
						 	 +exp(-ii*kt1_dx*SourceZ)*Green_A_phi_J_x_downward_1);
				}
				else if (Orientation=="y_directed")
				{
					phi_component *= exp(ii*ko_dx*(SourceX*sinTcosP+SourceY*sinTsinP))
							*(exp(ii*kt1_dx*SourceZ)*Green_A_phi_J_y_upward_1
						 	 +exp(-ii*kt1_dx*SourceZ)*Green_A_phi_J_y_downward_1);
				}
				else if (Orientation=="z_directed")
				{
					phi_component *= 0;	//z-directed dipole does not contribute to the phi-component of the radiated E-field
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

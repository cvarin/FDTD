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

//Definition of the abstract base class "Ctr_pd" for a PHASOR-DOMAIN near-field-to-far-field transformer
//This file includes the definitions of the phasor updates due to magnetic (M) and electric (J) surface currents

#include "headers.h"

#include "Ctr_pd.h"

//definition of Ctfsf needed
#include "tfsf/Ctfsf.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern double dt;

extern Array<double,3> Ex,Hy,Hz;
extern Array<double,3> Hx,Ey,Ez;

extern int rank;

extern int i,j,k;


void Ctr_pd::UpdateFarField(const int& n)
{
	 //update the equivalent-current J.M phasors on the transformer surface
	 UpdateElectric(n);
	 UpdateMagnetic(n);

	 //update the phasors for the scattered PW components
	 UpdatePWPhasors(n);	//uses the base-class method from Ctr_pd
}

void Ctr_pd::UpdatePWPhasors(const int& n)
{//updates the phasors for the scattered PW components (all referenced to the transformer origin)
	if (Data.TFSFPtr!=NULL)
	{//don't dereference TFSFPtr if it is NULL
		if (rank==0)
		{//only the master node has to do this, since it is the one that writes out "field_array"
			Data.TFSFPtr->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);	//get the instantaneous E-field amplitudes (x and y components) of the scattered PWs

			//update the PW phasors
			for (int l=0; l<L; l++)
			{
				for (int PW=0; PW<E_x_array.size(); PW++)	// could have also used "E_y_array.size()", since they are of the same size (num. of PWs)
				{
					//"dt" factor at the end approximates continuous-time Fourier transform
					//The 1/(2pi) factor is to get rid of the same factor in the inverse transform, such that the result represents the true weight of the desired frequency component.
					E_x_phasor_array(l,PW) += 1.0/(2*M_PI)*exp(-ii*ww(l)*(n*dt+get_initial_time_value()-origindelay_array(PW)))*E_x_array(PW)*dt;
					E_y_phasor_array(l,PW) += 1.0/(2*M_PI)*exp(-ii*ww(l)*(n*dt+get_initial_time_value()-origindelay_array(PW)))*E_y_array(PW)*dt;
					// for the time positioning information (reason for n*dt), see "timing_chart.txt"
				}
			}
			//resize the field arrays to 0
			E_x_array.resize(0);
			E_y_array.resize(0);
		}
	}
}

void Ctr_pd::UpdateMagnetic(const int& n)
{
	//Using the equivalence theorem, the E-values on the virtual surface are converted to Meq-values.
	//The Meq values are half-cell away from the virtual surface
	//It must be ensured that each E-value is used in only ONE NODE.

	//Back face of the virtual surface - Use Ey
	if (Back_Ey_PassesThroughNode)
	{
		i=Back_E_shell_index;
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=Ey_lowerlimit_in_node;k<=Ey_upperlimit_in_node;k++)
			{
				Mz = Ey(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_z_back(l,j,k) += Mz*phasors_Mt(l,n);
				}
			}
		}
	}

	//Back face of the virtual surface - Use Ez
	if (Back_Ez_PassesThroughNode)
	{
		i=Back_E_shell_index;
		for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				My = -Ez(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_y_back(l,j,k) += My*phasors_Mt(l,n);
				}
			}
		}
	}

	//Front face of the virtual surface - Use Ey
	if (Front_Ey_PassesThroughNode)
	{
		i=Front_E_shell_index;
		for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
		{
			for (k=Ey_lowerlimit_in_node;k<=Ey_upperlimit_in_node;k++)
			{
				Mz = -Ey(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_z_front(l,j,k) += Mz*phasors_Mt(l,n);
				}
			}
		}
	}

	//Front face of the virtual surface - Use Ez
	if (Front_Ez_PassesThroughNode)
	{
		i=Front_E_shell_index;
		for (j=Ez_leftlimit_in_node;j<=Ez_rightlimit_in_node;j++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				My = Ez(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_y_front(l,j,k) += My*phasors_Mt(l,n);
				}
			}
		}
	}

	//Left face of the virtual surface - Use Ex
	if (Left_Ex_PassesThroughNode)
	{
		j=Left_E_shell_index;
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=Ex_lowerlimit_in_node;k<=Ex_upperlimit_in_node;k++)
			{
				Mz = -Ex(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_z_left(l,i,k) += Mz*phasors_Mt(l,n);
				}
			}
		}
	}

	//Left face of the virtual surface - Use Ez
	if (Left_Ez_PassesThroughNode)
	{
		j=Left_E_shell_index;
		for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				Mx = Ez(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_x_left(l,i,k) += Mx*phasors_Mt(l,n);
				}
			}
		}
	}

	//Right face of the virtual surface - Use Ex
	if (Right_Ex_PassesThroughNode)
	{
		j=Right_E_shell_index;
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (k=Ex_lowerlimit_in_node;k<=Ex_upperlimit_in_node;k++)
			{
				Mz = Ex(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_z_right(l,i,k) += Mz*phasors_Mt(l,n);
				}
			}
		}
	}

	//Right face of the virtual surface - Use Ez
	if (Right_Ez_PassesThroughNode)
	{
		j=Right_E_shell_index;
		for (i=Ez_backlimit_in_node;i<=Ez_frontlimit_in_node;i++)
		{
			for (k=Ez_lowerlimit_in_node;k<=Ez_upperlimit_in_node;k++)
			{
				Mx = -Ez(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_x_right(l,i,k) += Mx*phasors_Mt(l,n);
				}
			}
		}
	}

	//Lower face of the virtual surface - Use Ex
	if (Lower_Ex_PassesThroughNode)
	{
		k=Lower_E_shell_index;
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
			{
				My = Ex(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_y_lower(l,i,j) += My*phasors_Mt(l,n);
				}
			}
		}
	}

	//Lower face of the virtual surface - Use Ey
	if (Lower_Ey_PassesThroughNode)
	{
		k=Lower_E_shell_index;
		for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
		{
			for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
			{
				Mx = -Ey(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_x_lower(l,i,j) += Mx*phasors_Mt(l,n);
				}
			}
		}
	}

	//Upper face of the virtual surface - Use Ex
	if (Upper_Ex_PassesThroughNode)
	{
		k=Upper_E_shell_index;
		for (i=Ex_backlimit_in_node;i<=Ex_frontlimit_in_node;i++)
		{
			for (j=Ex_leftlimit_in_node;j<=Ex_rightlimit_in_node;j++)
			{
				My = -Ex(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_y_upper(l,i,j) += My*phasors_Mt(l,n);
				}
			}
		}
	}

	//Upper face of the virtual surface - Use Ey
	if (Upper_Ey_PassesThroughNode)
	{
		k=Upper_E_shell_index;
		for (i=Ey_backlimit_in_node;i<=Ey_frontlimit_in_node;i++)
		{
			for (j=Ey_leftlimit_in_node;j<=Ey_rightlimit_in_node;j++)
			{
				Mx = Ey(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					M_x_upper(l,i,j) += Mx*phasors_Mt(l,n);
				}
			}
		}
	}
}

void Ctr_pd::UpdateElectric(const int& n)
{
	//Using the equivalence theorem, the H-values half-cell away from the virtual surface are converted to Jeq-values.
	//The Jeq values coincide with the E-field values on the virtual surface
	//It must be ensured that each H-value is used in only ONE NODE.

	//Back face of the virtual surface - Use Hy
	if (Back_Hy_PassesThroughNode)
	{
		i=Back_H_shell_index;
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=Hy_lowerlimit_in_node;k<=Hy_upperlimit_in_node;k++)
			{
				Jz = -Hy(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_z_back(l,j,k) += Jz*phasors_Jt(l,n);
				}
			}
		}
	}

	//Back face of the virtual surface - Use Hz
	if (Back_Hz_PassesThroughNode)
	{
		i=Back_H_shell_index;
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				Jy = Hz(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_y_back(l,j,k) += Jy*phasors_Jt(l,n);
				}
			}
		}
	}

	//Front face of the virtual surface - Use Hy
	if (Front_Hy_PassesThroughNode)
	{
		i=Front_H_shell_index;
		for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
		{
			for (k=Hy_lowerlimit_in_node;k<=Hy_upperlimit_in_node;k++)
			{
				Jz = Hy(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_z_front(l,j,k) += Jz*phasors_Jt(l,n);
				}
			}
		}
	}

	//Front face of the virtual surface - Use Hz
	if (Front_Hz_PassesThroughNode)
	{
		i=Front_H_shell_index;
		for (j=Hz_leftlimit_in_node;j<=Hz_rightlimit_in_node;j++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				Jy = -Hz(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_y_front(l,j,k) += Jy*phasors_Jt(l,n);
				}
			}
		}
	}

	//Left face of the virtual surface - Use Hx
	if (Left_Hx_PassesThroughNode)
	{
		j=Left_H_shell_index;
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=Hx_lowerlimit_in_node;k<=Hx_upperlimit_in_node;k++)
			{
				Jz = Hx(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_z_left(l,i,k) += Jz*phasors_Jt(l,n);
				}
			}
		}
	}

	//Left face of the virtual surface - Use Hz
	if (Left_Hz_PassesThroughNode)
	{
		j=Left_H_shell_index;
		for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				Jx = -Hz(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_x_left(l,i,k) += Jx*phasors_Jt(l,n);
				}
			}
		}
	}

	//Right face of the virtual surface - Use Hx
	if (Right_Hx_PassesThroughNode)
	{
		j=Right_H_shell_index;
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (k=Hx_lowerlimit_in_node;k<=Hx_upperlimit_in_node;k++)
			{
				Jz = -Hx(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_z_right(l,i,k) += Jz*phasors_Jt(l,n);
				}
			}
		}
	}

	//Right face of the virtual surface - Use Hz
	if (Right_Hz_PassesThroughNode)
	{
		j=Right_H_shell_index;
		for (i=Hz_backlimit_in_node;i<=Hz_frontlimit_in_node;i++)
		{
			for (k=Hz_lowerlimit_in_node;k<=Hz_upperlimit_in_node;k++)
			{
				Jx = Hz(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_x_right(l,i,k) += Jx*phasors_Jt(l,n);
				}
			}
		}
	}

	//Lower face of the virtual surface - Use Hx
	if (Lower_Hx_PassesThroughNode)
	{
		k=Lower_H_shell_index;
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
			{
				Jy = -Hx(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_y_lower(l,i,j) += Jy*phasors_Jt(l,n);
				}
			}
		}
	}

	//Lower face of the virtual surface - Use Hy
	if (Lower_Hy_PassesThroughNode)
	{
		k=Lower_H_shell_index;
		for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
		{
			for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
			{
				Jx = Hy(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_x_lower(l,i,j) += Jx*phasors_Jt(l,n);
				}
			}
		}
	}

	//Upper face of the virtual surface - Use Hx
	if (Upper_Hx_PassesThroughNode)
	{
		k=Upper_H_shell_index;
		for (i=Hx_backlimit_in_node;i<=Hx_frontlimit_in_node;i++)
		{
			for (j=Hx_leftlimit_in_node;j<=Hx_rightlimit_in_node;j++)
			{
				Jy = Hx(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_y_upper(l,i,j) += Jy*phasors_Jt(l,n);
				}
			}
		}
	}

	//Upper face of the virtual surface - Use Hy
	if (Upper_Hy_PassesThroughNode)
	{
		k=Upper_H_shell_index;
		for (i=Hy_backlimit_in_node;i<=Hy_frontlimit_in_node;i++)
		{
			for (j=Hy_leftlimit_in_node;j<=Hy_rightlimit_in_node;j++)
			{
				Jx = -Hy(i,j,k);
				for (l=0;l<ww.size(); l++)
				{
					J_x_upper(l,i,j) += Jx*phasors_Jt(l,n);
				}
			}
		}
	}
}

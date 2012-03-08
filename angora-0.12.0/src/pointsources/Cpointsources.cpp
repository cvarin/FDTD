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

#include "headers.h"

#include "Cpointsources.h"

//for the definition of MaterialId
#include "material_id.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern double dt,dx;
extern int NSTEPS;
extern int OriginX,OriginY,OriginZ;

extern Array<double,3> Ex,Ey,Ez;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<double,1> eps_x,eps_y,eps_z;

extern double epsilon_r_max,mu_r_max;

extern int GridIndex;
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;


int Cpointsources::AddElectricDipole(const int& xPos,  const int& yPos, const int& zPos,
		 const string& myOrientation, const double& myj0, const_Cwf_shared_ptr mywaveform)
{
	boost::shared_ptr<Celectricdipole> new_dipole_ptr(new Celectricdipole(xPos,yPos,zPos,myOrientation,myj0,mywaveform,NumberOfElectricDipoles()));
	ElectricDipoles.push_back(new_dipole_ptr);
}

void Cpointsources::ApplySources(const int& n)
{
	for (int i=0; i<NumberOfElectricDipoles(); i++)
	{
		ElectricDipoles[i]->ApplyElectricDipole(n);		//apply the electric dipole
	}
}

Celectricdipole::Celectricdipole(const int& xPos,  const int& yPos, const int& zPos,
		 const string& myOrientation, const double& myj0, const_Cwf_shared_ptr mywaveform, const int& Index)
	:	x(xPos),y(yPos),z(zPos),Orientation(myOrientation),j0(myj0),waveform(mywaveform),
		ElectricDipoleIndex(Index),ElectricDipoleBelongsToNode(false)
{
	if (Orientation=="x_directed")
	{
		if ((x>=iback)&&(x<=ifront)&&(y>=jleft)&&(y<=jright+1)&&(z>=klower)&&(z<=kupper+1))
		{
			ElectricDipoleBelongsToNode=true;
		}
	}
	else if (Orientation=="y_directed")
	{
		if ((x>=iback)&&(x<=ifront+1)&&(y>=jleft)&&(y<=jright)&&(z>=klower)&&(z<=kupper+1))
		{
			ElectricDipoleBelongsToNode=true;
		}
	}
	else if (Orientation=="z_directed")
	{
		if ((x>=iback)&&(x<=ifront+1)&&(y>=jleft)&&(y<=jright+1)&&(z>=klower)&&(z<=kupper))
		{
			ElectricDipoleBelongsToNode=true;
		}
	}
	else
	{
		if (rank==0)
		{
			cout << "Invalid dipole orientation (" << Orientation << ") for electric dipole "
				<< ElectricDipoleIndex << " in grid " << GridIndex << endl;
			exit(-1);
		}
	}

	double lambda_min = c/sqrt(epsilon_r_max*mu_r_max)/(waveform->w_max_40()/2/M_PI);		//minimum wavelength corresponding to w_max
	double L = lambda_min/dx;			//number of grid cells per minimum wavelength
	double L_req = 15;	//temporary value, may change the assignment method later
	if (rank==0)
	{
		if (L<L_req)
		{
			cout << endl << "Warning: Not enough cells per minimum wavelength (" << L << "<" << L_req << ") in grid "
				<< GridIndex << "." << endl;
		}
	}

	//if the there is a sudden jump in the excitation waveform, set the initial time value back
	if (waveform->starting_time()<get_initial_time_value())
	{
		set_initial_time_value(waveform->starting_time());
	}

	//determine the permittivity at the position of the source (used below for E0 & J0)
	SourceCell = x,y,z;	//the cell to which the source belongs

	//1/dx^3 converts from current moment (dimension: A*m) to 3D current density (dimension:A/m^2)
	//NOTE: This is different from 2-D, where the conversion factor was 1/dx^2.
	moment_to_current_density = 1/(dx*dx*dx);
}

double Celectricdipole::current_moment_value(const double& t)
{
	return j0*waveform->Value(t);
}

complex<double> Celectricdipole::current_moment_Fourier_transform(const double& w)
{
	return j0*waveform->FourierTransform(w);
}

complex<double> Celectricdipole::current_moment_Fourier_component(const double& w)
{
	return j0*waveform->FourierComponent(w);
}

Cwf_shared_ptr Celectricdipole::current_moment_waveform()
{
	Cwf_shared_ptr copied_wf_ptr(waveform->clone());
	copied_wf_ptr->multiply_amplitude_by(j0); //apply the extra current-moment amplitude
	return copied_wf_ptr;
}

double Celectricdipole::waveform_value(const double& n)
{
	return j0*(waveform->Value((n+0.5)*dt+get_initial_time_value()));
	//0.5dt delay is because J is 0.5dt ahead of the E-field component that it updates within a cycle.
}

void Celectricdipole::ApplyElectricDipole(const int& n)
{
	if (ElectricDipoleBelongsToNode)
	{
		if (Orientation=="x_directed")
		{
			Ex(SourceCell) += - Cb_X(Media_Ex(SourceCell))*dx*moment_to_current_density*waveform_value(n);
		}
		else if (Orientation=="y_directed")
		{
			Ey(SourceCell) += - Cb_Y(Media_Ey(SourceCell))*dx*moment_to_current_density*waveform_value(n);
		}
		else if (Orientation=="z_directed")
		{
			Ez(SourceCell) += - Cb_Z(Media_Ez(SourceCell))*dx*moment_to_current_density*waveform_value(n);
		}
	}
}

// /*********************************************************************/
// /**************  CURRENT MOMENT CALCULATOR  **************************/
// /*********************************************************************/
//
// //The following static method can be used to obtain the current moment needed to create an electric field maximum of Emax at the source point
// double Celectricdipole::CurrentMomentCalculator(const string& WaveShape, const double& tau, const double& dx,
// 							   const double& epsilon_r, const double& Emax)
// {
// 	if (WaveShape=="gaussian")
// 	{
// 		return pow(dx,3)*Emax*(3*epsilon_r*epsilon_0)/(sqrt(2*M_PI)*tau);
// 	}
// 	else if (WaveShape=="diffgaussian")
// 	{
// 		return pow(dx,3)*Emax*(3*epsilon_r*epsilon_0)*(exp(-0.5)/tau);
// 	}
// 	else if (WaveShape=="doublediffgaussian")
// 	{
// 		return pow(dx,3)*Emax*(3*epsilon_r*epsilon_0)*(exp(0.5)*tau)*(1/pow(tau,2));
// 	}
// 	else
// 	{
// 		if (rank==0)
// 		{
// 			cout << "Invalid wave shape (" << WaveShape << ") for ""CurrentMomentCalculator"" in node " << GridIndex << endl;
// 			exit(-1);
// 		}
// 	}
// }

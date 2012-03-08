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

//Definition of the abstract base class "Cpw" for a TF/SF plane-wave source

#include "headers.h"

#include "Cpw.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

extern double dx,courant;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern double c_upper;
extern double epsilon_r_max,mu_r_max;

extern int GridIndex;
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern double max_field_value;
extern bool max_field_value_set_in_configfile;


Cpw::Cpw(const PWDataType& MyData)
		:Data(MyData)
{
	THETA_INC = Data.THETA;		//can be anything
	PHI_INC = Data.PHI;		//can be anything

	PSI = Data.PSI; //(for LP incident field, if applicable) angle of linearly polarized incident electric field w.r.t.
					// (k_inc)x(z) in clockwise direction (k_inc points away from the clock surface)

	SinT = sin(THETA_INC); CosT = cos(THETA_INC);
	SinP = sin(PHI_INC); CosP = cos(PHI_INC);
	SinPsi = sin(PSI); CosPsi = cos(PSI);

	//Unit vector in the direction of incidence
	k_inc = (-CosP*SinT),(-SinP*SinT),(-CosT);
	//Lateral unit vector in the direction of incidence
	if (SinT>=0)	//theta=0 results in the same direction as theta=eps, where eps is positive and arbitrarily small
		k_inc_lateral = (-CosP),(-SinP),(0);
	else
		k_inc_lateral = (CosP),(SinP),(0);
	unit_z = (0.0),(0.0),(1.0);

	//components of k_inc
	k_inc_x = k_inc(firstDim);
	k_inc_y = k_inc(secondDim);
	k_inc_z = k_inc(thirdDim);

	//The TF/SF boundary is between the shell of E-field on the outer boundary of the box and the shell of H-field immediately outside the box.
	PWBackX = NPML+Data.PWMarginBackX+1;
	PWFrontX = NCELLS_X+NPML-Data.PWMarginFrontX;
	if (PWBackX>PWFrontX)
	{//the bounds are invalid, throw exception
		throw InvalidTFSFBounds();
	}
	PWLeftY = NPML+Data.PWMarginLeftY+1;
	PWRightY = NCELLS_Y+NPML-Data.PWMarginRightY;
	if (PWLeftY>PWRightY)
	{//the bounds are invalid, throw exception
		throw InvalidTFSFBounds();
	}
	PWLowerZ = NPML+Data.PWMarginLowerZ+1;
	PWUpperZ = NCELLS_Z+NPML-Data.PWMarginUpperZ;
	if (PWLowerZ>PWUpperZ)
	{//the bounds are invalid, throw exception
		throw InvalidTFSFBounds();
	}
	//Minimum and maximum indices cells contained in the TFSF box and in the current node
	TFSF_min_x = max(iback,PWBackX);
	TFSF_max_x = min(ifront,PWFrontX);
	TFSF_min_y = max(jleft,PWLeftY);
	TFSF_max_y = min(jright,PWRightY);
	TFSF_min_z = max(klower,PWLowerZ);
	TFSF_max_z = min(kupper,PWUpperZ);

	//Maximum field amplitude in the plane wave
	double PW_max_field_value = Data.waveform->A_max();	//amplitude of the electric-field waveform in the plane wave
	if (PW_max_field_value>max_field_value)
	{
		if (!max_field_value_set_in_configfile) max_field_value=PW_max_field_value;	//set the max. field amplitude, if not set in the config file
	}

	double lambda_min = c/sqrt(epsilon_r_max*mu_r_max)/(Data.waveform->w_max_40()/2/M_PI);		//minimum wavelength corresponding to w_max
	L = lambda_min/dx;			//number of grid cells per minimum wavelength
	if (rank==0)
	{
		if (L<Data.L_req)
		{
			if (Data.DisplayWarnings)
			{
				cout << endl << "Warning: Not enough cells per minimum wavelength (" << L << "<" << Data.L_req << ") in grid "
				<< GridIndex << "." << endl;
			}
			else
			{
				//implement exception later
			}
		}
		else if (Data.DisplayCellsPerLambda)
		{
			cout << endl << "Cells per minimum wavelength is (" << L << ">" << Data.L_req << ") in grid " << GridIndex << "." << endl;
		}
	}

	//the angular modulation frequency (if modulated)
	f_0 = Data.waveform->F_0();
	w_0 = Data.waveform->W_0();
}

double Cpw::GridVelocity(const double& theta, const double& phi)
{//calculates the grid velocity v(theta,phi), normalized by c
	if (w_0!=0)
	{
		double lambda_0 = c_upper/f_0;	//free-space wavelength
		double dx_norm = dx/lambda_0;	//normalized spatial step
		double N_lambda = 1/dx_norm;	//steps per wavelength
		//propagation direction cosines
		double k_x = -cos(phi)*sin(theta);
		double k_y = -sin(phi)*sin(theta);
		double k_z = -cos(theta);
		double A = dx_norm*k_x/2;
		double B = dx_norm*k_y/2;
		double C = dx_norm*k_z/2;
		double S = courant/sqrt(3.0);
		double D = 1/pow(S,2)*pow(sin(M_PI*S/N_lambda),2);

		double k_i = 2*M_PI;	//initial guess for k
		double correction = 1e6;
		double threshold = 2*M_PI/1e10;	//error threshold for convergence

		int iter=1;
		while (abs(correction)>threshold)
		{//see Taflove, eqn. (4.16a)
			correction = -(pow(sin(A*k_i),2)+pow(sin(B*k_i),2)+pow(sin(C*k_i),2)-D)/(A*sin(2*A*k_i)+B*sin(2*B*k_i)+C*sin(2*C*k_i));
			k_i = k_i + correction;
		}
		return 2*M_PI/k_i;	//the grid velocity normalized by c
	}
	else
	{
		return 1;	//if w_0=0, then grid velocity is approx. c
	}
}

void Cpw::CorrectE(const int& n)
{
	UpdateIncidentH(n);	//update Einc and Hinc until Hinc is at n+1/2
	ApplyCorrectionE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2) by projection
}

void Cpw::CorrectH(const int& n)
{
	UpdateIncidentE(n);	//update Einc and Hinc until Einc is at n+1
	ApplyCorrectionH(n);		//apply correction to the main-grid H-field (using Einc at n+1) by projection
}

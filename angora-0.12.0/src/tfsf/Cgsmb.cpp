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

//Defines the class "Cgsmb" for a TF/SF Hermite-Gaussian source in free space

#include "headers.h"

#include "Cgsmb.h"

#include "Chgb.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

extern int GridIndex;
extern int rank;

extern int number_of_layers;

//gauss-legendre quadrature rule generator
extern void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w);

//factorial function below
int factorial (int num);


Cgsmb::Cgsmb(const GSMBDataType& MyData, const int& Index)
		:Data(MyData), GSMPIndex(Index)
{
	total_num_of_modes = Data.n_mode_x*Data.n_mode_y;

	this_mode = (GridIndex % total_num_of_modes);	//index of this mode

	int global_mode_index = 0;
	for (int x_mode_index=0; x_mode_index<Data.n_mode_x; x_mode_index++)
	{
		for (int y_mode_index=0; y_mode_index<Data.n_mode_y; y_mode_index++)
		{
			if (global_mode_index==this_mode)
			{//if this is the right grid, add the Hermite-Gaussian beam mode
// cout << "Mode #" << this_mode << " is for x_mode_index=" << x_mode_index << ", y_mode_index=" << y_mode_index << endl;
				HGBDataType HGBData;

//				HGBData.n_sx = 15; //may change this later
//				HGBData.n_sy = 15; //may change this later

				HGBData.x_order = x_mode_index;
				HGBData.y_order = y_mode_index;

				HGBData.polarization = Data.polarization;
				HGBData.direction = Data.direction;

				aa = 1/(4*pow2(Data.beam_halfwidth));
				bb = 1/(2*pow2(Data.correlation_length));
				cc = sqrt(pow2(aa)+2*aa*bb);
				HGBData.beam_halfwidth = 1/sqrt(2*cc);

				HGBData.HGBMarginBackX = Data.GSMBMarginBackX;
				HGBData.HGBMarginFrontX = Data.GSMBMarginFrontX;
				HGBData.HGBMarginLeftY = Data.GSMBMarginLeftY;
				HGBData.HGBMarginRightY = Data.GSMBMarginRightY;
				HGBData.HGBMarginLowerZ = Data.GSMBMarginLowerZ;
				HGBData.HGBMarginUpperZ = Data.GSMBMarginUpperZ;
				HGBData.HGBOriginX = Data.GSMBOriginX;
				HGBData.HGBOriginY = Data.GSMBOriginY;
				HGBData.HGBOriginZ = Data.GSMBOriginZ;
				HGBData.L_req = Data.L_req;

				HGBData.DisplayWarnings = Data.DisplayWarnings;
				HGBData.DisplayCellsPerLambda = Data.DisplayCellsPerLambda;

				//see Eq. (2.15) in Starikov,Wolf (1982) and Eqs. (88)-(90) in Tervo,Setala,Friberg (2004)
				HGBData.E0 = Data.E0*sqrt(
							(2*cc/(aa+bb+cc))*
							pow(bb/(aa+bb+cc),x_mode_index+y_mode_index)*
							1.0/(pow(2.0,x_mode_index)*factorial(x_mode_index))*
							1.0/(pow(2.0,y_mode_index)*factorial(y_mode_index))
							);
				HGBData.waveform = Data.waveform;

				//Add the Hermite-Gaussian-beam source to the Gaussian-Schell-model-beam object
				boost::shared_ptr<Chgb> new_hgp_ptr(new Chgb(HGBData,NumberOfHermiteGaussianBeams()));
				HermiteGaussianBeams.push_back(new_hgp_ptr);
			}
			global_mode_index++;
		}
	}
}

void Cgsmb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the Gaussian-Schell-model beam.
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Cgsmb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the Gaussian-Schell-model beam.
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Cgsmb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWDirections(PW_THETA,PW_PHI);
	}
}

void Cgsmb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWDelaysFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Cgsmb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);
	}
}

//factorial function
int factorial(int num)
{
 if (num<0)
 {
   if (rank==0) cout << "Invalid factorial input (" << num << ")." << endl;
   exit(-1);
 }
 if (num<=1)
  return 1;
 return factorial(num-1)*num; // recursive call
}

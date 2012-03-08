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

#ifndef CTR_PD_FS_H
#define CTR_PD_FS_H

//Declaration of the class "Ctr_pd_fs" for a PHASOR-DOMAIN near-field-to-far-field transformer in free space
//Derived from the abstract base class "Ctr_pd".

#include "Ctr_pd.h"


class Ctr_pd_fs : public Ctr_pd
{
 public:
	 Ctr_pd_fs(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index);		//constructor

	 virtual void ConstructFarField();	//(virtual from Ctr_pd) Construct far-field from the stored equivalent-current phasors

 private:
 	 virtual void ConstructPotential_A(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the magnetic potential A
 	 virtual void ConstructPotential_F(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the electric potential F

	 double eps_0;	//relative permittivity of the homogeneous space
	 double mu_0;	//relative permeability of the homogeneous space
	 double c_0;	//phase velocity in the homogeneous space
	 double Z_0;	//wave impedance in the homogeneous space

	 //wavenumber array
	 Array<double,1> kk;	//in m^-1

	 //wavenumber corrected for grid-dispersion
	 //(first dim: frequency, second&third dims: theta and phi angles)
	 Array<double,3> kk_g;		//in m^-1

	 //prefactor in the far-field expressions
	 complex<double> K;

	 //Auxiliary far-field potential arrays (one for each cartesian component)
	 Array<complex<double>,3> A_x,A_y,A_z,F_x,F_y,F_z;

     //temporary variables for fast post-processing
     //sines and cosines
     double sinTcosP,sinTsinP,cosT;
     //wavenumber in the homogeneous space, multiplied by the grid spacing
     double k_dx;

	 //Theoretical far-field calculations (virtual functions inherited from Ctr)
	 virtual complex<double> TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical theta-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
	 virtual complex<double> TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical phi-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)

	 complex<double> theta_component,phi_component,total_theta_component,total_phi_component;	//field components are complex for freq.domain NFFFT
	 double temp;
	 string Orientation;
	 int x,y,z;
	 double SourceX,SourceY,SourceZ;
};

#endif

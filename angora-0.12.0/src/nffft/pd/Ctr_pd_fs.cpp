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

//Declaration of the class "Ctr_pd_fs" for a PHASOR-DOMAIN near-field-to-far-field transformer in free space

#include "headers.h"

#include "Ctr_pd_fs.h"

extern double dx;

extern int OriginX,OriginY,OriginZ;

extern int number_of_layers;

extern double epsilon_r_upper,epsilon_r_lower,mu_r_upper,mu_r_lower;
extern double c_upper;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int j,k;


Ctr_pd_fs::Ctr_pd_fs(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_pd(MyData,FarFieldFileName,Index)
{
	if (number_of_layers!=1)
	{//if the # of layers is not 1, there is a developing error, so display error message and exit (better error-reporting scheme later!?)
		if (rank==0)
		{
			cout << "Error: Ctr_pd_fs class cannot be used if number of layers is not 1." << endl;
		}
#ifndef MPI_DISABLE
		MPI_Barrier(MPI_CartSubComm);
#endif
		exit(-1);
	}

	//parameters pertaining to the homogeneous space
	eps_0 = epsilon_r_upper;	//relative permittivity in the homogeneous space
	mu_0 = mu_r_upper;	//relative permeability in the homogeneous space
	c_0 = c/sqrt(eps_0*mu_0);				//phase velocity in the homogeneous space
	Z_0 = eta_0*sqrt(mu_0/eps_0);			//wave impedance in the homogeneous space

	//construct wavenumber array
	kk.resize(L);
	kk = ww/c_0;		//c_0 converts from free-space to homogeneous medium

	//construct grid-dispersion-corrected wavenumber array
	kk_g.resize(L,D1,D2);
	for (int l=0; l<L; l++)
	{
		for (int d1=0; d1<D1; d1++)
		{
			for (int d2=0; d2<D2; d2++)
			{
				kk_g(l,d1,d2) = kk(l)/GridVelocity(ww(l),THETA(l,d1,d2),PHI(l,d1,d2));		//correct for the velocity factor
			}
		}
	}

	//prefactor K in the far-field expressions
	K = 1/(4*M_PI);

	//allocate and initialize the potential arrays (cartesian components)
	A_x.resize(L,D1,D2);
	A_x=0;
	A_y.resize(L,D1,D2);
	A_y=0;
	A_z.resize(L,D1,D2);
	A_z=0;
	F_x.resize(L,D1,D2);
	F_x=0;
	F_y.resize(L,D1,D2);
	F_y=0;
	F_z.resize(L,D1,D2);
	F_z=0;
}

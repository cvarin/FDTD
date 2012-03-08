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

//Definition of the class "Ctr_pd_2l" for a PHASOR-DOMAIN near-field-to-far-field transformer in a 2-layered medium

#include "headers.h"

#include "Ctr_pd_2l.h"

#include "material_id.h"

extern int number_of_layers;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;

extern Array<ElectricMaterialIndexType_X,1> LayerElectricMaterial_X;
extern Array<ElectricMaterialIndexType_Y,1> LayerElectricMaterial_Y;
extern Array<ElectricMaterialIndexType_Z,1> LayerElectricMaterial_Z;
extern Array<MagneticMaterialIndexType_X,1> LayerMagneticMaterial_X;
extern Array<MagneticMaterialIndexType_Y,1> LayerMagneticMaterial_Y;
extern Array<MagneticMaterialIndexType_Z,1> LayerMagneticMaterial_Z;

extern Array<int,1> LayerLowerZIndices;

extern Array<double,1> cond_e_x;	//this transformer works only for isotropic materials, may generalize later

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;


Ctr_pd_2l::Ctr_pd_2l(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_pd_ml(MyData,FarFieldFileName,Index)
{
	if (number_of_layers!=2)
	{//if the # of layers is not 2, there is a developing error, so display error message and exit (better error-reporting scheme later!?)
		if (rank==0)
		{
			cout << "Error: Ctr_pd_2l class cannot be used if number of layers is not 2." << endl;
		}
#ifndef MPI_DISABLE
		MPI_Barrier(MPI_CartSubComm);
#endif
		exit(-1);
	}

	//take the highest index in the lower half space from the layering info
	LowerHalfSpaceUpper = LayerLowerZIndices(1)-1;	//layer 1 is the uppermost layer

	//z-index of the cell used as the origin of the NFFFT until post-processing (inherited from Ctr_pd_ml)
	TemporaryFarFieldOriginZ = LowerHalfSpaceUpper+1;

	//These parameters may later be taken automatically from the layering structure. They are currently being taken from the global variables epsilon_r_upper,epsilon_r_lower.
	//relative permittivities of the upper and lower half spaces
	epsilon_r_0 = epsilon_r_upper;
	epsilon_r_1 = epsilon_r_lower;
	//relative permeabilities of the upper and lower half spaces
	mu_r_0 = mu_r_upper;
	mu_r_1 = mu_r_lower;
	//electrical conductivities of the upper and lower half spaces
	cond_e_0 = cond_e_x(LayerElectricMaterial_X(number_of_layers-1)); /** isotropy assumed! **/	//electrical conductivity of the upper layer;
	cond_e_1 = cond_e_x(LayerElectricMaterial_X(0)); /** isotropy assumed! **/	//electrical conductivity of the lower layer;
}

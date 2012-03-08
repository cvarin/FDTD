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

//Definition of the class "Ctr_pd_3l" for a PHASOR-DOMAIN near-field-to-far-field transformer in a 3-layered medium

#include "headers.h"

#include "Ctr_pd_3l.h"

#include "material_id.h"

extern int OriginZ;

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
extern Array<int,1> LayerThicknesses;
extern Array<bool,1> IsLayerGrounded;

extern Array<double,1> eps_x,mu_x,cond_e_x;	//this transformer works only for isotropic materials, may generalize later

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;


Ctr_pd_3l::Ctr_pd_3l(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_pd_ml(MyData,FarFieldFileName,Index)
{
	if (number_of_layers==1)
	{//if the # of layers is 1, issue a developer warning, since there is a more efficient class (Ctr_pd_fs) for free space
		if (rank==0)
		{
			cout << "Developer warning: Simulation space is homogeneous; consider using the Ctr_pd_fs class." << endl;
		}
	}
	if (number_of_layers==2)
	{//if the # of layers is 2, issue a developer warning, since there is a more efficient class (Ctr_pd_2l) for 2-layered media
		if (rank==0)
		{
			cout << "Developer warning: Only two layers present; consider using the Ctr_pd_2l class." << endl;
		}
	}
	if (number_of_layers>3)
	{//if the # of layers is not 3, isse a developer error and exit (better error-reporting scheme later!?)
		if (rank==0)
		{
			cout << "Developer error: Ctr_pd_3l class cannot be used if the number of layers is more than three." << endl;
		}
#ifndef MPI_DISABLE
		MPI_Barrier(MPI_CartSubComm);
#endif
		exit(-1);
	}

	//These parameters may later be taken automatically from the layering structure. They are currently being taken from the global variables epsilon_r_upper,epsilon_r_lower.
	//relative permittivities of the upper and lower half spaces
	epsilon_r_0 = epsilon_r_upper;
	epsilon_r_2 = epsilon_r_lower;
	//relative permeabilities of the upper and lower half spaces
	mu_r_0 = mu_r_upper;
	mu_r_2 = mu_r_lower;
	//electrical conductivities of the upper and lower half spaces
	cond_e_0 = cond_e_x(LayerElectricMaterial_X(number_of_layers-1)); /** isotropy assumed! **/	//electrical conductivity of the upper layer;
	cond_e_2 = cond_e_x(LayerElectricMaterial_X(0)); /** isotropy assumed! **/	//electrical conductivity of the lower layer;

	//take the position, thickness, permittivity and permeability of the middle layer from the layering info
	if (number_of_layers==3)
	{//if the # of layers is 3, the middle layer is a dielectric slab
		if (IsLayerGrounded(2))
		{//if the uppermost layer is grounded, then we basically have a 2-layered medium
			int MiddleLayerIndex = 2; // dummy middle layer has the same properties as the uppermost layer
			MiddleLayerPos = LayerLowerZIndices(MiddleLayerIndex);	//dummy middle layer of thickness 1 above the material interface: index of the highest (and lowest) cell in middle layer is the same as the lowest z-index cell in the upper layer
			MiddleLayerThickness = 1;	//dummy middle layer of thickness 1
			epsilon_r_1 = eps_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permittivity of the upper layer
			cond_e_1 = cond_e_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//electrical conductivity of the upper layer
			mu_r_1 = mu_x(LayerMagneticMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permeability of the upper layer
			grounded = true; //the dummy middle layer is grounded
		}
		else
		{//uppermost layer is not grounded (but the middle layer still could be)
			int MiddleLayerIndex = 1; // middle layer has index 1
			MiddleLayerPos = LayerLowerZIndices(MiddleLayerIndex+1)-1;	//index of the highest cell in the middle layer
			MiddleLayerThickness = LayerThicknesses(MiddleLayerIndex);	//thickness of the middle layer
			epsilon_r_1 = eps_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permittivity of the middle layer
			cond_e_1 = cond_e_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//electrical conductivity of the middle layer
			mu_r_1 = mu_x(LayerMagneticMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permeability of the middle layer
			grounded = IsLayerGrounded(MiddleLayerIndex);	//is the middle layer grounded?
		}
	}
	else if (number_of_layers==2)
	{
		//assign a "dummy" middle layer of thickness 1 above the material interface
		int MiddleLayerIndex = 1; // dummy middle layer has the same properties as the upper layer
		MiddleLayerPos = LayerLowerZIndices(MiddleLayerIndex);	//dummy middle layer of thickness 1 above the material interface: index of the highest (and lowest) cell in middle layer is the same as the lowest z-index cell in the upper layer
		MiddleLayerThickness = 1;	//dummy middle layer of thickness 1
		epsilon_r_1 = eps_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permittivity of the upper layer
		cond_e_1 = cond_e_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//electrical conductivity of the upper layer
		mu_r_1 = mu_x(LayerMagneticMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permittivity of the upper layer
		grounded = IsLayerGrounded(MiddleLayerIndex); //the upper layer may be grounded
	}
	else if (number_of_layers==1)
	{
		//assign a "dummy" middle layer of thickness 1 above the grid origin
		int MiddleLayerIndex = 0; // dummy middle layer has the same properties as the homogeneous space
		MiddleLayerPos = OriginZ;	//dummy middle layer of thickness 1 above the grid origin:  index of the highest (and lowest) cell in middle layer is the same as OriginZ
		MiddleLayerThickness = 1; //dummy middle layer of thickness 1
		epsilon_r_1 = eps_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permittivity of the homogeneous space
		cond_e_1 = cond_e_x(LayerElectricMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//electrical conductivity of the homogeneous space
		mu_r_1 = mu_x(LayerMagneticMaterial_X(MiddleLayerIndex)); /** isotropy assumed! **/	//permeability of the homogeneous space
		grounded = false;	//the dummy middle layer cannot be grounded, since there is only one layer
	}

	//z-index of the cell used as the origin of the NFFFT until post-processing (inherited from Ctr_pd_ml)
	TemporaryFarFieldOriginZ = MiddleLayerPos+1;
}

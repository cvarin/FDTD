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

//Definition of the class "Ctr_pd_ml" for a PHASOR-DOMAIN near-field-to-far-field transformer in a general N-layered medium

#include "Ctr_pd_ml.h"

#include "material_id.h"

extern double dx;

extern int number_of_layers;

extern Array<ElectricMaterialIndexType_X,1> LayerElectricMaterial_X;
extern Array<ElectricMaterialIndexType_Y,1> LayerElectricMaterial_Y;
extern Array<ElectricMaterialIndexType_Z,1> LayerElectricMaterial_Z;
extern Array<MagneticMaterialIndexType_X,1> LayerMagneticMaterial_X;
extern Array<MagneticMaterialIndexType_Y,1> LayerMagneticMaterial_Y;
extern Array<MagneticMaterialIndexType_Z,1> LayerMagneticMaterial_Z;

extern int number_of_layers;
extern Array<int,1> LayerLowerZIndices;
extern Array<int,1> LayerThicknesses;

extern Array<double,1> eps_x,mu_x,cond_e_x;	//this transformer works only for isotropic materials, may generalize later

namespace{
	int pln,interf;
};


Ctr_pd_ml::Ctr_pd_ml(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_pd(MyData,FarFieldFileName,Index)
{
	threshold = 1e-5;	//what is the threshold for sin(theta) that determines upper or lower half space?

	recursion_max_error = 1e-6;

	//number of interfaces in the layering structure
	number_of_interfaces = number_of_layers-1;

	//allocate the arrays
	epsilon_r.resize(number_of_layers);
	mu_r.resize(number_of_layers);
	cond_e.resize(number_of_layers);
	grounded.resize(number_of_layers);

	eps_r.resize(number_of_layers);
	m_r.resize(number_of_layers);

//	obliquityfactor.resize(number_of_layers);
	kt_g.resize(number_of_layers);

	eta_e.resize(number_of_layers);
	eta_h.resize(number_of_layers);

	/** these are defined at interfaces between layers, therefore they have size (#-of-layers-1) **/
	Z_e.resize(number_of_interfaces);
	Z_h.resize(number_of_interfaces);

	interface_coord.resize(Range(-1,number_of_interfaces)); //z_-1 and z_N-1 are defined, but they are pretty much arbitrary (they cancel out in the final results)
	/** these are defined at interfaces between layers, therefore they have size (#-of-layers-1) **/

	d_n.resize(number_of_layers);

	V_v_e_minus.resize(number_of_layers);
	V_v_e_plus.resize(number_of_layers);
	V_v_h_minus.resize(number_of_layers);
	V_v_h_plus.resize(number_of_layers);
	V_i_e_minus.resize(number_of_layers);
	V_i_e_plus.resize(number_of_layers);
	V_i_h_minus.resize(number_of_layers);
	V_i_h_plus.resize(number_of_layers);

	/** The following arrays *do* need to be arrays, since they are used in the post-processing step **/
	kt_dx.resize(number_of_layers);

	Green_A_theta_J_x_minus.resize(number_of_layers);
	Green_A_theta_J_x_plus.resize(number_of_layers);
	Green_A_theta_J_y_minus.resize(number_of_layers);
	Green_A_theta_J_y_plus.resize(number_of_layers);
	Green_A_theta_J_z_minus.resize(number_of_layers);
	Green_A_theta_J_z_plus.resize(number_of_layers);
	Green_A_phi_J_x_minus.resize(number_of_layers);
	Green_A_phi_J_x_plus.resize(number_of_layers);
	Green_A_phi_J_y_minus.resize(number_of_layers);
	Green_A_phi_J_y_plus.resize(number_of_layers);
	Green_F_theta_M_x_minus.resize(number_of_layers);
	Green_F_theta_M_x_plus.resize(number_of_layers);
	Green_F_theta_M_y_minus.resize(number_of_layers);
	Green_F_theta_M_y_plus.resize(number_of_layers);
	Green_F_theta_M_z_minus.resize(number_of_layers);
	Green_F_theta_M_z_plus.resize(number_of_layers);
	Green_F_phi_M_x_minus.resize(number_of_layers);
	Green_F_phi_M_x_plus.resize(number_of_layers);
	Green_F_phi_M_y_minus.resize(number_of_layers);
	Green_F_phi_M_y_plus.resize(number_of_layers);


	//store the constitutive parameters of the layers
	for (pln=0; pln<number_of_layers; pln++)
	{
		//relative permittivities of the layers
		epsilon_r(pln) = eps_x(LayerElectricMaterial_X(pln)); /** isotropy assumed! **/
		//relative permeabilities of the layers
		mu_r(pln) = mu_x(LayerMagneticMaterial_X(pln)); /** isotropy assumed! **/
		//electrical conductivities of the layers
		cond_e(pln) = cond_e_x(LayerElectricMaterial_X(pln)); /** isotropy assumed! **/
	}

	//z-index of the cell used as the origin of the NFFFT until post-processing
	TemporaryFarFieldOriginZ = LayerLowerZIndices(UPPERMOST_LAYER_INDEX); //cell with lowest z index contained in the uppermost layer

	//coordinates of the interfaces wrt the temporary NFFFT origin
	//the following two are arbitrary, but they have to be consistent with the post-processing step
	interface_coord(LOWERMOST_INTERFACE_INDEX-1) = LayerLowerZIndices(LOWERMOST_LAYER_INDEX)-TemporaryFarFieldOriginZ;
	interface_coord(UPPERMOST_INTERFACE_INDEX+1) = LayerLowerZIndices(UPPERMOST_LAYER_INDEX)+LayerThicknesses(UPPERMOST_LAYER_INDEX)-TemporaryFarFieldOriginZ;
	for (interf=0; interf<number_of_interfaces; interf++)
	{
		interface_coord(interf) = LayerLowerZIndices(interf+1)-TemporaryFarFieldOriginZ;
	}

	//thicknesses of each layer
	for (pln=0; pln<number_of_layers; pln++)
	{//values at pln=0 and pln=N-1 are arbitrary, but they should be consistent with the definitions of z_-1 and z_N-1 implied in Ctr_pd_ml_pp.cpp, hence the use of LayerThicknesses
		d_n(pln) = LayerThicknesses(pln)*dx;
	}
}

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

#ifndef MATERIAL_H
#define MATERIAL_H

//Includes declarations of types and data structures associated with material indices

//base Angora exception class
#include "angora_excp.h"

//for the C++ STL map class
#include <map>

//for definition of MaterialId
#include "material_id.h"


void AddIsotropicElectricMaterial(MaterialId& NewMaterialId, const double& epsilon_r, const double& sigma_e);
void AddDiagonalElectricMaterial(MaterialId& NewMaterialId,
								const double& epsilon_r_x, const double& epsilon_r_y, const double& epsilon_r_z,
								const double& sigma_e_x, const double& sigma_e_y, const double& sigma_e_z);
void AddIsotropicMagneticMaterial(MaterialId& NewMaterialId, const double& mu_r, const double& sigma_h);
void AddDiagonalMagneticMaterial(MaterialId& NewMaterialId,
								const double& mu_r_x, const double& mu_r_y, const double& mu_r_z,
								const double& sigma_h_x, const double& sigma_h_y, const double& sigma_h_z);
void AddIsotropicMaterial(MaterialId& NewMaterialId, const double& epsilon_r, const double& mu_r, const double& sigma_e, const double& sigma_h);
void AddDiagonalMaterial(MaterialId& NewMaterialId,
							const double& epsilon_r_x, const double& epsilon_r_y, const double& epsilon_r_z,
							const double& mu_r_x, const double& mu_r_y, const double& mu_r_z,
							const double& sigma_e_x, const double& sigma_e_y, const double& sigma_e_z,
							const double& sigma_h_x, const double& sigma_h_y, const double& sigma_h_z);

#endif // MATERIAL_H

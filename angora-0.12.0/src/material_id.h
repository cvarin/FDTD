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

#ifndef MATERIAL_ID_H
#define MATERIAL_ID_H

// For the material pointer types, short int (16-bit) is chosen for better range. (0-65,535)
// The memory difference between choosing char and short for both for E and H is ~11%
// The memory difference between choosing char and short for H is ~5%

//type for material indices with different (eps_x,cond_e_x) combinations in the x direction
typedef unsigned short ElectricMaterialIndexType_X;
//type for material indices with different (eps_y,cond_e_y) combinations in the y direction
typedef unsigned short ElectricMaterialIndexType_Y;
//type for material indices with different (eps_z,cond_e_z) combinations in the z direction
typedef unsigned short ElectricMaterialIndexType_Z;
//type for material indices with different (mu_x,cond_h_x) combinations in the x direction
typedef unsigned short MagneticMaterialIndexType_X;
//type for material indices with different (mu_y,cond_h_y) combinations in the y direction
typedef unsigned short MagneticMaterialIndexType_Y;
//type for material indices with different (mu_z,cond_h_z) combinations in the z direction
typedef unsigned short MagneticMaterialIndexType_Z;

//the maximum material indices representable by the material index types
#define FDTD_E_X_MAXINDEX numeric_limits<ElectricMaterialIndexType_X>::max()
#define FDTD_E_Y_MAXINDEX numeric_limits<ElectricMaterialIndexType_Y>::max()
#define FDTD_E_Z_MAXINDEX numeric_limits<ElectricMaterialIndexType_Z>::max()
#define FDTD_H_X_MAXINDEX numeric_limits<MagneticMaterialIndexType_X>::max()
#define FDTD_H_Y_MAXINDEX numeric_limits<MagneticMaterialIndexType_Y>::max()
#define FDTD_H_Z_MAXINDEX numeric_limits<MagneticMaterialIndexType_Z>::max()

//struct that represents a material with a distinct set of electric and magnetic properties in the x,y,z directions
struct MaterialId {ElectricMaterialIndexType_X ElectricIndex_X;
				 ElectricMaterialIndexType_Y ElectricIndex_Y;
				 ElectricMaterialIndexType_Z ElectricIndex_Z;
				 MagneticMaterialIndexType_X MagneticIndex_X;
				 MagneticMaterialIndexType_Y MagneticIndex_Y;
				 MagneticMaterialIndexType_Z MagneticIndex_Z;};

#endif // MATERIAL_ID_H

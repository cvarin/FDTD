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

#ifndef PLACEGEOM_H
#define PLACEGEOM_H

//for the definition of MaterialId
#include "material_id.h"


//void PlaceSlab(const string& tag, const double& start_position, const double& end_position); //Places dielectric slab of material with string tag
void PlaceSlab(const MaterialId& mat_id, const double& start_position, const double& end_position); //Places dielectric slab of material with given indices
void PlacePECSlab(const double& start_position, const double& end_position); //Places PEC slab
void PlaceGround(const int& GroundCoord);//Places PEC ground plane
void PlacePECBlock(//Places a square PEC block
			const int& BlockBack, const int& BlockFront,
			const int& BlockLeft, const int& BlockRight,
			const int& BlockLower, const int& BlockUpper);
void PlaceMaterialBlock(//Places a square block made of material with MaterialIndex
			const MaterialId& MaterialIdentifier, const int& BlockBack, const int& BlockFront,
			const int& BlockLeft, const int& BlockRight,
			const int& BlockLower, const int& BlockUpper);
void PlaceRotatedBlock(//Places a rotated square PEC block
			const double& CenterX, const double& CenterY, const double& CenterZ,
			const double& ThicknessX, const double& ThicknessY, const double& ThicknessZ,
			const double& phi_deg);
void PlaceCircularBlock(//Places a circular (in the xy-plane) PEC block
			const double& CenterX, const double& CenterY, const double& CenterZ,
			const double& Radius,
			const double& Height);
void PlaceSphere(//Places a PEC sphere
			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
			const double& Radius);
void PlaceMaterialSphere(//Places a sphere made of material with MaterialIndex
			const MaterialId& MaterialIdentifier,
			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
			const double& Radius);
void PlacePECMaskFromFile(//Reads a rectangular-prism-shaped boolean PEC mask from file and applies to grid
			const string& PECMaskFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor);
void PlaceSurfaceProfileFromFile(//Reads a rectangular homogeneous region with surface roughness from a file
			const string& SurfaceProfileFileName,
			const MaterialId& MaterialIdentifier,
			const int& xPos, const int& yPos, const int& zPos, const string& anchor);
void PlaceSurfaceEngravingProfileFromFile(//Reads a rectangular array specifying an engraving profile
			const string& SurfaceEngravingProfileFileName,
			const MaterialId& MaterialIdentifier,
			const int& xPos, const int& yPos, const int& zPos, const string& anchor);

/** these functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/
void analyze_layering();
void find_extremal_constitutive_params();
/** these functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/

#endif

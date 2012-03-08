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

#include "headers.h"

#include "place_obj.h"

#include "shape/Cshape.h"

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int FullIntPosMin_x,FullIntPosMax_x;
extern int FullIntPosMin_y,FullIntPosMax_y;
extern int FullIntPosMin_z,FullIntPosMax_z;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;


void place_obj(const MaterialId& mat_id, const Cshape* shapeptr)
{
	int i,j,k;

//	double BlockBack,BlockFront,BlockLeft,BlockRight,BlockLower,BlockUpper;

//	double BlockBack = shapeptr->bounding_box.back_limit;
//	double BlockFront = shapeptr->bounding_box.front_limit;
//	double BlockLeft = shapeptr->bounding_box.left_limit;
//	double BlockRight = shapeptr->bounding_box.right_limit;
//	double BlockLower = shapeptr->bounding_box.lower_limit;
//	double BlockUpper = shapeptr->bounding_box.upper_limit;

	//get the bounding box limits for quicker placement
	int BlockBack = shapeptr->bounding_box_back_cell();
	int BlockFront = shapeptr->bounding_box_front_cell();
	int BlockLeft = shapeptr->bounding_box_left_cell();
	int BlockRight = shapeptr->bounding_box_right_cell();
	int BlockLower = shapeptr->bounding_box_lower_cell();
	int BlockUpper = shapeptr->bounding_box_upper_cell();

	//place the bulk of the object
	//electric properties
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-1))
				{
					Media_Ex(i,j,k)=mat_id.ElectricIndex_X;
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-1))
				{
					Media_Ey(i,j,k)=mat_id.ElectricIndex_Y;
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-1,k-0.5))
				{
					Media_Ez(i,j,k)=mat_id.ElectricIndex_Z;
				}
			}
		}
	}
	//magnetic properties
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-0.5))
				{
					Media_Hx(i,j,k)=mat_id.MagneticIndex_X;
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-0.5))
				{
					Media_Hy(i,j,k)=mat_id.MagneticIndex_Y;
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-0.5,k-1))
				{
					Media_Hz(i,j,k)=mat_id.MagneticIndex_Z;
				}
			}
		}
	}

	/***************************************/
	/** TODO: Interpolation at boundaries **/
	/***************************************/
}

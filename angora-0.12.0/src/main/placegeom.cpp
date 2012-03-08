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

//Defines functions that place materials with different geometries in the grid.

#include "headers.h"

#include "placegeom.h"

#include "material.h"

#include <fstream>

extern double dx;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int rank;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern const int PEC;
extern const int vacuum;
extern Array<ElectricMaterialIndexType_X,1> Layering_e_x;
extern Array<ElectricMaterialIndexType_Y,1> Layering_e_y;
extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
extern Array<MagneticMaterialIndexType_X,1> Layering_h_x;
extern Array<MagneticMaterialIndexType_Y,1> Layering_h_y;
extern Array<MagneticMaterialIndexType_Z,1> Layering_h_z;

extern int number_of_layers;
extern Array<ElectricMaterialIndexType_X,1> LayerElectricMaterial_X;
extern Array<ElectricMaterialIndexType_Y,1> LayerElectricMaterial_Y;
extern Array<ElectricMaterialIndexType_Z,1> LayerElectricMaterial_Z;
extern Array<MagneticMaterialIndexType_X,1> LayerMagneticMaterial_X;
extern Array<MagneticMaterialIndexType_Y,1> LayerMagneticMaterial_Y;
extern Array<MagneticMaterialIndexType_Z,1> LayerMagneticMaterial_Z;
extern Array<int,1> LayerLowerZIndices;
extern Array<int,1> LayerThicknesses;
extern Array<bool,1> IsLayerGrounded;

extern Array<double,1> eps_x,eps_y,eps_z;
extern Array<double,1> mu_x,mu_y,mu_z;
extern Array<double,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<double,1> cond_h_x,cond_h_y,cond_h_z;

/** should be removed eventually **/
extern double epsilon_r_max_x,epsilon_r_min_x,epsilon_r_max_y,epsilon_r_min_y,epsilon_r_max_z,epsilon_r_min_z;
extern double mu_r_max_x,mu_r_min_x,mu_r_max_y,mu_r_min_y,mu_r_max_z,mu_r_min_z;
extern double epsilon_r_max,epsilon_r_min;
extern double mu_r_max,mu_r_min;
extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;
extern double c_upper,c_lower;
bool MaterialIsSame(const ElectricMaterialIndexType_X& Mat1_e_X,
						   const ElectricMaterialIndexType_Y& Mat1_e_Y,
						   const ElectricMaterialIndexType_Z& Mat1_e_Z,
						   const MagneticMaterialIndexType_X& Mat1_h_X,
						   const MagneticMaterialIndexType_Y& Mat1_h_Y,
						   const MagneticMaterialIndexType_Z& Mat1_h_Z,
						   const ElectricMaterialIndexType_X& Mat2_e_X,
						   const ElectricMaterialIndexType_Y& Mat2_e_Y,
						   const ElectricMaterialIndexType_Z& Mat2_e_Z,
						   const MagneticMaterialIndexType_X& Mat2_h_X,
						   const MagneticMaterialIndexType_Y& Mat2_h_Y,
						   const MagneticMaterialIndexType_Z& Mat2_h_Z);
/** should be removed eventually **/


//void PlaceSlab(const string& tag, const double& start_position, const double& end_position)
////Places an infinite dielectric slab (specified by the named material 'tag') at the specified position
////(For library use, not for config-file use)
//{
//	MaterialId mat_id;
//	if (!lookupMaterialWithTag(tag,mat_id))
//	{//string tag does not correspond to any material, throw an exception (DOCUMENT THIS!!)
//		throw TagNotFoundException();
//	}
//	else
//	{
//		PlaceSlab(mat_id,start_position,end_position); //use material that corresponds to the string tag
//	}
//}

void PlaceSlab(const MaterialId& mat_id, const double& start_position, const double& end_position)
//Places an infinite material slab at the specified position
//SlabLower is the z coordinate (w.r.t. the origin, in grid cells) of the lower interface
//SlabUpper is the z coordinate (w.r.t. the origin, in grid cells) of the upper interface
//Both SlabLower and SlabLower can be non-integer. Effective permittivities and permeabilities are used at interfaces for second order accuracy.
// The formulas for the permittivity and permeability are from
// K.-P. Hwang and A. C. Cangellaris, “Effective permittivities for second-order accurate FDTD equations at dielectric interfaces,” IEEE Microw. Wireless Compon. Lett., vol. 11, no. 4, pp. 158–60, Apr. 2001.
// The arithmetic averaging of the electric and magnetic conductivities is more or less empirical.
{
	int k; //counter for looping through z positions

	//Define the "uppermost" and "lowermost" grid voxel indices
	int SlabLower = (int)round(start_position)+1;//cell index is 1 more than the z coordinate of the x-y interface layer
	int SlabUpper = (int)round(end_position);//cell index is equal to the z coordinate of the x-y interface layer

	/** These are not used yet!!! start_position, end_position assumed integer! **/
	//the distance (in cell units) of the interfaces to the closest x-y grid interface layer
	double d_lower = start_position-SlabLower;
	double d_upper = end_position-SlabUpper;
	/** These are not used yet!!! start_position, end_position assumed integer! **/

	if (SlabUpper<SlabLower)
	{
		if (rank==0)
		{
			cout << "Error: Uppermost cell position (" << SlabUpper << ") in material slab is smaller than the lowermost cell position (" << SlabLower << ")." << endl;
		}
		exit(-1);
	}

	/** Place the bulk of the slab **/
	//Update eps_x,eps_y,mu_z
	for (k=max(klower,SlabLower); k<=min(kupper+1,SlabUpper+1); k++)
	{
		//x component of the E-field
		Media_Ex(Range::all(),Range::all(),k)=mat_id.ElectricIndex_X;
		//y component of the E-field
		Media_Ey(Range::all(),Range::all(),k)=mat_id.ElectricIndex_Y;
		//z component of the H-field
		Media_Hz(Range::all(),Range::all(),k)=mat_id.MagneticIndex_Z;
	}
	//Update eps_z,mu_x,mu_y
	for (k=max(klower,SlabLower); k<=min(kupper,SlabUpper); k++)
	{
		//z component of the E-field
		Media_Ez(Range::all(),Range::all(),k)=mat_id.ElectricIndex_Z;
		//x component of the H-field
		Media_Hx(Range::all(),Range::all(),k)=mat_id.MagneticIndex_X;
		//y component of the H-field
		Media_Hy(Range::all(),Range::all(),k)=mat_id.MagneticIndex_Y;
	}

	/** Place interface layers **/
	//create materials for the two planar interfaces with averaged constitutive parameters
	//Material indices for upper and lower interfaces
	MaterialId upper_interf_mat_id,lower_interf_mat_id;
	/** Place upper interface layer **/
	if ((SlabUpper<NCELLS_Z+2*NPML)&&(SlabUpper>=0))	//if the uppermost cell is outside the upper edge of the grid
														//then this is a half space, so do not place upper interface
														// (also check if the uppermost cell is beyond the lower edge)
	{
		double upper_interface_eps_x,upper_interface_eps_y,upper_interface_mu_z;
		double upper_interface_cond_e_x,upper_interface_cond_e_y,upper_interface_cond_h_z;
		//second permittivities in the permittivity averages are taken from the next higher z position
		upper_interface_eps_x = (eps_x(mat_id.ElectricIndex_X) + eps_x(Layering_e_x(SlabUpper+2)))/2;
		upper_interface_eps_y = (eps_y(mat_id.ElectricIndex_Y) + eps_y(Layering_e_y(SlabUpper+2)))/2;
		//second permeability in both the denominator and numerator below are taken from the next higher z position
		upper_interface_mu_z = 2*mu_z(mat_id.MagneticIndex_Z)*mu_z(Layering_h_z(SlabUpper+2))/(mu_z(mat_id.MagneticIndex_Z)+mu_z(Layering_h_z(SlabUpper+2)));
		//the following conductivity averages are empirical
		upper_interface_cond_e_x = (cond_e_x(mat_id.ElectricIndex_X) + cond_e_x(Layering_e_x(SlabUpper+2)))/2; //this is empirical
		upper_interface_cond_e_y = (cond_e_y(mat_id.ElectricIndex_Y) + cond_e_y(Layering_e_y(SlabUpper+2)))/2; //this is empirical
		upper_interface_cond_h_z = (cond_h_z(mat_id.MagneticIndex_Z) + cond_h_z(Layering_h_z(SlabUpper+2)))/2; //this is empirical

		// create the material using the averaged constitutive parameters
		AddDiagonalMaterial(upper_interf_mat_id,upper_interface_eps_x,upper_interface_eps_y,1,
												1,1,upper_interface_mu_z,
												upper_interface_cond_e_x,upper_interface_cond_e_y,0,
												0,0,upper_interface_cond_h_z);
		//Then, place the material index pointers
		if ((klower<=SlabUpper+1)&&(kupper>=SlabUpper)) //don't do anything if this interface layer does not belong to this node
		{
			//x component of the E-field
			Media_Ex(Range::all(),Range::all(),SlabUpper+1)=upper_interf_mat_id.ElectricIndex_X;
			//y component of the E-field
			Media_Ey(Range::all(),Range::all(),SlabUpper+1)=upper_interf_mat_id.ElectricIndex_Y;
			//z component of the H-field
			Media_Hz(Range::all(),Range::all(),SlabUpper+1)=upper_interf_mat_id.MagneticIndex_Z;
		}
	}

	/** Place lower interface layer **/
	if ((SlabLower>1)&&(SlabLower<=NCELLS_Z+2*NPML+1))	//if the lowermost cell is outside the lower edge of the grid
														//then this is a half space, so do not place lower interface
														// (also check if the lowermost cell is beyond the upper edge)
	{
		double lower_interface_eps_x,lower_interface_eps_y,lower_interface_mu_z;
		double lower_interface_cond_e_x,lower_interface_cond_e_y,lower_interface_cond_h_z;
		//second permittivities in the permittivity averages are taken from the next lower z position
		lower_interface_eps_x = (eps_x(mat_id.ElectricIndex_X) + eps_x(Layering_e_x(SlabLower-1)))/2;
		lower_interface_eps_y = (eps_y(mat_id.ElectricIndex_Y) + eps_y(Layering_e_y(SlabLower-1)))/2;
		//second permeability in both the denominator and numerator below are taken from the next lower z position
		lower_interface_mu_z = 2*mu_z(mat_id.MagneticIndex_Z)*mu_z(Layering_h_z(SlabLower-1))/(mu_z(mat_id.MagneticIndex_Z)+mu_z(Layering_h_z(SlabLower-1)));
		//the following conductivity averages are empirical
		lower_interface_cond_e_x = (cond_e_x(mat_id.ElectricIndex_X) + cond_e_x(Layering_e_x(SlabLower-1)))/2;
		lower_interface_cond_e_y = (cond_e_y(mat_id.ElectricIndex_Y) + cond_e_y(Layering_e_y(SlabLower-1)))/2;
		lower_interface_cond_h_z = (cond_h_z(mat_id.MagneticIndex_Z) + cond_h_z(Layering_h_z(SlabLower-1)))/2;

		// create the material using the averaged constitutive parameters
		AddDiagonalMaterial(lower_interf_mat_id,lower_interface_eps_x,lower_interface_eps_y,1,
												1,1,lower_interface_mu_z,
												lower_interface_cond_e_x,lower_interface_cond_e_y,0,
												0,0,lower_interface_cond_h_z);
		//Then, place the material index pointers
		if ((klower<=SlabLower)&&(kupper>=SlabLower-1)) //don't do anything if this interface layer does not belong to this node
		{
			//x component of the E-field
			Media_Ex(Range::all(),Range::all(),SlabLower)=lower_interf_mat_id.ElectricIndex_X;
			//y component of the E-field
			Media_Ey(Range::all(),Range::all(),SlabLower)=lower_interf_mat_id.ElectricIndex_Y;
			//z component of the H-field
			Media_Hz(Range::all(),Range::all(),SlabLower)=lower_interf_mat_id.MagneticIndex_Z;
		}
	}

	/** Update layering info **/
	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
	//x and y components
	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML+1,SlabUpper+1); k++)
	{
		Layering_e_x(k)=mat_id.ElectricIndex_X;
		Layering_e_y(k)=mat_id.ElectricIndex_Y;
		Layering_h_z(k)=mat_id.MagneticIndex_Z;
	}
	//z component
	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML,SlabUpper); k++)
	{
		Layering_e_z(k)=mat_id.ElectricIndex_Z;
		Layering_h_x(k)=mat_id.MagneticIndex_X;
		Layering_h_y(k)=mat_id.MagneticIndex_Y;
	}
	if ((SlabLower>1)&&(SlabLower<=NCELLS_Z+2*NPML+1))
	{
		Layering_e_x(SlabLower)=lower_interf_mat_id.ElectricIndex_X;
		Layering_e_y(SlabLower)=lower_interf_mat_id.ElectricIndex_Y;
		Layering_h_z(SlabLower)=lower_interf_mat_id.MagneticIndex_Z;
	}
	if ((SlabUpper<NCELLS_Z+2*NPML)&&(SlabUpper>=0))
	{
		Layering_e_x(SlabUpper+1)=upper_interf_mat_id.ElectricIndex_X;
		Layering_e_y(SlabUpper+1)=upper_interf_mat_id.ElectricIndex_Y;
		Layering_h_z(SlabUpper+1)=upper_interf_mat_id.MagneticIndex_Z;
	}
}

void PlacePECSlab(const double& start_position, const double& end_position)
//Places an infinite perfect-electric-conductor (PEC) slab at the specified position
{
	int k; //counter for looping through z positions

	//Define the "uppermost" and "lowermost" grid voxel indices
	int SlabLower = (int)round(start_position)+1;//cell index is 1 more than the z coordinate of the x-y interface layer
	int SlabUpper = (int)round(end_position);//cell index is equal to the z coordinate of the x-y interface layer

	/** These are not used yet!!! start_position, end_position assumed integer! **/
	//the distance (in cell units) of the interfaces to the closest x-y grid interface layer
	double d_lower = start_position-SlabLower;
	double d_upper = end_position-SlabUpper;
	/** These are not used yet!!! start_position, end_position assumed integer! **/

	//Update eps_x,eps_y
	for (k=max(klower,SlabLower); k<=min(kupper+1,SlabUpper+1); k++)
	{
		//x component
		Media_Ex(Range::all(),Range::all(),k)=PEC;
		//y component
		Media_Ey(Range::all(),Range::all(),k)=PEC;
	}
	//Update eps_z
	for (k=max(klower,SlabLower); k<=min(kupper,SlabUpper); k++)
	{
		//z component
		Media_Ez(Range::all(),Range::all(),k)=PEC;
	}

	//Update layering info
	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
	//x and y components
	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML+1,SlabUpper+1); k++)
	{
		Layering_e_x(k)=PEC;
		Layering_e_y(k)=PEC;
	}
	//z component
	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML,SlabUpper); k++)
	{
		Layering_e_z(k)=PEC;
	}

	//In PlaceSlab(), epsilon_r_upper, c_upper, etc. is updated at this point. (04/24/2011: Not anymore!) If the PEC slab extends to the upper or lower limit of the grid, it is the responsibility of the particular class that normally requires epsilon_r_upper, c_upper, etc. to be smart and accomodate the PEC layer at the top or bottom instead of using epsilon_r_upper, c_upper, etc.
}

void PlaceGround(const int& GroundCoord)
//Places an infinite PEC ground plane at the specified position
{
	//Define the "uppermost" and "lowermost" grid voxel indices
	int GroundCellIndex = GroundCoord+1;  //index of the cell right above the ground layer
//	//Define the "uppermost" and "lowermost" grid voxel indices
//	int GroundCellIndex = (int)round(GroundCoord)+1;  //index of the cell right above the ground layer

//	/** This is not used yet!!! GroundCoord assumed integer! **/
//	//the distance (in cell units) of the PEC ground layer to the next higher closest x-y grid grid interface
//	double d = GroundCoord-GroundCellIndex;
//	/** This is not used yet!!! GroundCoord assumed integer! **/

	//Place PEC ground
	if ((klower<=GroundCellIndex)&&(kupper>=GroundCellIndex-1))
	{
		//x component
		Media_Ex(Range::all(),Range::all(),GroundCellIndex)=PEC;
		//y component
		Media_Ey(Range::all(),Range::all(),GroundCellIndex)=PEC;
	}

	//Update layering info
	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
	if ((GroundCellIndex>=1)&&(GroundCellIndex<=NCELLS_Z+2*NPML+1))
	{
		Layering_e_x(GroundCellIndex)=PEC;
		Layering_e_y(GroundCellIndex)=PEC;
	}
}

void PlacePECBlock(
		const int& BlockBack, const int& BlockFront,	//indices of the rearmost and foremost cells in the block
		const int& BlockLeft, const int& BlockRight,	//indices of the leftmost and rightmost cells in the block
		const int& BlockLower, const int& BlockUpper)	//indices of the lowermost and uppermost cells in the block
//Places a PEC block at the specified position
{
	int i,j,k;
	//x component
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				Media_Ex(i,j,k)=PEC;
			}
		}
	}
	//y component
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				Media_Ey(i,j,k)=PEC;
			}
		}
	}
	//z component
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				Media_Ez(i,j,k)=PEC;
			}
		}
	}
}

void PlaceMaterialBlock(
			const MaterialId& MaterialIdentifier,
			const int& BlockBack, const int& BlockFront,	//indices of the rearmost and foremost cells in the block
			const int& BlockLeft, const int& BlockRight, 	//indices of the leftmost and rightmost cells in the block
			const int& BlockLower, const int& BlockUpper)	//indices of the lowermost and uppermost cells in the block
//Places a square block made of material with MaterialIndex at the specified position
{
	int i,j,k;
	//electric-field components
	//x component
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				Media_Ex(i,j,k)=MaterialIdentifier.ElectricIndex_X;
			}
		}
	}
	//y component
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				Media_Ey(i,j,k)=MaterialIdentifier.ElectricIndex_Y;
			}
		}
	}
	//z component
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				Media_Ez(i,j,k)=MaterialIdentifier.ElectricIndex_Z;
			}
		}
	}

	//magnetic-field components
	//x component
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				Media_Hx(i,j,k)=MaterialIdentifier.MagneticIndex_X;
			}
		}
	}
	//y component
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				Media_Hy(i,j,k)=MaterialIdentifier.MagneticIndex_Y;
			}
		}
	}
	//z component
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				Media_Hz(i,j,k)=MaterialIdentifier.MagneticIndex_Z;
			}
		}
	}
}

void PlaceRotatedBlock(
			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
			const double& ThicknessX, const double& ThicknessY, const double& ThicknessZ,	//thicknesses
			const double& phi_deg)			//rotation angle (in degrees)
//Places a rotated (in the xy-plane) PEC block at the specified position
{
	double x,y,z;	//x,y,z coordinates with respect to the center of the block
	double x_dist,y_dist,z_dist;	//distances of the point to the orthogonal faces of the cube

	double phi = phi_deg*M_PI/180;

	int i,j,k;
	//x component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i+0.5-CenterX);
				y = (j-CenterY);
				z = (k-CenterZ);
				x_dist = x*cos(phi)+y*sin(phi);
				y_dist = -x*sin(phi)+y*cos(phi);
				z_dist = z;
				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
				{
					Media_Ex(i,j,k)=PEC;
				}
			}
		}
	}
	//y component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i-CenterX);
				y = (j+0.5-CenterY);
				z = (k-CenterZ);
				x_dist = x*cos(phi)+y*sin(phi);
				y_dist = -x*sin(phi)+y*cos(phi);
				z_dist = z;
				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
				{
					Media_Ey(i,j,k)=PEC;
				}
			}
		}
	}
	//z component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i-CenterX);
				y = (j-CenterY);
				z = (k+0.5-CenterZ);
				x_dist = x*cos(phi)+y*sin(phi);
				y_dist = -x*sin(phi)+y*cos(phi);
				z_dist = z;
				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
				{
					Media_Ez(i,j,k)=PEC;
				}
			}
		}
	}
}

void PlaceCircularBlock(
			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
			const double& Radius, //radius of the bowtie (in cells)
			const double& Height)	//height of the block (in cells)
//Places a circular (in the xy-plane) PEC block
{
	double x,y,z;	//x,y,z coordinates with respect to the center of the block
	double dist_r;	//polar distance of the point to the center (in the xy plane)

	int i,j,k;
	//x component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i+0.5-CenterX);
				y = (j-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2));
				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
				{
					Media_Ex(i,j,k)=PEC;
				}
			}
		}
	}
	//y component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i-CenterX);
				y = (j+0.5-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2));
				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
				{
					Media_Ey(i,j,k)=PEC;
				}
			}
		}
	}
	//z component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i-CenterX);
				y = (j-CenterY);
				z = (k+0.5-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2));
				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
				{
					Media_Ez(i,j,k)=PEC;
				}
			}
		}
	}
}

void PlaceSphere(const double& CenterX, const double& CenterY, const double& CenterZ, const double& Radius)
//Places a PEC sphere
{
	double x,y,z;	//x,y,z coordinates with respect to the center of the sphere
	double dist_r;	//spherical distance of the point to the center of the sphere

	int i,j,k;
	//x component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i+0.5-CenterX);
				y = (j-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ex(i,j,k)=PEC;
				}
			}
		}
	}
	//y component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i-CenterX);
				y = (j+0.5-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ey(i,j,k)=PEC;
				}
			}
		}
	}
	//z component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i-CenterX);
				y = (j-CenterY);
				z = (k+0.5-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ez(i,j,k)=PEC;
				}
			}
		}
	}
}

void PlaceMaterialSphere(
			const MaterialId& MaterialIdentifier,
			const double& CenterX, const double& CenterY, const double& CenterZ,
			const double& Radius)
{//Places a sphere made of material with identifier MaterialIdentifier
	double x,y,z;	//x,y,z coordinates with respect to the center of the sphere
	double dist_r;	//spherical distance of the point to the center of the sphere

	int i,j,k;
	//electric-field components
	//x component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i+0.5-CenterX);
				y = (j-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ex(i,j,k)=MaterialIdentifier.ElectricIndex_X;
				}
			}
		}
	}
	//y component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i-CenterX);
				y = (j+0.5-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ey(i,j,k)=MaterialIdentifier.ElectricIndex_Y;
				}
			}
		}
	}
	//z component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i-CenterX);
				y = (j-CenterY);
				z = (k+0.5-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Ez(i,j,k)=MaterialIdentifier.ElectricIndex_Z;
				}
			}
		}
	}

	//magnetic-field components
	//x component
	for (i=iback; i<=ifront+1; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i-CenterX);
				y = (j+0.5-CenterY);
				z = (k+0.5-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Hx(i,j,k)=MaterialIdentifier.MagneticIndex_X;
				}
			}
		}
	}
	//y component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright+1; j++)
		{
			for (k=klower; k<=kupper; k++)
			{
				x = (i+0.5-CenterX);
				y = (j-CenterY);
				z = (k+0.5-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Hy(i,j,k)=MaterialIdentifier.MagneticIndex_Y;
				}
			}
		}
	}
	//z component
	for (i=iback; i<=ifront; i++)
	{
		for (j=jleft; j<=jright; j++)
		{
			for (k=klower; k<=kupper+1; k++)
			{
				x = (i+0.5-CenterX);
				y = (j+0.5-CenterY);
				z = (k-CenterZ);
				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
				if (dist_r<=Radius)
				{
					Media_Hz(i,j,k)=MaterialIdentifier.MagneticIndex_Z;
				}
			}
		}
	}
}

void PlacePECMaskFromFile(const string& PECMaskFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor)
{//Reads a rectangular-prism-shaped boolean PEC mask from file and applies to grid
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the applied mask (measured in cells from the back-left-lower corner of the grid)
	if (rank==0)
	{
		cout << endl << "Reading PEC mask from " << PECMaskFileName << " ..." << endl;
	}
	ifstream PECMaskFile;	//temporary ifstream object for reading the input
	PECMaskFile.open(PECMaskFileName.c_str(),ios::binary);	//open file for reading
	if (!PECMaskFile)
	{
		cout << "Error opening PEC mask input file " << PECMaskFileName << "." << endl << endl;
		exit(-1);
	}
	int xExtent,yExtent,zExtent;	//x, y and z extents of the material region (in cells)
	PECMaskFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
	PECMaskFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
	PECMaskFile.read((char*)&zExtent,sizeof(zExtent));	//read the z extent

	//calculate the coordinates of the back-left-lower corner of the region
	int xCornerPos=xPos;
	int yCornerPos=yPos;
	int zCornerPos=zPos;
	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
	{
		if (rank==0)
		{
			cout << "Invalid anchor point \"" << anchor << "\" for PEC mask input file " << PECMaskFileName << " in node " << rank << endl << endl;
			exit(-1);
		}
	}

	if (anchor=="center")
	{
		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
	}
	else if (anchor=="BLL")	//back-left-lower
	{
		//reference point is the back-left-lower corner by default
	}
	else if (anchor=="BLU")	//back-left-upper
	{
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="BRL")	//back-right-lower
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="BRU")	//back-right-upper
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FLL")	//front-left-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
	}
	else if (anchor=="FLU")	//front-left-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FRL")	//front-right-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="FRU")	//front-right-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}

	//indices of the cell at the back-left-lower corner
	//these are simply equal to the coordinates of the back-left-lower corner plus one
	int xCornerCell = xCornerPos+1;
	int yCornerCell = yCornerPos+1;
	int zCornerCell = zCornerPos+1;

	//every node has to read the file, and update the necessary portions of their grid
	Array<bool,1> PEC_mask_row(Range(xCornerCell,xCornerCell+xExtent-1));	//temporary array of size xExtent that will hold the boolean values belonging to a x-row in the 3-D PEC mask array in the file
	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
	int i,j,k;
	for (k=zCornerCell; k<=zCornerCell+zExtent-1; k++)
	{
		for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
		{
			PECMaskFile.read((char*)PEC_mask_row.data(),xExtent*sizeof(PEC_mask_row(0)));	//read the permittivity x-row into temporary array
			//Media_Ex
			if ((k>=klower)&&(k<=kupper+1))
			{
				if ((j>=jleft)&&(j<=jright+1))
				{
					for (i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtent-1); i++)
					{
						if (PEC_mask_row(i)==true)
						{
							//mark this position as perfect electric conductor
							Media_Ex(i,j,k) = PEC;
						}
					}
				}
			}

			//Media_Ey
			if ((k>=klower)&&(k<=kupper+1))
			{
				if ((j>=jleft)&&(j<=jright))
				{
					for (i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
					{
						if (PEC_mask_row(i)==true)
						{
							//mark this position as perfect electric conductor
							Media_Ey(i,j,k) = PEC;
						}
					}
				}
			}

			//Media_Ez
			if ((k>=klower)&&(k<=kupper))
			{
				if ((j>=jleft)&&(j<=jright+1))
				{
					for (i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
					{
						if (PEC_mask_row(i)==true)
						{
							//mark this position as perfect electric conductor
							Media_Ez(i,j,k) = PEC;
						}
					}
				}
			}
		}
	}

	//finally, close the file
	PECMaskFile.close();
	if (rank==0)
	{
		cout << "PEC mask applied to grid." << endl << endl;
	}
}

void PlaceSurfaceProfileFromFile(
			const string& SurfaceProfileFileName,
			const MaterialId& MaterialIdentifier,
			const int& xPos, const int& yPos, const int& zPos, const string& anchor)
{
//Reads a rectangular homogeneous region with surface roughness from a file.
//The region is completely specified by a 2D array indicating the surface fluctuations, the average thickness of the region, and the relative dielectric permittivity of the region (equal to square of the refractive index).
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the region (measured in cells from the back-left-lower corner of the grid)
	if (rank==0)
	{
		cout << endl << "Reading surface profile from " << SurfaceProfileFileName << " ..." << endl;
	}
	ifstream SurfaceProfileFile;	//temporary ifstream object for reading the input
	SurfaceProfileFile.open(SurfaceProfileFileName.c_str(),ios::binary);	//open file for reading
	if (!SurfaceProfileFile)
	{
		cout << "Error opening surface profile input file " << SurfaceProfileFileName << "." << endl << endl;
		exit(-1);
	}
	int xExtent,yExtent,zExtent;	//x, y and average z extents of the material region (in cells)
	SurfaceProfileFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
	SurfaceProfileFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
	SurfaceProfileFile.read((char*)&zExtent,sizeof(zExtent));	//read the average z extent

	//Following removed after ver.>= 0.12
/*	double eps_r; //relative dielectric permittivity of the homogeneous region
	//read the relative dielectric permittivity of the homogeneous region
	SurfaceProfileFile.read((char*)&eps_r,sizeof(eps_r));*/
// 	//add the new material to the material-index array
// 	int surfacematerialindex = AddMaterial(eps_r,0);	//relative permittivity of the new material is eps_r

	//calculate the coordinates of the back-left-lower corner of the region
	int xCornerPos=xPos;
	int yCornerPos=yPos;
	int zCornerPos=zPos;
	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
	{
		if (rank==0)
		{
			cout << "Invalid anchor point \"" << anchor << "\" for surface-profile input file " << SurfaceProfileFileName << " in node " << rank << endl << endl;
			exit(-1);
		}
	}

	if (anchor=="center")
	{
		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
	}
	else if (anchor=="BLL")	//back-left-lower
	{
		//reference point is the back-left-lower corner by default
	}
	else if (anchor=="BLU")	//back-left-upper
	{
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="BRL")	//back-right-lower
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="BRU")	//back-right-upper
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FLL")	//front-left-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
	}
	else if (anchor=="FLU")	//front-left-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FRL")	//front-right-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="FRU")	//front-right-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}

	//indices of the cell at the back-left-lower corner
	//these are simply equal to the coordinates of the back-left-lower corner plus one
	int xCornerCell = xCornerPos+1;
	int yCornerCell = yCornerPos+1;
	int zCornerCell = zCornerPos+1;

	//extra and absolute thickness values at each point of the surface
	int extra_thickness,thickness; //in cells

	//every node has to read the file, and update the necessary portions of their grid
	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
	int i,j,k;
	for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
	{
		for (i=xCornerCell; i<=xCornerCell+xExtent-1; i++)
		{
			SurfaceProfileFile.read((char*)&extra_thickness,sizeof(extra_thickness));	//read the surface-profile at the (x,y) position
			//calculate the thickness at this (x,y) position by adding the profile value to the average thickness (zExtent)
			//(extra_thickness might also be negative)
			thickness = zExtent+extra_thickness;
			for (k=zCornerCell; k<=zCornerCell+thickness-1; k++)
			{
				//Media_Ex
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						if ((i>=iback)&&(i<=ifront))
						{
								Media_Ex(i,j,k) = MaterialIdentifier.ElectricIndex_X;
						}
					}
				}

				//Media_Ey
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright))
					{
						if ((i>=iback)&&(i<=ifront+1))
						{
								Media_Ey(i,j,k) = MaterialIdentifier.ElectricIndex_Y;
						}
					}
				}

				//Media_Ez
				if ((k>=klower)&&(k<=kupper))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						if ((i>=iback)&&(i<=ifront+1))
						{
								Media_Ez(i,j,k) = MaterialIdentifier.ElectricIndex_Z;
						}
					}
				}
			}
		}
	}

	//finally, close the file
	SurfaceProfileFile.close();
	if (rank==0)
	{
		cout << "Surface profile read." << endl << endl;
	}
}

void PlaceSurfaceEngravingProfileFromFile(
			const string& SurfaceEngravingProfileFileName,
			const MaterialId& MaterialIdentifier,
			const int& xPos, const int& yPos, const int& zPos, const string& anchor)
{
//Reads a 2D array that specifies an engraving profile. The surface is engraved at the specified height values in the array (in grid cells), and filled with the material specified by "MaterialIndex".
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the base of the engraving profile (measured in cells from the back-left-lower corner of the grid)
	if (rank==0)
	{
		cout << endl << "Reading surface-engraving profile from " << SurfaceEngravingProfileFileName << " ..." << endl;
	}
	ifstream SurfaceEngravingProfileFile;	//temporary ifstream object for reading the input
	SurfaceEngravingProfileFile.open(SurfaceEngravingProfileFileName.c_str(),ios::binary);	//open file for reading
	if (!SurfaceEngravingProfileFile)
	{
		cout << "Error opening surface-engraving profile input file " << SurfaceEngravingProfileFileName << "." << endl << endl;
		exit(-1);
	}
	int xExtent,yExtent,zExtent;	//x, y and average z extents of the material region (in cells)
	SurfaceEngravingProfileFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
	SurfaceEngravingProfileFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
	SurfaceEngravingProfileFile.read((char*)&zExtent,sizeof(zExtent));	//read the average z extent

/*	double eps_r; //relative dielectric permittivity of the homogeneous region
	//read the relative dielectric permittivity of the homogeneous region
	SurfaceEngravingProfileFile.read((char*)&eps_r,sizeof(eps_r));*/

// 	//add the new material to the material-index array
// 	int surfacematerialindex = AddMaterial(eps_r,0);	//relative permittivity of the new material is eps_r

	//calculate the coordinates of the back-left-lower corner of the region
	int xCornerPos=xPos;
	int yCornerPos=yPos;
	int zCornerPos=zPos;
	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
	{
		if (rank==0)
		{
			cout << "Invalid anchor point \"" << anchor << "\" for surface-engraving profile input file " << SurfaceEngravingProfileFileName << " in node " << rank << endl << endl;
			exit(-1);
		}
	}

	if (anchor=="center")
	{
		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
	}
	else if (anchor=="BLL")	//back-left-lower
	{
		//reference point is the back-left-lower corner by default
	}
	else if (anchor=="BLU")	//back-left-upper
	{
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="BRL")	//back-right-lower
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="BRU")	//back-right-upper
	{
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FLL")	//front-left-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
	}
	else if (anchor=="FLU")	//front-left-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}
	else if (anchor=="FRL")	//front-right-lower
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
	}
	else if (anchor=="FRU")	//front-right-upper
	{
		xCornerPos = xPos - xExtent;	//shift reference to the front corner
		yCornerPos = yPos - yExtent;	//shift reference to the right corner
		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
	}

	//indices of the cell at the back-left-lower corner
	//these are simply equal to the coordinates of the back-left-lower corner plus one
	int xCornerCell = xCornerPos+1;
	int yCornerCell = yCornerPos+1;
	int zCornerCell = zCornerPos+1;

	//extra and absolute thickness values at each point of the surface
	int extra_thickness,thickness; //in cells
// cout << xCornerCell << "," << yCornerCell << "," << zCornerCell << endl;
	//every node has to read the file, and update the necessary portions of their grid
	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
	int i,j,k;
	for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
	{
		for (i=xCornerCell; i<=xCornerCell+xExtent-1; i++)
		{
			SurfaceEngravingProfileFile.read((char*)&extra_thickness,sizeof(extra_thickness));	//read the surface-profile at the (x,y) position
			//calculate the thickness at this (x,y) position by adding the profile value to the average thickness (zExtent)
			//(extra_thickness might also be negative)
			thickness = zExtent+extra_thickness;
// 			for (int k=zCornerCell; k<=zCornerCell+thickness-1; k++)
// 			cout << zExtent << "," << extra_thickness << endl;
// 			cout << xCornerCell << "," << yCornerCell << "," << zCornerCell << endl;
// 			exit(-1);
// cout << zCornerCell+thickness << endl;
			for (k=zCornerCell+thickness; k<=kupper+1; k++)  //from bottom of the engraving profile to the top of the grid
			{
// if ((j>=yCornerCell)&&(j<=yCornerCell+1)&&(i>=xCornerCell)&&(i<=xCornerCell+1)&&(k==(zCornerCell+thickness))){
// if ((j>=yCornerCell)&&(j<=yCornerCell+20)&&(i>=xCornerCell)&&(i<=xCornerCell)&&(k==(zCornerCell+thickness))){
// cout << i << "," << j << "," << k << endl;
// }
				//Media_Ex
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						if ((i>=iback)&&(i<=ifront))
						{
								Media_Ex(i,j,k) = MaterialIdentifier.ElectricIndex_X; //engrave with the specified material
						}
					}
				}

				//Media_Ey
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright))
					{
						if ((i>=iback)&&(i<=ifront+1))
						{
								Media_Ey(i,j,k) = MaterialIdentifier.ElectricIndex_Y; //engrave with the specified material
						}
					}
				}

				//Media_Ez
				if ((k>=klower)&&(k<=kupper))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						if ((i>=iback)&&(i<=ifront+1))
						{
								Media_Ez(i,j,k) = MaterialIdentifier.ElectricIndex_Z; //engrave with the specified material
						}
					}
				}
			}
		}
	}

	//finally, close the file
	SurfaceEngravingProfileFile.close();
	if (rank==0)
	{
		cout << "Surface-engraving profile read." << endl << endl;
	}
}

/****************************************************************************************************************/
/****************************************************************************************************************/
/** following functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/
/****************************************************************************************************************/
/****************************************************************************************************************/

void analyze_layering()
//Analyzes the layering structure, and builds the layering information arrays
{
	/** NOTE: **/
	/** This routine can only analyze isotropic layers in its current form. **/
	/** When the classes that use this routine are generalized to handle anisotropic layers, this routine will also be modified. **/
	/** Furthermore, the routine can only handle permittivity variations. Permeability support will be added in the near future. **/

	//the layers are numerically ordered from 0 to 'number_of_layers-1' beginning from the bottom (k=1) of the grid
	//initialize the number of layers and the layering arrays
	number_of_layers=1;
	LayerElectricMaterial_X.resize(1);
	LayerElectricMaterial_Y.resize(1);
	LayerElectricMaterial_Z.resize(1);
	LayerMagneticMaterial_X.resize(1);
	LayerMagneticMaterial_Y.resize(1);
	LayerMagneticMaterial_Z.resize(1);
	LayerLowerZIndices.resize(1);
	LayerThicknesses.resize(1);
	IsLayerGrounded.resize(1);
	//first layer is made of whatever is at the bottom of the grid
	LayerElectricMaterial_X(0) = Layering_e_z(1); /** isotropic layer assumed! **/
	LayerElectricMaterial_Y(0) = Layering_e_z(1); /** isotropic layer assumed! **/
	LayerElectricMaterial_Z(0) = Layering_e_z(1); /** isotropic layer assumed! **/
	LayerMagneticMaterial_X(0) = Layering_h_x(1); /** isotropic layer assumed! **/
	LayerMagneticMaterial_Y(0) = Layering_h_x(1); /** isotropic layer assumed! **/
	LayerMagneticMaterial_Z(0) = Layering_h_x(1); /** isotropic layer assumed! **/
	//lowermost index of the first layer is naturally 1
	LayerLowerZIndices(0) = 1;
	//thickness of the first layer is initialized to 1, and updated inside the loop below
	LayerThicknesses(0) = 1;
	//first layer is not grounded
	IsLayerGrounded(0)=false;

	if (NCELLS_Z+2*NPML>=2)
	{//sort of a silly check, but...
	for (int k=2; k<=NCELLS_Z+2*NPML; k++)
	{
		//is the material at k different than the one at k-1?
		if (!MaterialIsSame(Layering_e_z(k),Layering_e_z(k),Layering_e_z(k),  //electric properties
							Layering_h_x(k),Layering_h_x(k),Layering_h_x(k),  //magnetic properties
							Layering_e_z(k-1),Layering_e_z(k-1),Layering_e_z(k-1),  //electric properties
							Layering_h_x(k-1),Layering_h_x(k-1),Layering_h_x(k-1)))  //magnetic properties
							/** isotropy assumed! **/
		{//new material is encountered
			number_of_layers++;
			//increase the array sizes by 1
			LayerElectricMaterial_X.resizeAndPreserve(number_of_layers);
			LayerElectricMaterial_Y.resizeAndPreserve(number_of_layers);
			LayerElectricMaterial_Z.resizeAndPreserve(number_of_layers);
			LayerMagneticMaterial_X.resizeAndPreserve(number_of_layers);
			LayerMagneticMaterial_Y.resizeAndPreserve(number_of_layers);
			LayerMagneticMaterial_Z.resizeAndPreserve(number_of_layers);
			LayerLowerZIndices.resizeAndPreserve(number_of_layers);
			LayerThicknesses.resizeAndPreserve(number_of_layers);
			IsLayerGrounded.resizeAndPreserve(number_of_layers);
			//fill the info for the new layer
			LayerElectricMaterial_X(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
			LayerElectricMaterial_Y(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
			LayerElectricMaterial_Z(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
			LayerMagneticMaterial_X(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
			LayerMagneticMaterial_Y(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
			LayerMagneticMaterial_Z(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
			LayerLowerZIndices(number_of_layers-1) = k;
			LayerThicknesses(number_of_layers-1) = 1; //initialize the thickness to 1
			//check if the new layer is grounded or not
			if (Layering_e_x(k)==PEC)  /** isotropy assumed! **/
			{
				IsLayerGrounded(number_of_layers-1) = true;
			}
			else
			{
				IsLayerGrounded(number_of_layers-1) = false;
			}
		}
		else
		{//no new material at this z position, but there can be a PEC sheet between k and k-1
			if (Layering_e_x(k)==PEC)  /** isotropy assumed! **/
			{//there is a PEC sheet between k and k-1
				if (Layering_e_z(k-1)!=PEC)  /** isotropy assumed! **/
				{//if the previous z position was not PEC, start a new grounded layer
					number_of_layers++;
					//increase the array sizes by 1
					LayerElectricMaterial_X.resizeAndPreserve(number_of_layers);
					LayerElectricMaterial_Y.resizeAndPreserve(number_of_layers);
					LayerElectricMaterial_Z.resizeAndPreserve(number_of_layers);
					LayerMagneticMaterial_X.resizeAndPreserve(number_of_layers);
					LayerMagneticMaterial_Y.resizeAndPreserve(number_of_layers);
					LayerMagneticMaterial_Z.resizeAndPreserve(number_of_layers);
					LayerLowerZIndices.resizeAndPreserve(number_of_layers);
					LayerThicknesses.resizeAndPreserve(number_of_layers);
					IsLayerGrounded.resizeAndPreserve(number_of_layers);
					//fill the info for the new layer
					LayerElectricMaterial_X(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
					LayerElectricMaterial_Y(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
					LayerElectricMaterial_Z(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
					LayerMagneticMaterial_X(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
					LayerMagneticMaterial_Y(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
					LayerMagneticMaterial_Z(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
					LayerLowerZIndices(number_of_layers-1) = k;
					LayerThicknesses(number_of_layers-1) = 1; //initialize the thickness to 1
					IsLayerGrounded(number_of_layers-1) = true;
				}
				else
				{//the previous z position was also PEC, so never mind...
					LayerThicknesses(number_of_layers-1) += 1; //increase the thickness of the current layer by 1
				}
			}
			else
			{//no new material at this z position, AND there is no PEC sheet between k and k-1, so we are still in the same layer
				LayerThicknesses(number_of_layers-1) += 1; //increase the thickness of the current layer by 1
			}
		}
	}
	}
// cout << "Layering_e_z is " << Layering_e_z << endl;
// cout << "Layering_h_x is " << Layering_h_x << endl;
// cout << "number_of_layers is " << number_of_layers << endl;
//if (rank==0)  cout << "LayerElectricMaterial_X is " << LayerElectricMaterial_X << endl;
//if (rank==0) cout << "LayerElectricMaterial_Y is " << LayerElectricMaterial_Y << endl;
//if (rank==0)  cout << "LayerElectricMaterial_Z is " << LayerElectricMaterial_Z << endl;
//if (rank==0)  cout << "LayerMagneticMaterial_X is " << LayerMagneticMaterial_X << endl;
//if (rank==0)  cout << "LayerMagneticMaterial_Y is " << LayerMagneticMaterial_Y << endl;
//if (rank==0)  cout << "LayerMagneticMaterial_Z is " << LayerMagneticMaterial_Z << endl;
// cout << "LayerLowerZIndices is " << LayerLowerZIndices << endl;
// cout << "LayerThicknesses is " << LayerThicknesses << endl;
// cout << "IsLayerGrounded is " << IsLayerGrounded << endl;
}

void find_extremal_constitutive_params()
{
	//if the slab extends to the upper limit of the grid, then update the relative permittivity and velocity in the uppermost layer
	/** NOTE: **/
	/** epsilon_r_upper, mu_r_upper, c_upper, epsilon_r_lower, mu_r_lower, c_lower are only defined for an isotropic medium!!**/
	/** Using the x component, as the y and z components are assumed to be the same **/
	/** This can be changed when the classes using these variables are generalized to handle anisotropic media **/
	epsilon_r_upper = eps_x(LayerElectricMaterial_X(number_of_layers-1));
	mu_r_upper = mu_x(LayerMagneticMaterial_X(number_of_layers-1));
	c_upper = c/sqrt(epsilon_r_upper*mu_r_upper);

	//if the slab extends to the lower limit of the grid, then update the relative permittivity and velocity in the lowermost layer
	epsilon_r_lower = eps_x(LayerElectricMaterial_X(0));
	mu_r_lower = mu_x(LayerMagneticMaterial_X(0));
	c_lower = c/sqrt(epsilon_r_lower*mu_r_lower);

	//Find the maximum and minimum permittivity/permeability values in the grid
	//max/min permittivities
	Array<double,1> eps_x_noPEC = eps_x(Range(1,eps_x.size()-1)); //exclude the PEC material
	Array<double,1> eps_y_noPEC = eps_y(Range(1,eps_x.size()-1)); //exclude the PEC material
	Array<double,1> eps_z_noPEC = eps_z(Range(1,eps_x.size()-1)); //exclude the PEC material
	epsilon_r_max_x = max(eps_x_noPEC); epsilon_r_min_x = min(eps_x_noPEC);
	epsilon_r_max_y = max(eps_y_noPEC); epsilon_r_min_y = min(eps_y_noPEC);
	epsilon_r_max_z = max(eps_z_noPEC); epsilon_r_min_z = min(eps_z_noPEC);
	epsilon_r_max = max(max(epsilon_r_max_x,epsilon_r_max_y),epsilon_r_max_z);
	epsilon_r_min = min(min(epsilon_r_min_x,epsilon_r_min_y),epsilon_r_min_z);
	//max/min permeabilities
	Array<double,1> mu_x_noPEC = mu_x(Range(1,mu_x.size()-1)); //exclude the PEC material
	Array<double,1> mu_y_noPEC = mu_y(Range(1,mu_x.size()-1)); //exclude the PEC material
	Array<double,1> mu_z_noPEC = mu_z(Range(1,mu_x.size()-1)); //exclude the PEC material
	mu_r_max_x = max(mu_x_noPEC); mu_r_min_x = min(mu_x_noPEC);
	mu_r_max_y = max(mu_y_noPEC); mu_r_min_y = min(mu_y_noPEC);
	mu_r_max_z = max(mu_z_noPEC); mu_r_min_z = min(mu_z_noPEC);
	mu_r_max = max(max(mu_r_max_x,mu_r_max_y),mu_r_max_z);
	mu_r_min = min(min(mu_r_min_x,mu_r_min_y),mu_r_min_z);
}
/****************************************************************************************************************/
/****************************************************************************************************************/
/** the above functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/
/****************************************************************************************************************/
/****************************************************************************************************************/

inline bool MaterialIsSame(const ElectricMaterialIndexType_X& Mat1_e_X,
						   const ElectricMaterialIndexType_Y& Mat1_e_Y,
						   const ElectricMaterialIndexType_Z& Mat1_e_Z,
						   const MagneticMaterialIndexType_X& Mat1_h_X,
						   const MagneticMaterialIndexType_Y& Mat1_h_Y,
						   const MagneticMaterialIndexType_Z& Mat1_h_Z,
						   const ElectricMaterialIndexType_X& Mat2_e_X,
						   const ElectricMaterialIndexType_Y& Mat2_e_Y,
						   const ElectricMaterialIndexType_Z& Mat2_e_Z,
						   const MagneticMaterialIndexType_X& Mat2_h_X,
						   const MagneticMaterialIndexType_Y& Mat2_h_Y,
						   const MagneticMaterialIndexType_Z& Mat2_h_Z)
{
	//simplest version, but ignores the possibility of two material indices representing same material properties
// 	return Mat1==Mat2;
	//more stringent version, which requires the matching of the material properties to machine precision (in double format)
	const int relaxation_factor = 100;
	return ((abs(eps_x(Mat1_e_X)-eps_x(Mat2_e_X))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(eps_y(Mat1_e_Y)-eps_y(Mat2_e_Y))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(eps_z(Mat1_e_Z)-eps_z(Mat2_e_Z))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_x(Mat1_h_X)-mu_x(Mat2_h_X))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_y(Mat1_h_Y)-mu_y(Mat2_h_Y))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_z(Mat1_h_Z)-mu_z(Mat2_h_Z))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_x(Mat1_e_X)-cond_e_x(Mat2_e_X))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_y(Mat1_e_Y)-cond_e_y(Mat2_e_Y))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_z(Mat1_e_Z)-cond_e_z(Mat2_e_Z))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_x(Mat1_h_X)-cond_h_x(Mat2_h_X))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_y(Mat1_h_Y)-cond_h_y(Mat2_h_Y))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_z(Mat1_h_Z)-cond_h_z(Mat2_h_Z))<LIBSTD_DBL_EPSILON*relaxation_factor));
}
/****************************************************************************************************************/

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

#ifndef MATFILE_H
#define MATFILE_H

//template declaration for the function that reads a material region from a file

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

#include "material.h"

#include <fstream>

#ifndef FDTD_MAX_NEWMAT
#define FDTD_MAX_NEWMAT 1000
#endif

extern int rank;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern int NumOfElectricMaterials_X;
extern int NumOfElectricMaterials_Y;
extern int NumOfElectricMaterials_Z;
extern int NumOfMagneticMaterials_X;
extern int NumOfMagneticMaterials_Y;
extern int NumOfMagneticMaterials_Z;


template<typename MatType> //data type for the material property
void PlaceMaterialRegionFromFile(const string& MaterialFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor, const string& constitutive_param_type, const int& max_number_of_new_materials = FDTD_MAX_NEWMAT)
{//Reads rectangular-prism-shaped dielectric region from file and places into grid
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the region (measured in cells from the back-left-lower corner of the grid)

	if ((constitutive_param_type!="rel_permittivity")&&(constitutive_param_type!="rel_permeability")&&(constitutive_param_type!="electric_conductivity"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,constitutive_param_type,
			"(valid arguments are \"rel_permittivity\", \"rel_permeability\", or \"electric_conductivity\")");
	}

	ifstream MaterialFile;	//temporary ifstream object for reading the input
	MaterialFile.open(MaterialFileName.c_str(),ios::binary);	//open file for reading
	if (!MaterialFile)
	{
		cout << "Error opening material input file " << MaterialFileName << "." << endl << endl;
		exit(-1);
	}
	int xExtent,yExtent,zExtent;	//x, y and z extents of the material region (in cells)
	MaterialFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
	MaterialFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
	MaterialFile.read((char*)&zExtent,sizeof(zExtent));	//read the z extent

	//read through the file to determine the maximum and minumum constitutive parameter values
	int pos_saved = MaterialFile.tellg();	//first, save the current position of the read pointer
	MatType max_param=0;	//maximum constitutive parameter value
	MatType min_param=1e10;	//minimum constitutive parameter value
	MatType param_temp;	//constitutive parameter value that has been read
	double param_lower_limit; //lower limit of the constitutive parameter
	if ((constitutive_param_type=="rel_permittivity")||(constitutive_param_type=="rel_permeability"))
	{
		param_lower_limit = 1;
	}
	else if (constitutive_param_type=="electric_conductivity")
	{
		param_lower_limit = 0;
	}
	else
	{
		throw AngoraDeveloperException("unknown constitutive parameter type");
	}
	for (int i=1; i<=xExtent; i++)
	{
		for (int j=1; j<=yExtent; j++)
		{
			for (int k=1; k<=zExtent; k++)
			{
				MaterialFile.read((char*)&param_temp,sizeof(param_temp));	//read the constitutive parameter
				if (param_temp>=param_lower_limit)
				{//if the value is nonpositive, don't bother with it at all
					if (param_temp>max_param) max_param=param_temp;	//update the maximum constitutive parameter
					if (param_temp<min_param) min_param=param_temp;	//update the minimum constitutive parameter
				}
			}
		}
	}
	//maximum number of different material types that can be extracted from the region
//	int max_num_of_materials = 1000; 	//pretty random, may have to find a more efficient way in the future
	//minimum difference in constitutive parameter between different materials
	MatType param_step = (max_param-min_param)/(max_number_of_new_materials-1);
	//add the new materials to the material list
	//before increasing the number of materials, save the current maximum material indices
	int material_index_saved_ex = NumOfElectricMaterials_X-1;
	int material_index_saved_ey = NumOfElectricMaterials_Y-1;
	int material_index_saved_ez = NumOfElectricMaterials_Z-1;
	int material_index_saved_hx = NumOfMagneticMaterials_X-1;
	int material_index_saved_hy = NumOfMagneticMaterials_Y-1;
	int material_index_saved_hz = NumOfMagneticMaterials_Z-1;
	//dummy material index
	MaterialId NewMaterialId;
	if ((constitutive_param_type=="rel_permittivity"))
	{//the relative permittivity of the new material is min_param+(i-1)*param_step, electric_conductivity is 0
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			AddIsotropicElectricMaterial(NewMaterialId, min_param+(i-1)*param_step,0);
		}
	}
	if ((constitutive_param_type=="rel_permeability"))
	{//the relative permeability of the new material is min_param+(i-1)*param_step, magnetic conductivity is 0
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			AddIsotropicMagneticMaterial(NewMaterialId, min_param+(i-1)*param_step,0);
		}
	}
	if ((constitutive_param_type=="electric_conductivity"))
	{//the electric_conductivity (in S units) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			AddIsotropicElectricMaterial(NewMaterialId, 1, min_param+(i-1)*param_step); //relative permittivity is 1
		}
	}

	MaterialFile.seekg(pos_saved,ios::beg);	//return to the saved position in the file

	//now, read through the file again, determine the material index for each point, and update the material indices in the main grid
	int material_offset;	//offset of the current material index beginning from the material index that was saved before the creation of new materials
	int material_index_ex,material_index_ey,material_index_ez;		//absolute material indices (for Ex,Ey,Ez components) in the current material list
	int material_index_hx,material_index_hy,material_index_hz;		//absolute material indices (for Hx,Hy,Hz components) in the current material list

	//calculate the coordinates of the back-left-lower corner of the region
	int xCornerPos=xPos;
	int yCornerPos=yPos;
	int zCornerPos=zPos;
	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
	{
		if (rank==0)
		{
			cout << "Invalid anchor point \"" << anchor << "\" for material input file " << MaterialFileName << " in node " << rank << endl << endl;
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
	Array<MatType,1> material_row(Range(xCornerCell,xCornerCell+xExtent-1));	//temporary array of size xExtent that will hold the constitutive parameter values belonging to a x-row in the 2-D material property array in the file
	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.

	//if the file stores electric properties, update Media_Ex,Media_Ey,Media_Ez
	if ((constitutive_param_type=="rel_permittivity")||(constitutive_param_type=="electric_conductivity"))
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtent-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
			{
				MaterialFile.read((char*)material_row.data(),xExtent*sizeof(material_row(0)));	//read the constitutive parameter x-row into temporary array
				//Media_Ex
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_ex = material_index_saved_ex + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								Media_Ex(i,j,k) = material_index_ex;
							}
						}
					}
				}

				//Media_Ey
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_ey = material_index_saved_ey + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								Media_Ey(i,j,k) = material_index_ey;
							}
						}
					}
				}

				//Media_Ez
				if ((k>=klower)&&(k<=kupper))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_ez = material_index_saved_ez + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								Media_Ez(i,j,k) = material_index_ez;
							}
						}
					}
				}
			}
		}
	}

	//if the file stores magnetic properties, update Media_Hx,Media_Hy,Media_Hz
	if (constitutive_param_type=="rel_permeability")
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtent-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
			{
				//Media_Hx
				if ((k>=klower)&&(k<=kupper))
				{
					if ((j>=jleft)&&(j<=jright))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_hx = material_index_saved_hx + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								Media_Hx(i,j,k) = material_index_hx;
							}
						}
					}
				}

				//Media_Hy
				if ((k>=klower)&&(k<=kupper))
				{
					if ((j>=jleft)&&(j<=jright+1))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_hy = material_index_saved_hy + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								Media_Hy(i,j,k) = material_index_hy;
							}
						}
					}
				}

				//Media_Hz
				if ((k>=klower)&&(k<=kupper+1))
				{
					if ((j>=jleft)&&(j<=jright))
					{
						for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtent-1); i++)
						{
							param_temp = material_row(i);
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_hz = material_index_saved_hz + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								Media_Hz(i,j,k) = material_index_hz;
							}
						}
					}
				}
			}
		}
	}

	//finally, close the file
	MaterialFile.close();
//	if (rank==0)
//	{
//		cout << "Material region read." << endl << endl;
//	}
}

#endif // MATFILE_H

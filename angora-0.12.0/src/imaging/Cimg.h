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

#ifndef CIMG_H
#define CIMG_H

//Declaration of the class "Cimg" for synthesizing an optical image

#include <fstream>
//for the vector STL class
#include <vector>
//for the "find" algorithm from C++ STL
#include<algorithm>

//only the declaration of Ctr_pd needed: use forward declaration
class Ctr_pd;
//only the declaration of Ctfsf needed: use forward declaration
class Ctfsf;

extern int OriginX,OriginY,OriginZ;

const string admissible_list_members_strbuf[] = {"E_x_tot","E_y_tot","E_z_tot","E_x_sca","E_y_sca","E_z_sca","E_x_inc","E_y_inc","E_z_inc","intensity_tot","intensity_sca","intensity_inc"};

class ImgDataType
{
	public:
	ImgDataType(const Array<double,1>& mylambda):
		lambda(mylambda),
		ap_half_angle(M_PI/2),
		magnification(1),
		img_expansionfactor_x(1),
		img_expansionfactor_y(1),
		img_oversamplingrate_x(1),
		img_oversamplingrate_y(1),
		choose_smallest_power_of_2(true),
		img_space_refr_index(1),
		coll_half_space("upper"),
		NFFFTMarginBackX(3),
		NFFFTMarginFrontX(3),
		NFFFTMarginLeftY(3),
		NFFFTMarginRightY(3),
		NFFFTMarginLowerZ(3),
		NFFFTMarginUpperZ(3),
		ImgOriginX(OriginX),
		ImgOriginY(OriginY),
		ImgOriginZ(OriginZ),
//		xz_symmetry(false),
//		yz_symmetry(false),
		admissible_list_members(admissible_list_members_strbuf,
						admissible_list_members_strbuf+sizeof(admissible_list_members_strbuf)/sizeof(string)),
	   	TFSFPtr(NULL)
	{};

	//adds the output data named "output_data_item_name" to the list, throws AngoraInvalidArgumentExceptionWithType<string> if string is not valid
	void add_output_data_item_to_list(const string& output_data_item_name);

	bool output_data_item_exists(const string& output_data_item_name) const
	{//returns true if the data string named "output_data_item_name" is in the list
		return (find(image_array_list.begin(),
					 image_array_list.end(),
					 output_data_item_name)
				!=image_array_list.end());
	}

	const Array<double,1> lambda;	//wavelengths at which the image is calculated (in meters)
	double ap_half_angle;  //half-angle subtended by the collection aperture (in radians)
	double magnification; //magnification of the objective
	double img_expansionfactor_x,img_expansionfactor_y; //ratios of the estimated maximum x and y extents of the field distribution at the object space to the actual x and y extents of the grid [default: 1]. The higher these are, the denser the far-field directions, and the higher the computational cost. However, aliasing is reduced and image quality is improved.
	double img_oversamplingrate_x,img_oversamplingrate_y;	//how much is the final image oversampled spatially? [default: 1]
	bool choose_smallest_power_of_2; //if true, the FFT lengths are the smallest powers of 2 that give sampling periods larger than what are requested
	double img_space_refr_index;	//refractive index of the image space [default:1]
	string coll_half_space; //collection half space (upper or lower)
	int NFFFTMarginBackX,NFFFTMarginFrontX,NFFFTMarginLeftY,NFFFTMarginRightY,NFFFTMarginLowerZ,NFFFTMarginUpperZ;	//distances (in cells) of the NFFFT box from the PML boundary
	double ImgOriginX,ImgOriginY,ImgOriginZ;	 //the index of the cell assigned as the origin of the image (focal point of the objective)
	 //(its rear-left-lower corner is the origin)
	 //(can be fractional, therefore defined as "double")
//	bool xz_symmetry,yz_symmetry;	//true if there is xz (or yz) symmetry in the problem: if so, reflections of the final image is added to itself

	const Ctfsf* TFSFPtr;	//points to a non-changeable CTFSF object

	private:
	bool is_output_data_item_valid(const string& output_data_item_name) const
	{//returns true if the data string named "output_data_item_name" is valid
		return (find(admissible_list_members.begin(),
					admissible_list_members.end(),
					output_data_item_name)
				!=admissible_list_members.end());
	}

	const vector<string> admissible_list_members;
	vector<string> image_array_list; //list of image arrays to be computed and written into file
};

class Cimg
{
 public:
	 Cimg(const ImgDataType& MyData, const string& ImgOutputFileName, const int& Index);	//constructor
	 virtual ~Cimg();	//virtual destructor is needed when deleting derived objects using a base class pointer

	 void UpdateFarField(const int& n);		//Update far-field arrays

	 void ConstructImage();	//Construct far-field from the far-field arrays

 private:
	 const int ImgIndex;	//index of the current image in Cimgs

	 const ImgDataType Data;	//data type that holds the imaging data

	 Ctr_pd* Transformer;	//pointer to the transformer object

	 Array<double,1> lambda;	//wavelength array
	 double theta_max; 	//aperture half-angle (in radians)
	 double SinThetaMax; 	//sine of the aperture half-angle
	 double NA_obj; //numerical aperture on the collection side = obj_space_refr_index*SinThetaMax
	 double NA_img; //numerical aperture on the image side = img_space_refr_index*SinThetaMax
	 double nimg,nobj; //refractive indices of the image and object spaces, respectively
	 Array<double,1> sx,sy;	//direction-cosine arrays

	 double M,Mpr; //spatial magnification and angular demagnification

	 string FarFieldOutputFileName; //full name of the far-field file (including path)
//	 ifstream FarFieldFile;		//far-field file stream (read-only)


	 const string ImageFileName; //name of the image output file
//	 ofstream ImageFile;		//image file stream
};

#endif

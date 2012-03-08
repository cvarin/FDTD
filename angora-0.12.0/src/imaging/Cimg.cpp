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

//Definitions for the class "Cimg" for synthesizing an optical image
//This file includes the constructor method that adds a transformer for the collection step, the on-the-fly FFT updating method and the destructor.

#include "headers.h"

#include "Cimg.h"

//base Angora exception class
#include "angora_excp.h"

//definition of Ctr_pd needed
#include "nffft/pd/Ctr_pd.h"

#include "nffft/pd/Ctr_pd_fs.h"
#include "nffft/pd/Ctr_pd_2l.h"
#include "nffft/pd/Ctr_pd_3l.h"

extern double dx;
extern int NCELLS_X,NCELLS_Y;

extern int number_of_layers;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;

extern int GridIndex;
extern int rank;

extern void MPI_exit(const int& exitcode);


void ImgDataType::add_output_data_item_to_list(const string& output_data_item_name)
{//adds the output data named "output_data_item_name" to the list, returns false if string is not valid
	//const_iterator needed because "admissible_list_members" is a const vector
	vector<string>::const_iterator it = find(admissible_list_members.begin(),
											 admissible_list_members.end(),
											 output_data_item_name);
	if (it==admissible_list_members.end()) //not found in the admissible member list?
	{
//			return false;
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name, output_data_item_name);
	}
	else
	{
		image_array_list.push_back(output_data_item_name);
	}
}

Cimg::Cimg(const ImgDataType& MyData, const string& ImgOutputFileName, const int& Index)
	:ImageFileName(ImgOutputFileName), Data(MyData), ImgIndex(Index)
{
	lambda.resize(Data.lambda.size());
	lambda = Data.lambda;
	theta_max = Data.ap_half_angle;

	M = Data.magnification;
	nimg = Data.img_space_refr_index;

	if ((cos(theta_max)<0)&&(sin(theta_max)<0))
	{
		if (rank==0)
		{
			cout << "Error: Theta half-range (ap_half_angle) in optical image " << ImgIndex << " should be between 0 and 90deg." << endl;
		}
		MPI_exit(-1);
	}

	SinThetaMax = abs(sin(theta_max));

	string directionspec;

	if (Data.coll_half_space=="upper")
	{
		directionspec = "dircosx-dircosy-upper";
		nobj = sqrt(epsilon_r_upper*mu_r_upper);
	}
	else if (Data.coll_half_space=="lower")
	{
		directionspec = "dircosx-dircosy-lower";
		nobj = sqrt(epsilon_r_lower*mu_r_lower);
	}
	else
	{
		if (rank==0)
		{
			cout << "Developer error: Invalid observation half space (" << Data.coll_half_space << ") for optical image " << ImgIndex << " in grid " << GridIndex << endl;
		}
		MPI_exit(-1);
	}

	double lambda_min,lambda_max;
	lambda_min = min(lambda)/nobj;
	lambda_max = max(lambda)/nobj;

	NA_obj = nobj*SinThetaMax;
	NA_img = nimg*SinThetaMax;

	//maximum x and y extents of the field distribution (used to determine the direction-cosine spacing)
	//x and y extents are determined either by the beam waist or the size of the simulation space
	double minimum_aliasing_distance_in_beam_width = 4;
	double beam_waist = minimum_aliasing_distance_in_beam_width*lambda_max/SinThetaMax;  //FIXME:should find a better way for this later

	if (minimum_aliasing_distance_in_beam_width<3)
		if (rank==0) cout << "Warning: There may be aliasing!!!!" << endl;

	double W_x = max(beam_waist,Data.img_expansionfactor_x*(NCELLS_X*dx));
	double W_y = max(beam_waist,Data.img_expansionfactor_y*(NCELLS_Y*dx));

	//use sampling theorem for default spacing
	double dsx,dsy;
	int N_sx,N_sy;
	dsx = lambda_min/W_x;
	dsy = lambda_min/W_y;
	N_sx = int(2*SinThetaMax/dsx)+1;	//sx is at the midpoint of each patch
	N_sy = int(2*SinThetaMax/dsy)+1;	//sy is at the midpoint of each patch
	Array<double,1> sx,sy;
	sx.resize(N_sx);
	sy.resize(N_sy);

// cout << W_x << "," << W_y << endl;
// cout << lambda_min << "," << lambda_max << endl;
// cout << dsx << "," << dsy << endl;
// cout << N_sx << "," << N_sy << endl;
// cout << SinThetaMax << endl;
//
// sx = -SinThetaMax*(1-1.0/N_X1) + sx_index*dsx;
// sy = -SinThetaMax*(1-1.0/N_X2) + sy_index*dsy;

	//construct direction-cosine arrays
	for (int sxindex=0; sxindex<N_sx; sxindex++)
	{
		sx(sxindex) = dsx*(sxindex-(N_sx-1.0)/2.0);
	}
	for (int syindex=0; syindex<N_sy; syindex++)
	{
		sy(syindex) = dsy*(syindex-(N_sy-1.0)/2.0);
	}

// cout << sx << "," << sy << endl;
// cout << SinThetaMax << endl;

	//transformer data object
	TrDataType_pd TransformerData(lambda,sx,sy,SinThetaMax,directionspec);

	/** VERY IMPORTANT: Scale the direction cosines with wavelength for broadband imaging **/

	TransformerData.scale_with_wavelength = true;

	/** VERY IMPORTANT: Scale the direction cosines with wavelength for broadband imaging **/

	TransformerData.NFFFTMarginBackX = Data.NFFFTMarginBackX;
	TransformerData.NFFFTMarginFrontX = Data.NFFFTMarginFrontX;
	TransformerData.NFFFTMarginLeftY = Data.NFFFTMarginLeftY;
	TransformerData.NFFFTMarginRightY = Data.NFFFTMarginRightY;
	TransformerData.NFFFTMarginLowerZ = Data.NFFFTMarginLowerZ;
	TransformerData.NFFFTMarginUpperZ = Data.NFFFTMarginUpperZ;

	TransformerData.NFFFTOriginX = Data.ImgOriginX;
	TransformerData.NFFFTOriginY = Data.ImgOriginY;
	TransformerData.NFFFTOriginZ = Data.ImgOriginZ;

	TransformerData.PointSourcesPtr = NULL;
	TransformerData.TFSFPtr = Data.TFSFPtr;


	ostringstream FarFieldOutputFileNameStream;
	FarFieldOutputFileNameStream << ImageFileName << ".ff";
	FarFieldOutputFileName = FarFieldOutputFileNameStream.str();

	try{
	if (number_of_layers==1)
		{//attach a free-space transformer
			Transformer = new Ctr_pd_fs(TransformerData,FarFieldOutputFileName,0);
		}
		else if (number_of_layers==2)
		{//attach a 2-layered medium transformer
			Transformer = new Ctr_pd_2l(TransformerData,FarFieldOutputFileName,0);
		}
		else if (number_of_layers==3)
		{//attach a 3-layered medium transformer
			Transformer = new Ctr_pd_3l(TransformerData,FarFieldOutputFileName,0);
		}
	}
	catch(bad_alloc& e)
	{
		/** TODO: Handle memory exception later **/
	};
}

void Cimg::UpdateFarField(const int& n)
{
	//update the far-field arrays in (*Transformer)
	Transformer->UpdateFarField(n);
}

Cimg::~Cimg()
{//release memory allocated to Transformer
	delete Transformer;
}

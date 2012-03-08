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

#ifndef CKOHLER_H
#define CKOHLER_H

//Declaration of the class "Ckohler" for a TF/SF Kohler-beam source

//only the declaration of Cpw needed: use forward declaration
class Cpw;

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

#define UNINITIALIZED_RANDOM_SEED -1

extern int OriginX,OriginY,OriginZ;


struct KBDataType
{
	KBDataType():
		simulation_type("deterministic"),
		angular_discretization("cartesian"),
// 		theta_max(23.5782*M_PI/180),
		n_1(-1), //-1 used as a marker to indicate that the value has not been initialized
		n_2(-1), //-1 used as a marker to indicate that the value has not been initialized
// 		f(0.01),
		pw_pol(0),
		KBMarginBackX(6),
		KBMarginFrontX(6),
		KBMarginLeftY(6),
		KBMarginRightY(6),
		KBMarginLowerZ(6),
		KBMarginUpperZ(6),
		KBOriginX(OriginX),
		KBOriginY(OriginY),
		KBOriginZ(OriginZ),
		L_req(15),

		DisplayWarnings(true),
		DisplayCellsPerLambda(false),

		E0(1)
	{};

	string simulation_type;	//method for dividing the stochastic beam into deterministic FDTD simulations
	string angular_discretization;	//method for discretizing the plane wave directions inside the NA
	double theta_max;	//maximum theta value of the incidence cone (in radians)
	int n_1; //for "cartesian": number of x-direction cosines;  For "radial": number of theta between -theta_max (not included) and theta_max (not included)
	int n_2; //for "cartesian": number of y-direction cosines;  For "radial": number of phi between 0 (not included) and pi (not included)
//	double f;	//focal length (in m)
	double pw_pol;	//polarization angle for the incident plane wave (in radians) (normally incident in the -z direction, psi measured ccw from +x axis -- like the phi azimuthal angle in spherical coordinates)

	//distances (in cells) of the TF/SF box from the PML boundary
	int KBMarginBackX,KBMarginFrontX,KBMarginLeftY,KBMarginRightY,KBMarginLowerZ,KBMarginUpperZ;
	//the coordinates of the origin of the TF/SF box //(measured from the rear-left-lower corner of the grid)
	double KBOriginX,KBOriginY,KBOriginZ;
	double L_req;	//minimum number of grid cells required per minimum wavelength (checked against if warnings are enabled)

	bool DisplayWarnings;	//are warnings displayed when (grid cells/wavelength) is low?
	bool DisplayCellsPerLambda;		//is the number of cells per lambda displayed?

// 	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the beam origin
// 	//	defined by KBOriginY,KBOriginZ
	const_Cwf_shared_ptr waveform;	//pointer to the waveform object Cwf, representing the beam waveform. The Cwf object cannot be changed using this pointer.
	double E0;	//extra amplitude factor applied to the waveform (1 by default)
};

class Ckohler{
 public:
	 Ckohler(const KBDataType& MyData, const int& Index);	//constructor

	 //Applies field corrections on the TF/SF box
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

	 int NumberOfPlaneWaves() const
	 {
		 return PlaneWaves.size();
	 };

 	 void WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Ckohler object
	 void WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Ckohler object
	 void WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Ckohler object

 private:
	 const int KBIndex;	//index of the current beam in the TF/SF object Ctfsf

	 const KBDataType Data;	//Kohler-beam parameters

//	 double d_theta,d_phi;	//uniform spacings between discrete positions (if applicable)
//	 double theta,phi,psi;	//incidence angles and the polarization angle for a single plane-wave component

	 Array<double,2> pwfactor;	//factor multiplied by each plane wave (dsx1*dsx2 for Cartesian, |sin(theta)|*w_theta*dphi for radial
	 Array<bool,2> angle_within_ill_cone;	//incidence angles and the polarization angle for a single plane-wave component
	 Array<double,2> theta_array,phi_array;	//incidence angles and the polarization angle for a single plane-wave component
	 double theta,phi,psi;	//incidence angles and polarization angle for a single plane-wave component
	 int N_X1, N_X2;//these correspond to n_1 and n_2 in FBDataType
	 //cartesian placement:
// 	 int Nsx,Nsy;	//these correspond to n_1 and n_2 in FBDataType
	 double dsx,dsy;	//uniform spacings between direction cosines
	 double sx,sy;	//direction cosines for Cartesian placement
	 //radial placement:
// 	 int Ntheta,Nphi;	//these correspond to n_1 and n_2 in FBDataType
	 double d_phi;	//uniform spacing between phi values
	 Array<double,1> GLx,GLw;	//Gauss-Legendre quadrature points and weights between [-1,1]

	 int total_num_of_pws; //total number of plane waves
	 int this_pw; //index of current plane wave, depending on the grid being simulated

	 //*********** Plane wave components ***************//
	 vector<boost::shared_ptr<Cpw> > PlaneWaves;	//array of pointers to plane wave objects
	 //*********** Plane wave components ***************//
//	 boost::shared_ptr<Cpw> PlaneWave;
};

#endif

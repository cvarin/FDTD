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

#ifndef CHGP_H
#define CHGP_H

//Declaration of the class "Chgb" for a TF/SF Hermite-Gaussian-beam source

//only the declaration of Cpw needed: use forward declaration
class Cpw;

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

extern int OriginX,OriginY,OriginZ;


struct HGBDataType
{
	HGBDataType():
//		n_sx(15),
//		n_sy(15),
		x_order(0),
		y_order(0),
		direction("+z"),
		beam_halfwidth(1),
		polarization(0),
		HGBMarginBackX(6),
		HGBMarginFrontX(6),
		HGBMarginLeftY(6),
		HGBMarginRightY(6),
		HGBMarginLowerZ(6),
		HGBMarginUpperZ(6),
		HGBOriginX(OriginX),
		HGBOriginY(OriginY),
		HGBOriginZ(OriginZ),
		L_req(15),

		DisplayWarnings(true),
		DisplayCellsPerLambda(false),

		E0(1)
	{};

//	int n_sx;	//number of x-direction cosines
//	int n_sy;	//number of y-direction cosines
	int x_order; //order of the Hermite polynomial in the x direction
	int y_order; //order of the Hermite polynomial in the y direction
	string direction; //direction of the beam ("+z" or "-z")
	double beam_halfwidth; //spatial half-width of the FIELD-AMPLITUDE profile at the center frequency of the excitation waveform (similar to tau in waveforms) (in m)  Since the Hermite-Gaussian beam is mostly a monochromatic beam, the behavior at off-center frequencies is of secondary interest. It is assumed here that the beam width is PROPORTIONAL to the wavelength. This allows the representation of the beam by a fixed set of plane waves in the angle space.
	double polarization; //polarization angle of the E field w.r.t. the x-axis (positive in +y direction) (in radians)

	//distances (in cells) of the TF/SF box from the PML boundary
	int HGBMarginBackX,HGBMarginFrontX,HGBMarginLeftY,HGBMarginRightY,HGBMarginLowerZ,HGBMarginUpperZ;
	//the coordinates of the origin of the TF/SF box //(measured from the rear-left-lower corner of the grid)
	double HGBOriginX,HGBOriginY,HGBOriginZ;
	double L_req;	//minimum number of grid cells required per minimum wavelength (checked against if warnings are enabled)

	bool DisplayWarnings;	//are warnings displayed when (grid cells/wavelength) is low?
	bool DisplayCellsPerLambda;		//is the number of cells per lambda displayed?

// 	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the Hermite-Gaussian-beam origin
// 	//	defined by HGBOriginY,HGBOriginZ
	const_Cwf_shared_ptr waveform;	//pointer to the waveform object Cwf, representing the Hermite-Gaussian-beam waveform. The Cwf object cannot be changed using this pointer.
	double E0;	//extra amplitude factor applied to the waveform (1 by default)

// 	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the plane-wave origin
// 	//	defined by PWOriginX,PWOriginY,PWOriginZ
// 	double E0;	//amplitude of the plane wave
// 	string WaveShape;	//shape of the plane-wave waveform
// 	double tau;	//time constant for the Gaussian waveform
// 	double w_0;	//angular modulation frequency (if modulated)
// 	double delay;	//delay of the waveform (in tau) (should be increased from within main(), since OriginDelay is always positive)
};

class Chgb
{
 public:
	 Chgb(const HGBDataType& MyData, const int& Index);	//constructor

	 //Applies field corrections on the TF/SF box
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

	 int NumberOfPlaneWaves() const
	 {
		 return PlaneWaves.size();
	 };

 	 void WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Chgb object
	 void WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Chgb object
	 void WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Chgb object

 private:
	 const int HGPIndex;	//index of the current Hermite-Gaussian beam in the TF/SF object Ctfsf	(1st, 2nd,.. etc.)

	 const HGBDataType Data;	//Hermite-Gaussian-beam parameters

//	 int n_sx;	//number of x-direction cosines
//	 int n_sy;	//number of y-direction cosines
	 int N_X1, N_X2; //number of x and y direction cosines

	 double dsx,dsy;	//uniform direction-cosine spacings
	 double theta,phi,psi;	//incidence angles and the polarization angle for a single plane-wave component

	 //*********** Plane wave components ***************//
	 vector<boost::shared_ptr<Cpw> > PlaneWaves;	//array of pointers to plane wave objects
	 //*********** Plane wave components ***************//
};

#endif

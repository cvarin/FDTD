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

#ifndef CGSMP_H
#define CGSMP_H

//Declaration of the class "Cgsmb" for a TF/SF isotropic Gaussian Schell-model beam

//only the declaration of Chgb needed: use forward declaration
class Chgb;

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the vector STL class
#include <vector>

//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

extern int OriginX,OriginY,OriginZ;


struct GSMBDataType
{
	GSMBDataType():
		n_mode_x(0),
		n_mode_y(0),
		direction("+z"),
		beam_halfwidth(1),
		correlation_length(1),
		polarization(0),
		GSMBMarginBackX(6),
		GSMBMarginFrontX(6),
		GSMBMarginLeftY(6),
		GSMBMarginRightY(6),
		GSMBMarginLowerZ(6),
		GSMBMarginUpperZ(6),
		GSMBOriginX(OriginX),
		GSMBOriginY(OriginY),
		GSMBOriginZ(OriginZ),
		L_req(15),

		DisplayWarnings(true),
		DisplayCellsPerLambda(false),

		E0(1)
	{};

	int n_mode_x;	//number of modes in the x direction
	int n_mode_y;	//number of modes in the y direction
	string direction; //direction of the beam ("+z" or "-z")
	double beam_halfwidth; //half-width of the Gaussian part of the INTENSITY profile (similar to tau in waveforms) (in m)
	double correlation_length; //correlation length of the Gaussian correlation function (similar to lc in exponential) (in m)
	double polarization; //polarization angle of the E field w.r.t. the x-axis (positive in +y direction) (in radians)

	//distances (in cells) of the TF/SF box from the PML boundary
	int GSMBMarginBackX,GSMBMarginFrontX,GSMBMarginLeftY,GSMBMarginRightY,GSMBMarginLowerZ,GSMBMarginUpperZ;
	//the coordinates of the origin of the TF/SF box //(measured from the rear-left-lower corner of the grid)
	double GSMBOriginX,GSMBOriginY,GSMBOriginZ;
	double L_req;	//minimum number of grid cells required per minimum wavelength (checked against if warnings are enabled)

	bool DisplayWarnings;	//are warnings displayed when (grid cells/wavelength) is low?
	bool DisplayCellsPerLambda;		//is the number of cells per lambda displayed?

// 	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the Gaussian-Schell-model-beam origin
// 	//	defined by GSMBOriginY,GSMBOriginZ
	const_Cwf_shared_ptr waveform;	//pointer to the waveform object Cwf, representing the Gaussian-Schell-model-beam waveform. The Cwf object cannot be changed using this pointer.
	double E0;	//extra amplitude factor applied to the waveform (1 by default)

// 	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the plane-wave origin
// 	//	defined by PWOriginX,PWOriginY,PWOriginZ
// 	double E0;	//amplitude of the plane wave
// 	string WaveShape;	//shape of the plane-wave waveform
// 	double tau;	//time constant for the Gaussian waveform
// 	double w_0;	//angular modulation frequency (if modulated)
// 	double delay;	//delay of the waveform (in tau) (should be increased from within main(), since OriginDelay is always positive)
};

class Cgsmb
{
 public:
	 Cgsmb(const GSMBDataType& MyData, const int& Index);	//constructor

	 //Applies field corrections on the TF/SF box
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

	 int NumberOfHermiteGaussianBeams() const
	 {
		 return HermiteGaussianBeams.size();
	 };

 	 void WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Cgsmb object
	 void WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Cgsmb object
	 void WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Cgsmb object

 private:
	 const int GSMPIndex;	//index of the current Gaussian-Schell-model beam in the TF/SF object Ctfsf	(1st, 2nd,.. etc.)

	 const GSMBDataType Data;	//Gaussian-Schell-model-beam parameters

	 int total_num_of_modes; //total number of modes
	 int this_mode; //index of current mode, depending on the grid being simulated

	 double aa,bb; //the a and b parameters, where a=1/(4*beam_halfwidth^2), b=1/(2*correlation_length^2)
	 double cc; //the c-parameter = (a^2+2ab)^1/2

	 //*********** Plane wave components ***************//
	 vector<boost::shared_ptr<Chgb> > HermiteGaussianBeams;	//array of pointers to Hermite-Gaussian beams
	 //*********** Plane wave components ***************//
};

#endif

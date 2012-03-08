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

#ifndef CPW_H
#define CPW_H

//Declaration of the abstract base class "Cpw" for a TF/SF plane-wave source

//base Angora exception class
#include "angora_excp.h"

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

#include <fstream>

class InvalidTFSFBounds: public AngoraException
{// the exception class for invalid TF/SF box bounds
public:
  virtual const string getError() const
  {
    return "bounds of the TF/SF surface are invalid";
  }
};


extern int OriginX,OriginY,OriginZ;

struct PWDataType
{
	PWDataType():
// 		THETA(0*M_PI/180),
// 		PHI(0*M_PI/180),
// 		PSI(45*M_PI/180),
		PWMarginBackX(6),
		PWMarginFrontX(6),
		PWMarginLeftY(6),
		PWMarginRightY(6),
		PWMarginLowerZ(6),
		PWMarginUpperZ(6),
		PWOriginX(OriginX),
		PWOriginY(OriginY),
		PWOriginZ(OriginZ),
		L_req(15),

		DisplayWarnings(true),
		DisplayCellsPerLambda(false),

		E0(1),
		extra_phase(0)
// 		differentiate(false),
// 		integrate(false)
	{};

	double THETA,PHI,PSI;	//incidence angles
	//distances (in cells) of the TF/SF box from the PML boundary
	int PWMarginBackX,PWMarginFrontX,PWMarginLeftY,PWMarginRightY,PWMarginLowerZ,PWMarginUpperZ;
	//the coordinates of the origin of the TF/SF box //(measured from the rear-left-lower corner of the grid)
	double PWOriginX,PWOriginY,PWOriginZ;

	double L_req;	//minimum number of grid cells required per minimum wavelength (checked against if warnings are enabled)

	bool DisplayWarnings;	//are warnings displayed when (grid cells/wavelength) is low?
	bool DisplayCellsPerLambda;		//is the number of cells per lambda displayed?

	//the following parameters for the incident waveform are always referenced to (i.e. measured at) the plane-wave origin
	//	defined by PWOriginX,PWOriginY,PWOriginZ
	const_Cwf_shared_ptr waveform;	//boost shared pointer to the waveform object Cwf, representing the plane-wave waveform. The Cwf object cannot be changed using this pointer.
	double E0;	//extra amplitude factor applied to the waveform (1 by default)
	double extra_phase; //extra phase (in rad) applied to the waveform (using the analytic waveform obtained via Hilbert transform -- this is only viable for waveforms having no DC (w=0) power)
// 	bool differentiate;		//take the time derivative of the waveform
// 	bool integrate;		//take the time integral of the waveform
};

class Cpw
{
 public:
	 Cpw(const PWDataType& MyData);	//constructor
 	 virtual ~Cpw(){};	//virtual destructor for cleaning up some objects (e.g. ofstream object)
						//this is needed when deleting derived objects using a base class pointer

	 //Applies field corrections on the TF/SF box
	 virtual void CorrectE(const int& n);	//these are overloaded in Cpw_2l
	 virtual void CorrectH(const int& n);	//these are overloaded in Cpw_2l

	 //method that returns the number of cells/min.lambda
	 double CellsPerLambda()
	 {
	 	return L;
	 }

	 virtual void WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const=0;	//cannot modify the Cpw object
	 virtual void WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const=0; //cannot modify the Cpw object
	 virtual void WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const=0;	//cannot modify the Cpw object

 protected:

	 //These virtual methods need to be implemented by the specific plane-wave source (free-space, layered, etc)
	 //These methods update the auxiliary grids used for generating the incident field
	 virtual void UpdateIncidentE(const int& n)=0;
	 virtual void UpdateIncidentH(const int& n)=0;
	 //These methods apply field corrections in the main grid (by projection onto the auxiliary grid or any other method)
	 //Previously, they were defined as base methods, and the incident field at each point on the TFSF box was calculated via a virtual function call.
	 //The function call couldn't be inlined, since it was virtual. Now, these functions are copied in each derived class at the expense of code duplication.
	 //This enables the inlining of "IncidentEx" etc. in each derived class, and up to 25% increase in speed.
	 virtual void ApplyCorrectionE(const int& n)=0;
	 virtual void ApplyCorrectionH(const int& n)=0;

	 const PWDataType Data;	//plane-wave parameters

	 //incident field information
	 double THETA_INC, PHI_INC;		//direction of the incident field (z-component always downward-directed)
	 double SinT,CosT,SinP,CosP,SinPsi,CosPsi;	//necessary sines and cosines

	 double PSI;	//(for LP incident field, if applicable) angle of linearly polarized incident electric field w.r.t.
					// (k_inc)x(z) in clockwise direction (k_inc points away from the clock surface)

	 //Indices of the cells right inside the plane-wave boundary
	 int PWBackX,PWFrontX,PWLeftY,PWRightY,PWLowerZ,PWUpperZ;
	 //Minimum and maximum indices cells contained in the TFSF box and in the current node
	 int TFSF_min_x,TFSF_max_x,TFSF_min_y,TFSF_max_y,TFSF_min_z,TFSF_max_z;

	 TinyVector<double,3> k_inc;	//unit vector in the direction of incidence
	 TinyVector<double,3> k_inc_lateral;	//lateral unit vector in the direction of incidence
	 TinyVector<double,3> unit_z;		//unit vector in the z-direction
	 double k_inc_x,k_inc_y,k_inc_z;	//x,y,z components of k_inc

	 TinyVector<double,3> PWOriginationPoint;	//coordinates of the point that the plane wave originates from
						//(measured from the back-left-bottom point of the grid)
						//(not the same as PWContactPoint, since the z component of PWOriginationPoint is fixed)
	 double PWOriginationPoint_x,PWOriginationPoint_y,PWOriginationPoint_z;	//x,y,z components of k_inc

	 TinyVector<double,3> PWOrigin;	//coordinates of the origin of the TF/SF box
						//(measured from the back-left-bottom point of the grid)
	 double OriginDelay;	//the time it takes for the plane wave to reach the origin

	 double Hinc, Einc;	//temporary incident field values

	 double L;		//cells per minimum lambda in the incident waveform

	 //grid velocity at (THETA_INC,PHI_INC)
	 double GridVelocity(const double& theta, const double& phi);

	 ///this may be removed later
	 //the angular modulation frequency (if modulated)
	 double w_0;
	 //the modulation frequency
	 double f_0;
	 ///this may be removed later
};

#endif

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

#ifndef CPW_2L_H
#define CPW_2L_H

//Declaration of the class "Cpw_2l" for a TF/SF plane-wave source in a lossless 2-layered medium

//base Angora exception class
#include "angora_excp.h"

#include "Cpw_fs.h"

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class EvanescentPWException: public AngoraException
{// the exception class for signaling the presence of an evanescent plane wave
public:
  virtual const string getError() const
  {
    return "Warning: Evanescent plane waves cannot be created with Cpw_2l, using Cpw_ml instead.";
  }
};

//Derived class of the free-space plane-wave class, defined only in a half space
class Cpw_2l_hs : public Cpw_fs
{
 public:
	 Cpw_2l_hs(const PWDataType& MyData, const string& myhalfspace, const int& myHighestIndexInLowerHalfSpace, const double& epsilon_r_space);

 private:
	 //Given the coordinates of the field component on the TF/SF boundary, the functions below calculate the incident field by projection onto the 1-D grid (only if the coordinates fall within the desired half space)
	 inline double IncidentHx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHz(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEz(const int& i, const int& j, const int& k, const int& n);

	 //these are copied from Cpw_fs, to be able to define IncidentEx,etc. inline (by having them in the same file as ApplyCorrectionE,H)
	 void ApplyCorrectionE(const int& n);
	 void ApplyCorrectionH(const int& n);

	 const string halfspace;	//the half space over which the plane wave is defined ("upper" or "lower")
	 const int HighestIndexInLowerHalfSpace; 	//highest cell index in the lower half space in the 2-layer medium

	 //is the given field component in the desired half space?
	 Array<bool,1> IsFullIntegerComponentInHalfSpace;		//for components at full-integer z-positions
	 Array<bool,1> IsHalfIntegerComponentInHalfSpace;		//for components at half-integer z-positions
};

class Cpw_2l : public Cpw
{
 public:
	 Cpw_2l(const PWDataType& MyData);

	 //Applies field corrections on the TF/SF box (overloaded from the base class Cpw)
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

 	 void WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Cpw_2l object
	 void WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Cpw_2l object
	 void WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Cpw_2l object

 private:
	 //**********************************************************************************************************************//
	 //These are defined (since they are pure virtual in Cpw), but left blank                                                //
	 // What is done by these functions is instead done by CorrectE, CorrectH (virtual Cpw functions, redefined in Cpw_2l)   //
	 //**********************************************************************************************************************//
	 void UpdateIncidentE(const int& n);
	 void UpdateIncidentH(const int& n);
	 void ApplyCorrectionE(const int& n);
	 void ApplyCorrectionH(const int& n);
	//**********************************************************************************************************************//
	//**********************************************************************************************************************//
	//**********************************************************************************************************************//

	 string incidence_layer;	//the layer from which the plane wave is incident ("upper" or "lower")

	 double eps_r_i;	//relative permittivity of the layer from which the plane wave is incident
	 double eps_r_t;	//relative permittivity of the layer into which the plane wave is transmitted

	 double CosT_i,CosT_t;	//cosines of the propagation angles at the incidence and transmission layers

	 double Refl_TE,Trans_TE;	//TE refl./trans. coefficients
	 double Refl_TM,Trans_TM;	//TM refl./trans. coefficients

	 double E_TE,E_TM;	//TE and TM components of the incident E-field vector
	 double E_TE_refl,E_TM_refl;	//TE and TM components of the reflected E-field vector
	 double E_TE_trans,E_TM_trans;	//TE and TM components of the transmitted E-field vector
	 double E_refl,E_trans;	//amplitudes of the reflected and transmitted E-field vectors

	 double THETA_refl,PHI_refl,THETA_trans,PHI_trans;	//direction angles of the reflected and transmitted waves
	 double PSI_refl,PSI_trans;	//polarization angles of the reflected and transmitted E-field vectors

	 //highest cell index in the lower half space
	 int HighestIndexInLowerHalfSpace;

	 //*********** Plane wave components ***************//
	 vector<boost::shared_ptr<Cpw_2l_hs> > PlaneWaves;	//array of pointers to plane wave objects (incident, reflected, transmitted)
	 //*********** Plane wave components ***************//
};

#endif

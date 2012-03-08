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

#ifndef CPOINTSOURCES_H
#define CPOINTSOURCES_H

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class Celectricdipole
{
 public:
	 Celectricdipole(const int& xPos,  const int& yPos, const int& zPos,
		 const string& myOrientation, const double& myj0, const_Cwf_shared_ptr mywaveform, const int& Index);
	 void ApplyElectricDipole(const int& n);

	 bool ElectricDipoleBelongsToNode;		//is the electric dipole within the node?

	 const int x,y,z;	//position of the electric dipole
	 const double j0; // extra amplitude applied to the current-moment waveform (in A*m)
	 const string Orientation;	//orientation

	 //value of the current moment at time t
	 double current_moment_value(const double& t);

	 //Fourier transform of the current moment at (radian) frequency w
	 complex<double> current_moment_Fourier_transform(const double& w);

	 //Fourier component of the current moment at (radian) frequency w [= Fourier trans./(2pi)]
	 complex<double> current_moment_Fourier_component(const double& w);

	 //smart pointer to a new Cwf object representing the current-moment waveform
	 Cwf_shared_ptr current_moment_waveform();

 private:
	 int ElectricDipoleIndex;		//index of the electric dipole

	 const_Cwf_shared_ptr waveform;	//pointer to the waveform object Cwf, representing the source waveform. The Cwf object cannot be changed using this pointer.

	 double epsilon;	//permittivity at the source position
	 TinyVector<int,3> SourceCell;
//	 Array<double,1> waveformarray;	//source waveform
	 double moment_to_current_density; //factor to convert from current moment to 3D current density
	 double waveform_value(const double& n); //value of the waveform at the n'th time step
};

class Cpointsources
{
 public:
	 int AddElectricDipole(const int& xPos,  const int& yPos, const int& zPos,
		 const string& myOrientation, const double& myj0, const_Cwf_shared_ptr mywaveform);

	 void ApplySources(const int& n);

	 boost::shared_ptr<Celectricdipole> ElectricDipole(const int& ElectricDipoleIndex) const
	 {
		 return ElectricDipoles[ElectricDipoleIndex];
	 }

	 int NumberOfElectricDipoles() const
	 {
		 return ElectricDipoles.size();
	 }

 private:
	 vector<boost::shared_ptr<Celectricdipole> > ElectricDipoles;	//array of pointers to the electric dipoles
};


#endif

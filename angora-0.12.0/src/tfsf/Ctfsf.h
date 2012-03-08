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

#ifndef CTFSF_H
#define CTFSF_H

//Declaration of the class "Ctfsf" for a collection of total-field/scattered-field (TF/SF) incident wave sources.

//only the declaration of Cpw needed: use forward declaration
class Cpw;
//only the declaration of Cfb needed: use forward declaration
class Cfb;
//only the declaration of Chgb needed: use forward declaration
class Chgb;
//only the declaration of Cgsmb needed: use forward declaration
class Cgsmb;
//only the declaration of Ckohler needed: use forward declaration
class Ckohler;
//only the declaration of PWDataType needed: use forward declaration
class PWDataType;
//only the declaration of FBDataType needed: use forward declaration
class FBDataType;
//only the declaration of HGBDataType needed: use forward declaration
class HGBDataType;
//only the declaration of GSMBDataType needed: use forward declaration
class GSMBDataType;
//only the declaration of KBDataType needed: use forward declaration
class KBDataType;

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class Ctfsf
{
 public:
	 //Applies field corrections on the TF/SF box
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

	 //*********** Plane wave sources ***************//
	 int AddPlaneWave(const PWDataType& MyData);	//adds a plane wave
	 int NumberOfPlaneWaves() const
	 {
		 return PlaneWaves.size();
	 };
	 //*********** Plane wave sources ***************//

	 //*********** Focused-beam source ***************//
	 int AddFocusedBeam(const FBDataType& MyData);	//adds a focused beam
	 int NumberOfFocusedBeams() const
	 {
		 return FocusedBeams.size();
	 };
	 //*********** Focused-beam source ***************//

 	 //*********** Hermite-Gaussian-beam source ***************//
	 int AddHermiteGaussianBeam(const HGBDataType& MyData);	//adds a Hermite-Gaussian beam
	 int NumberOfHermiteGaussianBeams() const
	 {
		 return HermiteGaussianBeams.size();
	 };
	 //*********** Hermite-Gaussian-beam source ***************//

 	 //*********** Gaussian-Schell-model source ***************//
	 int AddGaussianSchellModelBeam(const GSMBDataType& MyData);	//adds a Gaussian-Schell-model beam
	 int NumberOfGaussianSchellModelBeams() const
	 {
		 return GaussianSchellModelBeams.size();
	 };
	 //*********** Gaussian-Schell-model source ***************//

	 //*********** Kohler-beam source ***************//
	 int AddKohlerBeam(const KBDataType& MyData);	//adds Kohler beam
	 int NumberOfKohlerBeams() const
	 {
		 return KohlerBeams.size();
	 };	 //*********** Kohler-beam source ***************//

	 void WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Ctfsf object
	 void WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Ctfsf object
	 void WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Ctfsf object

 private:

	 //*********** Plane wave sources ***************//
	 vector<boost::shared_ptr<Cpw> > PlaneWaves;	//array of pointers to plane wave objects
	 //*********** Plane wave sources ***************//

	 //*********** Focused-beam source ***************//
	 vector<boost::shared_ptr<Cfb> > FocusedBeams;	//array of pointers to focused-beam objects
	 //*********** Focused-beam source ***************//

	 //*********** Hermite-Gaussian-beam source ***************//
	 vector<boost::shared_ptr<Chgb> > HermiteGaussianBeams;	//array of pointers to Hermite-Gaussian-beam objects
	 //*********** Hermite-Gaussian-beam source ***************//

	 //*********** Hermite-Gaussian-beam source ***************//
	 vector<boost::shared_ptr<Cgsmb> > GaussianSchellModelBeams;	//array of pointers to Gaussian-Schell-model-beam objects
	 //*********** Hermite-Gaussian-beam source ***************//

	 //*********** Kohler-beam sources ***************//
	 vector<boost::shared_ptr<Ckohler> > KohlerBeams;	//array of pointers to the Kohler-beam objects
	 //*********** Kohler-beam sources ***************//
};

#endif

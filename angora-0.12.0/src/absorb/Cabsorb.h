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

#ifndef CABSORB_H
#define CABSORB_H

//declaration of the class that implements absorbing boundary conditions

//only the declaration of Cpw needed: use forward declaration
class Cpml;

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class Cabsorb
{
 public:
	 void UpdateE();		//Do the necessary E-field updates at the boundaries
	 void UpdateH();		//Do the necessary H-field updates at the boundaries

	 //*********** Convolution PMLs ***************//
	 int AddCPMLs(const int& pml_thickness, const double& SizeOfScatterer);	//adds CPML layers all around
	 int NumberOfCPMLs() const
	 {
		 return CPMLs.size();
	 };
	 //*********** Convolution PMLs ***************//

 private:

	 //*********** Convolution PMLs ***************//
	 vector<boost::shared_ptr<Cpml> > CPMLs;	//array of pointers to convolution PML objects
	 //*********** Convolution PMLs ***************//

};

#endif // CABSORB_H

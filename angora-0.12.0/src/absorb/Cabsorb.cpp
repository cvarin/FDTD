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

//definition of the class that implements absorbing boundary conditions

#include "headers.h"

#include "Cabsorb.h"

#include "CPML.h"


void Cabsorb::UpdateE()
{
	//*********** Convolution PMLs ***************//
	for (int i=0; i<NumberOfCPMLs(); i++)
	{
		CPMLs[i]->UpdateE();		//apply CPML updates
	}
	//*********** Convolution PMLs ***************//
}

void Cabsorb::UpdateH()
{
	//*********** Convolution PMLs ***************//
	for (int i=0; i<NumberOfCPMLs(); i++)
	{
		CPMLs[i]->UpdateH();		//apply CPML updates
	}
	//*********** Convolution PMLs ***************//
}

int Cabsorb::AddCPMLs(const int& pml_thickness, const double& SizeOfScatterer)
{//adds CPML layers on all faces of the grid boundary, returns the CPML index in the collection
	boost::shared_ptr<Cpml> new_cmpl_ptr(new Cpml(pml_thickness, SizeOfScatterer));
	CPMLs.push_back(new_cmpl_ptr);
	return NumberOfCPMLs()-1;
}

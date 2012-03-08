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

#ifndef CREC_H
#define CREC_H

//Declaration and definition of the abstract recorder base class "Crec"

//base Angora exception class
#include "angora_excp.h"

class Crec
{//recorder abstract base class
 public:
 	 virtual ~Crec(){};	//virtual destructor for cleaning up some objects (e.g. ofstream object)
						//this is needed when deleting derived objects using a base class pointer
	 virtual void Record() =0;
};

#endif

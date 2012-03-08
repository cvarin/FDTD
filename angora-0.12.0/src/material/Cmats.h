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

#ifndef CMATS_H
#define CMATS_H

//Declaration of the class "Cmats" that holds the pointers to all the Cmat (and derived) objects
//A good reason to define this as a class is the automatic deletion of dynamically-allocated "Cmat" objects upon destruction.

#include "headers.h"

#include "Cmat_excp.h"

//for the C++ STL map class
#include <map>

//for definition of MaterialId
#include "material_id.h"


class Cmats
{
 public:
	 void CreateMaterial(const MaterialId& mat_id, const string& material_tag);

	 bool lookupMaterialWithTag(const string& material_tag, MaterialId& mat_id) const;  //copies the material identifier that corresponds to tag into mat_id, returns false if tag is not found

	 const MaterialId operator[] (const string& material_tag) const; //returns the material ID corresponding to the material tag, throws exception if not found

	 int NumberOfMaterials() const
	 {
		 return NamedMaterials.size();
	 }

private:
	 //C++ STL map object that holds the named-material information
	 map<string,MaterialId> NamedMaterials;

	 bool MaterialTagExists(const string& material_tag);
};

#endif // CMATS_H

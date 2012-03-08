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

#ifndef CSHAPES_H
#define CSHAPES_H

//Declaration of the class "Cshapes" that holds the pointers to named (tagged) Cshape objects

#include "Cshapes_excp.h"

//only the declaration of Cshape needed: use forward declaration
class Cshape;

//for the C++ STL map class
#include <map>
////for the vector STL class
//#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

typedef boost::shared_ptr<Cshape> Cshape_ptr;


class Cshapes
{
 public:
	 void CreatePlanarLayer(const double& my_min_coord, const double& my_max_coord,const string& my_orientation, const string& shape_tag);
	 void CreateRectBox(const double& my_back_x, const double& my_front_x,
						   const double& my_left_y, const double& my_right_y,
						   const double& my_lower_z, const double& my_upper_z, const string& shape_tag);
	 void CreateSphere(const double& my_center_x, const double& my_center_y, const double& my_center_z, const double& my_radius, const string& shape_tag);

//	 bool lookupShapeWithTag(const string& shape_tag, const Cshape* shape_id) const;  //copies the shape identifier that corresponds to tag into mat_id, returns false if tag is not found

	 const Cshape* operator[] (const string& shape_tag) const; //returns the smart pointer to the shape corresponding to the shape tag, throws exception if not found

private:
	 //C++ STL map object that holds the named-shape information
	 map<string,Cshape_ptr> NamedShapes;

	 bool ShapeTagExists(const string& shape_tag);
};

#endif // CSHAPES_H

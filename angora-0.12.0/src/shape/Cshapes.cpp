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

//Definition of the class "Cshapes" that holds the pointers to all the Cshape (and derived) objects

#include "Cshapes.h"

#include "Cshape.h"


bool Cshapes::ShapeTagExists(const string& shape_tag)
{//returns true if a shape with the given tag exists.
	return (NamedShapes.find(shape_tag)!=NamedShapes.end());
}

//bool Cshapes::lookupShapeWithTag(const string& shape_tag, ShapeId& shape_id) const
//{//copies the identifier of the shape with tag into mat_id if it exists.
////returns false if the string tag is not found
//	map<string,ShapeId>::const_iterator map_it = NamedShapes.find(shape_tag);
//	if (map_it==NamedShapes.end())
//	{//string tag does not correspond to any shape, return false
//		return false;
//	}
//	else
//	{
//		shape_id = map_it->second;  //second element of the pair<string,ShapeId> object is the shape identifier
//		return true;
//	}
//}

const Cshape* Cshapes::operator[] (const string& shape_tag) const
{
	map<string,Cshape_ptr>::const_iterator map_it = NamedShapes.find(shape_tag);
	if (map_it==NamedShapes.end())
	{//string tag does not correspond to any shape, throw exception
		throw NamedShapeNotFoundException(shape_tag);
	}
	else
	{
		return map_it->second.get();
	}
}

void Cshapes::CreatePlanarLayer(const double& my_min_coord, const double& my_max_coord,const string& my_orientation, const string& shape_tag)
{//creates a planar layer with the given string tag
	if (!ShapeTagExists(shape_tag))
	{
		//create shape
		Cshape_ptr new_shape_ptr(new Cplanarlayer(my_min_coord,my_max_coord,my_orientation));
//		ShapePtrs.push_back(new_shape_ptr);
		//add the shape with the given tag
		NamedShapes.insert(pair<string,Cshape_ptr>(shape_tag,new_shape_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedShapeExistsException(shape_tag);
	}
}

void Cshapes::CreateRectBox(const double& my_back_x, const double& my_front_x,
						   const double& my_left_y, const double& my_right_y,
						   const double& my_lower_z, const double& my_upper_z, const string& shape_tag)
{//creates a rectangular box with the given string tag
	if (!ShapeTagExists(shape_tag))
	{
		//create shape
		Cshape_ptr new_shape_ptr(new Crectbox(my_back_x,my_front_x,my_left_y,my_right_y,my_lower_z,my_upper_z));
		//add the shape with the given tag
		NamedShapes.insert(pair<string,Cshape_ptr>(shape_tag,new_shape_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedShapeExistsException(shape_tag);
	}
}

void Cshapes::CreateSphere(const double& my_center_x, const double& my_center_y, const double& my_center_z, const double& my_radius, const string& shape_tag)
{//creates a sphere with the given string tag
	if (!ShapeTagExists(shape_tag))
	{
		//create shape
		Cshape_ptr new_shape_ptr(new Csphere(my_center_x,my_center_y,my_center_z,my_radius));
		//add the shape with the given tag
		NamedShapes.insert(pair<string,Cshape_ptr>(shape_tag,new_shape_ptr));
	}
	else
	{//string tag already exists, throw exception
		throw NamedShapeExistsException(shape_tag);
	}
}

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

#include "Cshape.h"

#include "float_comp.h"

namespace{
	const double tolerance_factor = 10;
	const double double_eps = tolerance_factor*LIBSTD_DBL_EPSILON;
}

extern double dx;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;


/**********************/
/** Cdirection class **/
/**********************/
bool Cdirection::IsParallel(const double& theta, const double& phi) const
{
	return approximatelyEqual(dot_product(theta,phi),1.0,double_eps);
}

bool Cdirection::IsParallel(const double& rx, const double& ry, const double& rz) const
{
	return approximatelyEqual(dot_product(rx,ry,rz),1.0,double_eps);
}

bool Cdirection::IsNormal(const double& theta, const double& phi) const
{
	return approximatelyEqual(dot_product(theta,phi),0.0,double_eps);
}

bool Cdirection::IsNormal(const double& rx, const double& ry, const double& rz) const
{
	return approximatelyEqual(dot_product(rx,ry,rz),0.0,double_eps);
}
//
///************************/
///** Cplanarsheet class **/
///************************/
//Cplanarsheet::Cplanarsheet(const double& my_coord, const string& my_orientation)
//		: coord(my_coord), orientation(my_orientation)
//{
//	if (orientation=="yz")
//	{
//		axis_of_symmetry = X;
//	}
//	else if (orientation=="xz")
//	{
//		axis_of_symmetry = Y;
//	}
//	else if (orientation=="xy")
//	{
//		axis_of_symmetry = Z;
//	}
//	else
//	{
//#ifdef __GNUG__
////GNU C++ compiler is being used, use the nice predefined variables for the function name
////		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
//		string func_name = __FUNCTION__;
//#else
//		string func_name = "";
//#endif
//		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,orientation,
//			"(valid arguments are \"yz\", \"xz\", or \"xy\")");
//	}
//};
//
//double Cplanarsheet::shortest_distance_to_surface(const double& x, const double& y, const double& z)
//{
//	/** TODO : Implement this later **/
//	return 0;
//	/** TODO : Implement this later **/
//}

/************************/
/** Cplanarlayer class **/
/************************/
Cplanarlayer::Cplanarlayer(const double& my_min_coord, const double& my_max_coord,const string& my_orientation)
		: min_coord(my_min_coord), max_coord(my_max_coord), orientation(my_orientation)
{
	if (orientation=="yz")
	{
		axis_of_symmetry = X;
		grid_min_coord_in_cells = (min_coord+_origin_x)/dx;
		grid_max_coord_in_cells = (max_coord+_origin_x)/dx;
	}
	else if (orientation=="xz")
	{
		axis_of_symmetry = Y;
		grid_min_coord_in_cells = (min_coord+_origin_y)/dx;
		grid_max_coord_in_cells = (max_coord+_origin_y)/dx;
	}
	else if (orientation=="xy")
	{
		axis_of_symmetry = Z;
		grid_min_coord_in_cells = (min_coord+_origin_z)/dx;
		grid_max_coord_in_cells = (max_coord+_origin_z)/dx;
	}
	else
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,orientation,
			"(valid arguments are \"yz\", \"xz\", or \"xy\")");
	}

	if (definitelyGreaterThan(min_coord,max_coord,double_eps))
	{
		throw AngoraInvalidBoundsException("planar layer",min_coord,max_coord);
	}
};

bool Cplanarlayer::IsInside(const double& x, const double& y, const double& z) const
{//returns true if the point with coordinate (x,y,z) with respect to the back-left-lower corner of the grid (in grid cells) is "inside" the layer (or within the lowest distance to the boundary representable by the double type)
	if (axis_of_symmetry==X)
	{
		return (lessThanOrEqual(x,grid_max_coord_in_cells,double_eps)
			  &&greaterThanOrEqual(x,grid_min_coord_in_cells,double_eps));
	}
	else if (axis_of_symmetry==Y)
	{
		return (lessThanOrEqual(y,grid_max_coord_in_cells,double_eps)
			  &&greaterThanOrEqual(y,grid_min_coord_in_cells,double_eps));
	}
	else if (axis_of_symmetry==Z)
	{
		return (lessThanOrEqual(z,grid_max_coord_in_cells,double_eps)
			  &&greaterThanOrEqual(z,grid_min_coord_in_cells,double_eps));
	}
}

bool Cplanarlayer::IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const
{//only works for planar layers that conform to the grid (i.e., perpendicular to either x,y, or z)
	if (axis_of_symmetry==X)
	{
		return ( (greaterThanOrEqual(x+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(x-0.5,grid_max_coord_in_cells,double_eps))
			   ||(lessThanOrEqual(x-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(x+0.5,grid_min_coord_in_cells,double_eps)));
	}
	else if (axis_of_symmetry==Y)
	{
		return ( (greaterThanOrEqual(y+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(y-0.5,grid_max_coord_in_cells,double_eps))
			   ||(lessThanOrEqual(y-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(y+0.5,grid_min_coord_in_cells,double_eps)));
	}
	else if (axis_of_symmetry==Z)
	{
		return ( (greaterThanOrEqual(z+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(z-0.5,grid_max_coord_in_cells,double_eps))
			   ||(lessThanOrEqual(z-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(z+0.5,grid_min_coord_in_cells,double_eps)));
	}
}

double Cplanarlayer::ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const
{
	if (axis_of_symmetry==X)
	{
		if (greaterThanOrEqual(x+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(x-0.5,grid_max_coord_in_cells,double_eps))
		{//upper boundary
			return (grid_max_coord_in_cells-max(grid_min_coord_in_cells,x-0.5))/dx; //the layer may be sub-cell, hence the max function
		}
		else if (lessThanOrEqual(x-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(x+0.5,grid_min_coord_in_cells,double_eps))
		{//lower boundary
			return (min(grid_max_coord_in_cells,x+0.5)-grid_min_coord_in_cells)/dx; //the layer may be sub-cell, hence the min function
		}
		else if (definitelyLessThan(x+0.5,grid_max_coord_in_cells,double_eps)
			      &&definitelyGreaterThan(x-0.5,grid_min_coord_in_cells,double_eps))
		{//completely inside
			return 1;
		}
		else
		{//completely outside
			return 0;
		}
	}
	else if (axis_of_symmetry==Y)
	{
		if (greaterThanOrEqual(y+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(y-0.5,grid_max_coord_in_cells,double_eps))
		{//upper boundary
			return (grid_max_coord_in_cells-max(grid_min_coord_in_cells,y-0.5))/dx; //the layer may be sub-cell, hence the max function
		}
		else if (lessThanOrEqual(y-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(y+0.5,grid_min_coord_in_cells,double_eps))
		{//lower boundary
			return (min(grid_max_coord_in_cells,y+0.5)-grid_min_coord_in_cells)/dx; //the layer may be sub-cell, hence the min function
		}
		else if (definitelyLessThan(y+0.5,grid_max_coord_in_cells,double_eps)
			      &&definitelyGreaterThan(y-0.5,grid_min_coord_in_cells,double_eps))
		{//completely inside
			return 1;
		}
		else
		{//completely outside
			return 0;
		}
	}
	else if (axis_of_symmetry==Z)
	{
		if (greaterThanOrEqual(z+0.5,grid_max_coord_in_cells,double_eps)
			      &&lessThanOrEqual(z-0.5,grid_max_coord_in_cells,double_eps))
		{//upper boundary
			return (grid_max_coord_in_cells-max(grid_min_coord_in_cells,z-0.5))/dx; //the layer may be sub-cell, hence the max function
		}
		else if (lessThanOrEqual(z-0.5,grid_min_coord_in_cells,double_eps)
			      &&greaterThanOrEqual(z+0.5,grid_min_coord_in_cells,double_eps))
		{//lower boundary
			return (min(grid_max_coord_in_cells,z+0.5)-grid_min_coord_in_cells)/dx; //the layer may be sub-cell, hence the min function
		}
		else if (definitelyLessThan(z+0.5,grid_max_coord_in_cells,double_eps)
			      &&definitelyGreaterThan(z-0.5,grid_min_coord_in_cells,double_eps))
		{//completely inside
			return 1;
		}
		else
		{//completely outside
			return 0;
		}
	}
}

int Cplanarlayer::bounding_box_back_cell() const
{
	if (orientation=="yz")
	{
		return (int)ceil(grid_min_coord_in_cells);
	}
	else if (orientation=="xz")
	{
		return 1;
	}
	else if (orientation=="xy")
	{
		return 1;
	}
}

int Cplanarlayer::bounding_box_front_cell() const
{
	if (orientation=="yz")
	{
		return (int)ceil(grid_max_coord_in_cells);
	}
	else if (orientation=="xz")
	{
		return NCELLS_X+2*NPML;
	}
	else if (orientation=="xy")
	{
		return NCELLS_X+2*NPML;
	}
}

int Cplanarlayer::bounding_box_left_cell() const
{
	if (orientation=="yz")
	{
		return 1;
	}
	else if (orientation=="xz")
	{
		return (int)ceil(grid_min_coord_in_cells);
	}
	else if (orientation=="xy")
	{
		return 1;
	}
}

int Cplanarlayer::bounding_box_right_cell() const
{
	if (orientation=="yz")
	{
		return NCELLS_Y+2*NPML;
	}
	else if (orientation=="xz")
	{
		return (int)ceil(grid_max_coord_in_cells);
	}
	else if (orientation=="xy")
	{
		return NCELLS_Y+2*NPML;
	}
}

int Cplanarlayer::bounding_box_lower_cell() const
{
	if (orientation=="yz")
	{
		return 1;
	}
	else if (orientation=="xz")
	{
		return 1;
	}
	else if (orientation=="xy")
	{
		return (int)ceil(grid_min_coord_in_cells);
	}
}

int Cplanarlayer::bounding_box_upper_cell() const
{
	if (orientation=="yz")
	{
		return NCELLS_Z+2*NPML;
	}
	else if (orientation=="xz")
	{
		return NCELLS_Z+2*NPML;
	}
	else if (orientation=="xy")
	{
		return (int)ceil(grid_max_coord_in_cells);
	}
}

/********************/
/** Crectbox class **/
/********************/
Crectbox::Crectbox(const double& my_back_x, const double& my_front_x,
		 const double& my_left_y, const double& my_right_y,
		 const double& my_lower_z, const double& my_upper_z)
		: back_x(my_back_x), front_x(my_front_x), left_y(my_left_y), right_y(my_right_y), lower_z(my_lower_z), upper_z(my_upper_z),
		  grid_back_x_in_cells((back_x+_origin_x)/dx), grid_front_x_in_cells((front_x+_origin_x)/dx), grid_left_y_in_cells((left_y+_origin_y)/dx), grid_right_y_in_cells((right_y+_origin_y)/dx), grid_lower_z_in_cells((lower_z+_origin_z)/dx), grid_upper_z_in_cells((upper_z+_origin_z)/dx)
//		  x_extent(front_x-back_x), y_extent(right_y-left_y), z_extent(upper_z-lower_z)
//		  center_x((front_x+back_x)/2), center_y((right_y+left_y)/2), center_z((upper_z+lower_z)/2)
{
	if (definitelyGreaterThan(back_x,front_x,double_eps))
	{
		throw AngoraInvalidBoundsException("rectangular box",back_x,front_x);
	}
	if (definitelyGreaterThan(left_y,right_y,double_eps))
	{
		throw AngoraInvalidBoundsException("rectangular box",left_y,right_y);
	}
	if (definitelyGreaterThan(lower_z,upper_z,double_eps))
	{
		throw AngoraInvalidBoundsException("rectangular box",upper_z,front_x);
	}
};

bool Crectbox::IsInside(const double& x, const double& y, const double& z) const
{//returns true if the point with coordinate (x,y,z) with respect to the back-left-lower corner of the grid (in grid cells) is "inside" the rectangular box (or within the lowest distance to the boundary representable by the double type)
	return definitelyLessThan(x,grid_front_x_in_cells,double_eps)&&definitelyGreaterThan(x,grid_back_x_in_cells,double_eps)
		 &&definitelyLessThan(y,grid_right_y_in_cells,double_eps)&&definitelyGreaterThan(y,grid_left_y_in_cells,double_eps)
		 &&definitelyLessThan(z,grid_upper_z_in_cells,double_eps)&&definitelyGreaterThan(z,grid_lower_z_in_cells,double_eps);
}

bool Crectbox::IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const
{
	/** TODO: Implement an efficient way to compute this **/
	return false;
	/** TODO: Implement an efficient way to compute this **/
}

double Crectbox::ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const
{
	/** TODO: Implement an efficient way to compute this **/
	return 0;
	/** TODO: Implement an efficient way to compute this **/
}

int Crectbox::bounding_box_back_cell() const
{
	return (int)ceil(grid_back_x_in_cells);
}

int Crectbox::bounding_box_front_cell() const
{
	return (int)ceil(grid_front_x_in_cells);
}

int Crectbox::bounding_box_left_cell() const
{
	return (int)ceil(grid_left_y_in_cells);
}

int Crectbox::bounding_box_right_cell() const
{
	return (int)ceil(grid_right_y_in_cells);
}

int Crectbox::bounding_box_lower_cell() const
{
	return (int)ceil(grid_lower_z_in_cells);
}

int Crectbox::bounding_box_upper_cell() const
{
	return (int)ceil(grid_upper_z_in_cells);
}

/*******************/
/** Csphere class **/
/*******************/
Csphere::Csphere(const double& my_center_x, const double& my_center_y, const double& my_center_z, const double& my_radius)
		: center_x(my_center_x), center_y(my_center_y), center_z(my_center_z), radius(my_radius),
		  grid_center_x_in_cells((center_x+_origin_x)/dx), grid_center_y_in_cells((center_y+_origin_y)/dx), grid_center_z_in_cells((center_z+_origin_z)/dx), radius_in_cells(radius/dx)
{
	if (definitelyGreaterThan(0.0,radius,double_eps))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<double>(func_name,radius,
			"(radius should be positive)");
	}
	radius_squared = radius*radius;
	radius_in_cells_squared = radius_in_cells*radius_in_cells;
};

bool Csphere::IsInside(const double& x, const double& y, const double& z) const
{//returns true if the point with coordinate (x,y,z) with respect to the back-left-lower corner of the grid (in grid cells) is "inside" the sphere (or within the lowest distance to the boundary representable by the double type)
	return (lessThanOrEqual((pow2(x-grid_center_x_in_cells)+pow2(y-grid_center_y_in_cells)+pow2(z-grid_center_z_in_cells)),radius_in_cells_squared,double_eps));
}

bool Csphere::IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const
{
	/** TODO: Implement an efficient way to compute this **/
	return false;
	/** TODO: Implement an efficient way to compute this **/
}

double Csphere::ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const
{
	/** TODO: Implement an efficient way to compute this **/
	return 0;
	/** TODO: Implement an efficient way to compute this **/
}

int Csphere::bounding_box_back_cell() const
{
	return (int)ceil(grid_center_x_in_cells-radius_in_cells);
}

int Csphere::bounding_box_front_cell() const
{
	return (int)ceil(grid_center_x_in_cells+radius_in_cells);
}

int Csphere::bounding_box_left_cell() const
{
	return (int)ceil(grid_center_y_in_cells-radius_in_cells);
}

int Csphere::bounding_box_right_cell() const
{
	return (int)ceil(grid_center_y_in_cells+radius_in_cells);
}

int Csphere::bounding_box_lower_cell() const
{
	return (int)ceil(grid_center_z_in_cells-radius_in_cells);
}

int Csphere::bounding_box_upper_cell() const
{
	return (int)ceil(grid_center_z_in_cells+radius_in_cells);
}

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

#ifndef CSHAPE_H
#define CSHAPE_H

#include "headers.h"

#include "Cshape_excp.h"

//for the vector STL class
#include <vector>

extern double dx;
extern int OriginX,OriginY,OriginZ;


class Cshape
{//base class for all kinds of shapes
	public:
		/** REMOVE LATER !!!! **/
		Cshape() : _origin_x((OriginX-1)*dx), _origin_y((OriginY-1)*dx), _origin_z((OriginZ-1)*dx) {};
		/** REMOVE LATER !!!! **/

		virtual ~Cshape() {};

		//is the point inside the volume? (pure virtual)
		virtual bool IsInside(const double& x, const double& y, const double& z) const =0;

		//does the boundary of the volume cut through a (dx,dx,dx) voxel centered at this position? (pure virtual)
		virtual bool IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const =0;
		//volumetric ratio of a (dx,dx,dx) voxel centered at this position that falls inside the volume [0 if completely outside, 1 if completely inside] (pure virtual)
		virtual double ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const =0;

		virtual int bounding_box_back_cell() const=0;
		virtual int bounding_box_front_cell() const=0;
		virtual int bounding_box_left_cell() const=0;
		virtual int bounding_box_right_cell() const=0;
		virtual int bounding_box_lower_cell() const=0;
		virtual int bounding_box_upper_cell() const=0;

	protected:
		/** REMOVE LATER !!!! **/
		//coordinates of the grid origin w.r.t. to the back-left-lower corner of the grid (in m)
		double _origin_x,_origin_y,_origin_z;
		/** REMOVE LATER !!!! **/
	private:
};

class Cdirection
{//class that represents a direction in space
	public:
		//ctor #1
		Cdirection(const double& theta, //in radians
				   const double& phi) //in radians
				: sx(sin(theta)*cos(phi)), sy(sin(theta)*sin(phi)), sz(cos(theta))
			{};
		//ctor #2
		Cdirection(const double& rx, const double& ry, const double& rz)
				: sx(rx/sqrt(pow2(rx)+pow2(ry)+pow2(rz))), sy(ry/sqrt(pow2(rx)+pow2(ry)+pow2(rz))), sz(rz/sqrt(pow2(rx)+pow2(ry)+pow2(rz)))
			{};
		virtual ~Cdirection() {};

		double dot_product(const double& theta, const double& phi) const
		{
			return (sx*sin(theta)*cos(phi)+sy*sin(theta)*sin(phi)+sz*cos(theta));
		}
		double dot_product(const double& rx, const double& ry, const double& rz) const
		{
			double r = sqrt(pow2(rx)+pow2(ry)+pow2(rz));
			return (sx*rx/r+sy*ry/r+sz*rz/r);
		}
		bool IsParallel(const double& theta, const double& phi) const;
		bool IsParallel(const double& rx, const double& ry, const double& rz) const;
		bool IsNormal(const double& theta, const double& phi) const;
		bool IsNormal(const double& rx, const double& ry, const double& rz) const;

	protected:
	private:
	 double sx,sy,sz;
};

//class Cvolume : public Cshape
//{//base class for all kinds of 3-D volumes
//	public:
////		Cvolume();
//		virtual ~Cvolume() {};
//
//		//does the boundary of the volume cut through a (dx,dx,dx) voxel centered at this position? (pure virtual)
//		virtual bool IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const =0;
//		//volumetric ratio of a (dx,dx,dx) voxel centered at this position that falls inside the volume [0 if completely outside, 1 if completely inside] (pure virtual)
//		virtual double ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const =0;
//
//	protected:
//	private:
//};

//class Csurface : public Cshape
//{//base class for all kinds of 2-D surfaces
//	public:
////		Csurface();
//		virtual ~Csurface() {};
//
//		//is the point inside the volume? (pure virtual)
//		virtual bool IsInside(const double& x, const double& y, const double& z) const
//		{
//			return false; //there is no "inside" of the surface
//		};
//
////		//does the surface cut through the cell with the given index? (pure virtual)
////		virtual bool CutsThroughCell(const int& m, const int& n, const int& p) const =0;
//		//shortest distance from a given point to the surface (pure virtual)
////		virtual double shortest_distance_to_surface(const double& x, const double& y, const double& z) const =0;
//
//	protected:
//	private:
//};
//
//class Cplanarsheet : public Csurface
//{//class that represents an infinitely-thin planar sheet that is oriented along the x,y, or z axes
//	public:
//		Cplanarsheet(const double& my_coord, const string& my_orientation);
//
//		virtual double shortest_distance_to_surface(const double& x, const double& y, const double& z);
//
//	private:
//	 enum {X=0,Y=1,Z=2}; //axis-of-invariance specifiers
//	 int axis_of_symmetry; //either 0,1, or 2
//	 const double coord;
//	 const string orientation;
//};

class Cplanarlayer : public Cshape
{//class that represents a planar layer oriented along one of the principal axes (x,y, or z)
	public:
		Cplanarlayer(const double& my_min_coord, const double& my_max_coord,const string& my_orientation);

		virtual bool IsInside(const double& x, const double& y, const double& z) const;
		virtual bool IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const;
		virtual double ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const;

		virtual int bounding_box_back_cell() const;
		virtual int bounding_box_front_cell() const;
		virtual int bounding_box_left_cell() const;
		virtual int bounding_box_right_cell() const;
		virtual int bounding_box_lower_cell() const;
		virtual int bounding_box_upper_cell() const;

	private:
	 enum {X=0,Y=1,Z=2}; //axis-of-invariance specifiers
	 int axis_of_symmetry; //either 0,1, or 2
	 const double min_coord; //(x, y or z) coordinate of the lower interface (w.r.t. the grid origin, in m)
	 const double max_coord; //(x, y or z) coordinate of the upper interface (w.r.t. the grid origin, in m)
	 const string orientation;
	 // coordinates w.r.t to the back-left-lower corner of the grid, in grid cells
	 double grid_min_coord_in_cells,grid_max_coord_in_cells;
};

class Crectbox : public Cshape
{//class that represents a rectangular box-shaped volume
	public:
		Crectbox(const double& my_back_x, const double& my_front_x,
				 const double& my_left_y, const double& my_right_y,
				 const double& my_lower_z, const double& my_upper_z);

		virtual bool IsInside(const double& x, const double& y, const double& z) const;
		virtual bool IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const;
		virtual double ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const;

		virtual int bounding_box_back_cell() const;
		virtual int bounding_box_front_cell() const;
		virtual int bounding_box_left_cell() const;
		virtual int bounding_box_right_cell() const;
		virtual int bounding_box_lower_cell() const;
		virtual int bounding_box_upper_cell() const;

	private:
	 //coordinates of the faces w.r.t. the grid origin, in m
	 const double back_x,front_x,left_y,right_y,lower_z,upper_z;
	 // coordinates of the faces w.r.t to the back-left-lower corner of the grid, in grid cells
	 const double grid_back_x_in_cells,grid_front_x_in_cells,grid_left_y_in_cells,grid_right_y_in_cells,grid_lower_z_in_cells,grid_upper_z_in_cells;
//	 //extents of the box in the x, y and z directions, in m
//	 const double x_extent,y_extent,z_extent;
//	 // extents of the box in grid cells
//	 const double grid_x_extent_in_cells,grid_y_extent_in_cells,grid_z_extent_in_cells;
//	 const double center_x,center_y,center_z;
};

class Csphere : public Cshape
{//class that represents a spherical volume
	public:
		Csphere(const double& my_center_x, const double& my_center_y, const double& my_center_z, const double& my_radius);

		virtual bool IsInside(const double& x, const double& y, const double& z) const;
		virtual bool IsAtBoundary(const double& x, const double& y, const double& z, const double& dx) const;
		virtual double ratio_of_voxel_inside_volume(const double& x, const double& y, const double& z, const double& dx) const;

		virtual int bounding_box_back_cell() const;
		virtual int bounding_box_front_cell() const;
		virtual int bounding_box_left_cell() const;
		virtual int bounding_box_right_cell() const;
		virtual int bounding_box_lower_cell() const;
		virtual int bounding_box_upper_cell() const;

	private:
	 //coordinates of the center w.r.t. the grid origin, in m
	 const double center_x,center_y,center_z;
	 // coordinates of the center w.r.t to the back-left-lower corner of the grid, in grid cells
	 const double grid_center_x_in_cells,grid_center_y_in_cells,grid_center_z_in_cells;
	 //radius in m
	 const double radius;
	 double radius_squared;
	 //radius in grid cells
	 const double radius_in_cells;
	 double radius_in_cells_squared;
};

#endif // CSHAPE_H

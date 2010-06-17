/*
  Copyright (C) 2006 - 2008  Jeffrey M. McMahon

  This file is part of JSCIENCE.

  JSCIENCE is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  JSCIENCE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with JSCIENCE.  If not, see <http://www.gnu.org/licenses/>.
*/


//************************************************************************
//------------------------------------------------------------------------
//
// NAME: JMATH_GEOMETRY
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_geometry.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 11/20/08
//
//	DESC:
//
//	NOTES: 	i. 
//
//
//	LIST OF SUBROUTINES:
//
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATH_GEOMETRY_H
#define JMATH_GEOMETRY_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include <cstdlib>
#include <cstring>

#include <fstream> 

#include <cmath>
#include <complex>

// CUSTOM FILES
#include "jmath_vector.h"

using namespace std;


//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************

typedef struct VERT3
{
  int v[3];
} vert3_t;

typedef struct VERT4
{
  int v[4];
} vert4_t;

// TETRAHEDRON
typedef struct TETRAHEDRON
{
  VECTOR_DOUBLE v[4];
} tet_t;

// TRIANGLE
// i. !!! note the bad use of 4 positions and not VERT3
typedef struct TRIANGLE
{
  // triangle vertices
  VECTOR_DOUBLE v[4];
} tri_t;

// 3D LINE
struct LINE_DOUBLE
{
  VECTOR_DOUBLE pt1_pos, pt2_pos;
};


//************************************************************************
//			EXTERNAL VARIABLES
//************************************************************************
extern int gauss_npts_2d[4], gauss_npts_3d[9];
extern double gauss_coord_1d[20][20], gauss_weight_1d[20][20], gauss_coord_2d[4][7][3], gauss_weight_2d[4][7], gauss_coord_3d[9][45][4], gauss_weight_3d[9][45];
extern int igauss_generated;

//************************************************************************
//			SUBROUTINES
//************************************************************************

// TETRAHEDRA
double get_volume(TETRAHEDRON tet_info);
void get_tet_info(TETRAHEDRON tet, double ae[4], double be[4], double ce[4], double de[4], double &volume);
void get_simplex_coord_tet(TETRAHEDRON tet, VECTOR_DOUBLE r, double zeta[4]);
void get_cart_coord_tet(TETRAHEDRON tet, double zeta[4], VECTOR_DOUBLE &r);
bool point_in_tet(tet_t tet, VECTOR_DOUBLE r);
void get_centroid_tet(TETRAHEDRON tet, VECTOR_DOUBLE &r);

// TRIANGLE
double get_area(TRIANGLE tri);
VECTOR_DOUBLE get_tri_nhat(TRIANGLE tri_info);
void get_tri_info(TRIANGLE tri, double ae[3], double be[3], double ce[3], double &area, VECTOR_DOUBLE &tri_nhat, VECTOR_DOUBLE &vec_u, VECTOR_DOUBLE &vec_v);
void get_simplex_coord_tri(TRIANGLE tri, VECTOR_DOUBLE r, double &zeta1, double &zeta2, double &zeta3);
void get_cart_coord_tri(TRIANGLE tri, double L[3], VECTOR_DOUBLE &r);
void get_projpt_tri(TRIANGLE tri, VECTOR_DOUBLE r, VECTOR_DOUBLE &r0);

double get_tri_circumcircle_radius(TRIANGLE tri);
void get_tri_circumcircle_center(TRIANGLE tri, VECTOR_DOUBLE &center);

// GAUSSIAN-LEGENDRE INTEGRATION
void generate_gauss();

#endif


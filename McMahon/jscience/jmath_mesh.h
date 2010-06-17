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
// NAME: JMESH: MAIN MESH CLASS FOR 3D OBJECTS
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_mesh.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 9/27/08
//
//------------------------------------------------------------------------
//************************************************************************
//
// NOTES:	i.
//
//
// TODO:	i. 
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATH_MESH_H	
#define JMATH_MESH_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include "jmath_geometry.h"

//************************************************************************
//			CLASSES / STRUCTS
//************************************************************************

// i. at first the structures vert3_t and vert4_t, as well as the definition of
// tetmesh_t may seem a little odd.  You may initially suspect the most natural
// definition of tetmesh_t would be in terms of tet_t and tri_t, but it is actually
// easiest to store all of the node positions and define the tetrahedra and
// triangles in terms of global node numbers


class TETRAHEDRAL_MESH
{
  public:
    int nnodes, ntets, ntris;
    VECTOR_DOUBLE *node;
    vert4_t *tet;
    vert3_t *tri;

    TETRAHEDRAL_MESH(int npts, int nele, int nface);
    ~TETRAHEDRAL_MESH(void);
};


typedef struct ELEMENT_TRI
{
  int attr;
  int v[3];
} element_tri_t;

typedef struct ELEMENT_TET
{
  int attr;
  int vert[4];
} element_tet_t;

typedef struct ELEMENT_LINE
{
  int attr;
  int n[2];
} element_line_t;

#endif	

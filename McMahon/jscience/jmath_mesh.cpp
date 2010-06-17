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
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_mesh.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 10/02/08
//
//
// TODO:	i. 
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

using namespace std;

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "jmath_mesh.h"


//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			SUBROUTINES
//************************************************************************

// i. starting at 1 is not a bad use of memory, but actually the way things are defined / created
TETRAHEDRAL_MESH::TETRAHEDRAL_MESH(int npts, int nele, int nface)
{
  nnodes = npts;
  ntets = nele;
  ntris = nface;

  node = new VECTOR_DOUBLE [nnodes*3];
  tet = new vert4_t[ntets+1]; 
  tri = new vert3_t[ntris+1];
}

TETRAHEDRAL_MESH::~TETRAHEDRAL_MESH(void)
{
  delete [] node;
  delete [] tet;
  delete [] tri;
}


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
// NAME: JSCIENCE VECTOR MATH SUBROUTINES
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_vector.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// last checked/updated on 11/22/08
//
//
//	DESC:	Geometric subroutines are for simple geometry, mainly for
//		Delaunay triangulations of 2D points.
//
//
//	LIST OF SUBROUTINES:
//
//
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATH_VECTOR_H
#define JMATH_VECTOR_H

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include <iostream>

#include <cmath>
#include <complex>

#include "jconsts.h"

using namespace std;

//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************

//=========================================================
// VECTOR_DOUBLE
//=========================================================
class VECTOR_DOUBLE
{ 

public:

  // CONSTRUCTORS
  VECTOR_DOUBLE();


  // operator overload

  VECTOR_DOUBLE operator+=(VECTOR_DOUBLE);
  VECTOR_DOUBLE operator+(VECTOR_DOUBLE);
//   COMPLEX_DOUBLE const operator+=(double const &d);
//   COMPLEX_DOUBLE const operator+=(COMPLEX_DOUBLE const &cd2);

  VECTOR_DOUBLE operator-=(VECTOR_DOUBLE);
  VECTOR_DOUBLE operator-(VECTOR_DOUBLE);
//   COMPLEX_DOUBLE const operator-=(double const &d);

//   COMPLEX_DOUBLE const operator*=(COMPLEX_DOUBLE const &cd2);
  VECTOR_DOUBLE operator*=(double);
  VECTOR_DOUBLE operator*(double);
  //VECTOR_DOUBLE operator*(complex<double>);

  VECTOR_DOUBLE operator/=(double);
  VECTOR_DOUBLE operator/(double);

  double x, y, z;
};


//=========================================================
// VECTOR_CD
//=========================================================
class VECTOR_CD
{ 

public:

  // CONSTRUCTORS
   VECTOR_CD();

   // operator overload
//   complex<double> const operator+=(double const &d);
//   complex<double> const operator+=(complex<double> const &cd2);

   VECTOR_CD operator+=(VECTOR_CD);
   VECTOR_CD operator+(VECTOR_CD);
//   complex<double> const operator-=(double const &d);

  //VECTOR_DOUBLE operator-=(VECTOR_DOUBLE);
  VECTOR_CD operator-(VECTOR_CD);

  //VECTOR_DOUBLE operator*=(double);
  VECTOR_CD operator*(double);
  VECTOR_CD operator*(complex<double>);

//   complex<double> const operator*=(double const &d);

//   complex<double> const operator/=(double const &d);
//   complex<double> const operator/=(complex<double> const &cd2);

   complex<double> x, y, z;
};


VECTOR_CD operator*(VECTOR_DOUBLE, complex<double>);

//************************************************************************
//			SUBROUTINES
//************************************************************************

// SINGLE VECTOR OPERATIONS
void zero(VECTOR_DOUBLE &vector);
void zero(VECTOR_CD &vector);
void normalize(VECTOR_DOUBLE &vector);
double get_length(VECTOR_DOUBLE vector);
VECTOR_CD conj(VECTOR_CD vector);

// SWAP VECTORS
void vector_swap(VECTOR_DOUBLE &vector1, VECTOR_DOUBLE &vector2);

// DOT PRODUCT OF VECTORS
double dot_product(VECTOR_DOUBLE vector1, VECTOR_DOUBLE vector2);
complex<double> dot_product(VECTOR_DOUBLE vector_d, VECTOR_CD vector_cd);
complex<double> dot_product(VECTOR_CD vector_cd, VECTOR_DOUBLE vector_d);
complex<double> dot_product(VECTOR_CD vector_cd1, VECTOR_CD vector_cd2);

// CROSS PRODUCT OF VECTORS
VECTOR_DOUBLE cross_product(VECTOR_DOUBLE vector1, VECTOR_DOUBLE vector2);
VECTOR_CD cross_product(VECTOR_CD vector1, VECTOR_DOUBLE vector2);
VECTOR_CD cross_product(VECTOR_DOUBLE vector1, VECTOR_CD vector2);
VECTOR_CD cross_product(VECTOR_CD vector1, VECTOR_CD vector2);



#endif


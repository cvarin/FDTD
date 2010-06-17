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
// last updated on 11/20/08
//
//	DESC:	Geometric subroutines are for simple geometry, mainly for
//		Delaunay triangulations of 2D points.
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



//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include "jmath_vector.h"


//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************

//=========================================================
// VECTOR_DOUBLE
//=========================================================


// default constructor
VECTOR_DOUBLE::VECTOR_DOUBLE() 
{
    x = y = z = 0.0;
}

//---------------------------------------------------------
// ADDITION
//---------------------------------------------------------
// Operator overloaded using a member function


/*
VECTOR_DOUBLE& VECTOR_DOUBLE::operator=(const VECTOR_DOUBLE &rhs) 
{
  if(this == &rhs)
  {
    return *this;
  }

  // ???

  return *this;
}
*/

VECTOR_DOUBLE VECTOR_DOUBLE::operator+=(VECTOR_DOUBLE rhs) 
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;

  return (*this);
}

VECTOR_DOUBLE VECTOR_DOUBLE::operator+(VECTOR_DOUBLE other) 
{
  VECTOR_DOUBLE temp;
  temp.x = x + other.x;
  temp.y = y + other.y;
  temp.z = z + other.z;

  return (temp);
}


//---------------------------------------------------------
// SUBTRACTION
//---------------------------------------------------------
VECTOR_DOUBLE VECTOR_DOUBLE::operator-=(VECTOR_DOUBLE rhs) 
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return (*this);
}


VECTOR_DOUBLE VECTOR_DOUBLE::operator-(VECTOR_DOUBLE other)
{
  VECTOR_DOUBLE temp;
  temp.x = x - other.x;
  temp.y = y - other.y;
  temp.z = z - other.z;

  return (temp);
}

VECTOR_CD VECTOR_CD::operator-(VECTOR_CD other)
{
  VECTOR_CD temp;
  temp.x = x - other.x;
  temp.y = y - other.y;
  temp.z = z - other.z;

  return (temp);
}



//---------------------------------------------------------
// MULTIPLICATION
//---------------------------------------------------------
VECTOR_DOUBLE VECTOR_DOUBLE::operator*=(double rhs) 
{
  x *= rhs;
  y *= rhs;
  z *= rhs;

  return (*this);
}

VECTOR_DOUBLE VECTOR_DOUBLE::operator*(double other)
{
  VECTOR_DOUBLE temp;
  temp.x = x*other;
  temp.y = y*other;
  temp.z = z*other;

  return (temp);
}

VECTOR_CD VECTOR_CD::operator*(double other)
{
  VECTOR_CD temp;
  temp.x = x*other;
  temp.y = y*other;
  temp.z = z*other;

  return (temp);
}

VECTOR_CD VECTOR_CD::operator*(complex<double> other)
{
  VECTOR_CD temp;
  temp.x = x*other;
  temp.y = y*other;
  temp.z = z*other;

  return (temp);
}

VECTOR_CD operator*(VECTOR_DOUBLE vector, complex<double> other)
{
  VECTOR_CD temp;
  temp.x = vector.x*other;
  temp.y = vector.y*other;
  temp.z = vector.z*other;

  return (temp);
}

//---------------------------------------------------------
// DIVISION
//---------------------------------------------------------
VECTOR_DOUBLE VECTOR_DOUBLE::operator/=(double rhs) 
{
  x /= rhs;
  y /= rhs;
  z /= rhs;

  return (*this);
}

VECTOR_DOUBLE VECTOR_DOUBLE::operator/(double other) 
{
  VECTOR_DOUBLE temp;
  temp.x = x/other;
  temp.y = y/other;
  temp.z = z/other;

  return (temp);
}

//=========================================================
// VECTOR_CD
//=========================================================

// default constructor
VECTOR_CD::VECTOR_CD() 
{
  x = y = z = COMPLEXZERO;
}


//---------------------------------------------------------
// ADDITION
//---------------------------------------------------------
VECTOR_CD VECTOR_CD::operator+=(VECTOR_CD vcd2) 
{
   x += vcd2.x;
   y += vcd2.y;
   z += vcd2.z;

   return (*this);
}

VECTOR_CD VECTOR_CD::operator+(VECTOR_CD other)
{
  VECTOR_CD temp;
  temp.x = x + other.x;
  temp.y = y + other.y;
  temp.z = z + other.z;

  return (temp);
}


//========================================================================
//========================================================================
//
//	NAME:	void vector_swap(VECTOR_DOUBLE &vector1, VECTOR_DOUBLE &vector2)
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
void vector_swap(VECTOR_DOUBLE &vector1, VECTOR_DOUBLE &vector2)
{

  double tempx, tempy, tempz;

  tempx = vector1.x;
  tempy = vector1.y;
  tempz = vector1.z;

  vector1.x = vector2.x;
  vector1.y = vector2.y;
  vector1.z = vector2.z;

  vector2.x = tempx;
  vector2.y = tempy;
  vector2.z = tempz;

  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
void normalize(VECTOR_DOUBLE &vector)
{

  double length = sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);

  vector.x /= length;
  vector.y /= length;
  vector.z /= length;

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
double get_length(VECTOR_DOUBLE vector)
{
  return sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
}

//========================================================================
//========================================================================
//
//	NAME:	void zero(VECTOR_DOUBLE &vector)
//		void zero(VECTOR_CD &vector)
//
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
void zero(VECTOR_DOUBLE &vector)
{
  vector.x = vector.y = vector.z = 0.0;

  return;
}

void zero(VECTOR_CD &vector)
{
  vector.x = vector.y = vector.z = COMPLEXZERO;

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
VECTOR_CD conj(VECTOR_CD vector)
{
  VECTOR_CD temp;
  temp.x = conj(vector.x); 
  temp.y = conj(vector.y); 
  temp.z = conj(vector.z); 

  return (temp);
}


//========================================================================
//========================================================================
//
//	NAME:	double dot_product(VECTOR_DOUBLE vector1, VECTOR_DOUBLE vector2)
//		complex<double> dot_product(VECTOR_CD vector_cd, VECTOR_DOUBLE vector_d)
//	DESC:	
//
//	NOTES: 	i. 
//
//
//========================================================================
//========================================================================
double dot_product(VECTOR_DOUBLE vector1, VECTOR_DOUBLE vector2)
{
  return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z;
}

complex<double> dot_product(VECTOR_DOUBLE vector_d, VECTOR_CD vector_cd)
{
  return vector_d.x*vector_cd.x + vector_d.y*vector_cd.y + vector_d.z*vector_cd.z;
}

complex<double> dot_product(VECTOR_CD vector_cd, VECTOR_DOUBLE vector_d)
{
  return vector_cd.x*vector_d.x + vector_cd.y*vector_d.y + vector_cd.z*vector_d.z;
}

complex<double> dot_product(VECTOR_CD vector_cd1, VECTOR_CD vector_cd2)
{
  return vector_cd1.x*vector_cd2.x + vector_cd1.y*vector_cd2.y + vector_cd1.z*vector_cd2.z;
}

//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	User-defined parameters
//
//	NOTES: 	i. recall that everything is in meters, seconds, etc.
//
//
//========================================================================
//========================================================================
VECTOR_DOUBLE cross_product(VECTOR_DOUBLE vector1, VECTOR_DOUBLE vector2)
{

  VECTOR_DOUBLE result;

  result.x = vector1.y*vector2.z - vector1.z*vector2.y;
  result.y = vector1.z*vector2.x - vector1.x*vector2.z;
  result.z = vector1.x*vector2.y - vector1.y*vector2.x;

  return result;
}

VECTOR_CD cross_product(VECTOR_CD vector1, VECTOR_DOUBLE vector2)
{

  VECTOR_CD result;

  result.x = vector1.y*vector2.z - vector1.z*vector2.y;
  result.y = vector1.z*vector2.x - vector1.x*vector2.z;
  result.z = vector1.x*vector2.y - vector1.y*vector2.x;

  return result;
}

VECTOR_CD cross_product(VECTOR_DOUBLE vector1, VECTOR_CD vector2)
{

  VECTOR_CD result;

  result.x = vector1.y*vector2.z - vector1.z*vector2.y;
  result.y = vector1.z*vector2.x - vector1.x*vector2.z;
  result.z = vector1.x*vector2.y - vector1.y*vector2.x;

  return result;
}

VECTOR_CD cross_product(VECTOR_CD vector1, VECTOR_CD vector2)
{

  VECTOR_CD result;

  result.x = vector1.y*vector2.z - vector1.z*vector2.y;
  result.y = vector1.z*vector2.x - vector1.x*vector2.z;
  result.z = vector1.x*vector2.y - vector1.y*vector2.x;

  return result;
}






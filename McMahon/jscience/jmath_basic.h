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
// NAME: JSCIENCE BASIC MATH SUBROUTINES
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_basic.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 12/18/08
//
//	DESC:	Basic math subroutines are for simple number comparisons,
//		number manipulations, and for determining the accuracy 
//		limits of a particular machine (by determining the round 
//		off unit for floating point arithmetic).
//
//
//	LIST OF SUBROUTINES:
//
//		NUMBER COMPARISON:
//
//			1) dmin(double d1, double d2) :: returns the minimum 
//				value of two real numbers
//			2) dmax(double d1, double d2) :: returns the maximum 
//				value of two real numbers
//			3) imin(int i1, int i2) :: returns the minimum value 
//				of two integers
//			4) imax(int i1, int i2) :: returns the maximum value 
//				of two integers
//			5) isign(int i) :: returns the sign of i as +1 or -1
//
//		NUMBER MANIPULATION:
//
//
//			7) int ifactorial(int n) :: returns the factorial of
//				n as an integer
//
//			8) double dfactorial(int n):: returns the factorial of
//				n as a real number
//
//		MACHINE ACCURACY DETERMINATION:
//
//			9) depsilon() :: returns the roundoff unit for floating 
//				point arithmetic.
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATH_BASIC_H
#define JMATH_BASIC_H



//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include <cstdlib>
#include <cmath>	
// screen I/O	
#include <iostream>

#include "jconsts.h"

using namespace std;	

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			SUBROUTINES
//************************************************************************

// NUMBER COMPARISON
double dmin(double d1, double d2);
double dmax(double d1, double d2);

int imax(int i1, int i2);
int imin(int i1, int i2);

int isign(int i);

// NUMBER MANIPULATION
int ifactorial(int n);
double dfactorial(int n);
double dtwofac(int n);

// MACHINE ACCURACY DETERMINATION
double depsilon(void);



#endif

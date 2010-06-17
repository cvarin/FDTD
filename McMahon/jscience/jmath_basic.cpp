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
// FILE: jmath_basic.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 04/19/08
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


//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include "jmath_basic.h"



//************************************************************************
//			SUBROUTINES
//************************************************************************



//========================================================================
//========================================================================
//
//	NAME:	double dmin(double d1, double d2)
//	DESC:	Returns the minimum value of two real numbers
//
//	INPUT:
//		double d1:: real number 1
//		double d2:: real number 2
//	
//	OUTPUT:
//		double d1, d2:: the minimum value of two real numbers
//		
//	NOTES:
//
//========================================================================
//========================================================================
double dmin(double d1, double d2)
{

  if(d1 < d2)
  {
    return d1;
  } 
  else
  {
    return d2;
  }

}


//========================================================================
//========================================================================
//
//	NAME:	double dmax(double d1, double d2)
//	DESC:	Returns the maximum value of two real numbers
//
//	INPUT:
//		double d1:: real number 1
//		double d2:: real number 2
//	
//	OUTPUT:
//		double d1, d2:: the maximum value of two real numbers
//		
//	NOTES:
//
//========================================================================
//========================================================================
double dmax(double d1, double d2)
{

  if(d1 > d2)
  {
    return d1;
  } 
  else
  {
    return d2;
  }

}


//========================================================================
//========================================================================
//
//	NAME:	int imin(int i1, int i2)
//	DESC:	Returns the minimum value of two integers
//
//	INPUT:
//		int i1:: integer 1
//		int i2:: integer 2
//	
//	OUTPUT:
//		int i1, i2:: the minimum value of two integers
//		
//	NOTES:
//
//========================================================================
//========================================================================
int imin(int i1, int i2)
{

  if(i1 < i2)
  {
    return i1;
  }
  else
  {
    return i2;
  }

}



//========================================================================
//========================================================================
//
//	NAME:	int imax(int i1, int i2)
//	DESC:	Returns the maximum value of two integers
//
//	INPUT:
//		int i1:: integer 1
//		int i2:: integer 2
//	
//	OUTPUT:
//		int i1, i2:: the maximum value of two integers
//
//	NOTES:
//
//========================================================================
//========================================================================
int imax(int i1, int i2)
{

  if(i1 > i2)
  {
    return i1;
  }
  else
  {
    return i2;
  }

}


//========================================================================
//========================================================================
//
//	NAME:	int isign(int i)
//	DESC:	Returns the sign of int i as +1 or -1
//
//	INPUT:
//		int i:: integer
//	
//	OUTPUT:
//		int sign:: +1 for i>=0, -1 for i<0
//
//	NOTES:	1) The sign of 0 is taken to be +
//
//========================================================================
//========================================================================
int isign(int i)
{

  if(i < 0) 
  {
    return -1;
  }
  else
  {
    return 1;
  }

}


//========================================================================
//========================================================================
//
//	NAME:	double depsilon(void)
//
//	DESC:	Returns the roundoff unit for floating point arithmetic
//
//	INPUT:
//	
//	OUTPUT:
//		double epsilon:: the roundoff unit for floating point arithmetic
//
//	NOTES:
//		1) The algorithm used is that the round-off unit is a number, 
//		r, such that:
//			1.0 < (1.0 + r)
//			1.0 = (1.0 + r/2.0)
//
//		ii. ! I saw this routine somewhere, but I don't have the reference for it
//
//========================================================================
//========================================================================
double depsilon(void)
{

  double epsilon;

  double r = 1.0;

  while(1.0 < (1.0 + r))
  {
    r = r/2.0;
  }

  epsilon = 2.0*r;

  return epsilon;

}



//========================================================================
//========================================================================
//
//	NAME:	int ifactorial(int n)
//	DESC:	Calculates the factorial of n and returns the value as an 
//		integer
//
//	INPUT:
//		int n == integer
//
//	OUTPUT:
//		int factorial == n!
//
//	NOTES:	1) This will only works to a certain max n (max int)
//
//
//========================================================================
//========================================================================
int ifactorial(int n)
{
  return static_cast<int>(dfactorial(n));
}


//========================================================================
//========================================================================
//
//	NAME:	double dfactorial(int n)
//	DESC:	Calculates the factorial of n and returns the value as a 
//		real number
//
//	INPUT:
//		int n == integer
//
//	OUTPUT:
//		double factorial == double(n!)
//
//	NOTES:	1) !!!
//		2) There is no support is factorial > double
//
//
//========================================================================
//========================================================================
double dfactorial(int n)
{

  double factorial = 1.0;

  if(n < 2)
  {
    factorial = 1.0;
  }
  else if(n > 30)
  {
    // !!!
    double dn = static_cast<double>(n);
    factorial = sqrt(2.0*PI)*pow(dn, n)*sqrt(dn)*exp(-dn);
  }
  else
  {
    for(int i = n; i >= 2; --i)
    {
      factorial = factorial*static_cast<double>(i);
    }
  }


  return factorial;
}


//========================================================================
//========================================================================
//
//	NAME:	double dtwofac(int l)
//	DESC:	Calculates the double factorial (2l-1)!!=1*3*5*...*(2l-1) for
//		argument l, which can be defined as:
//
//			(2l-1)!!=(2l)!/((2^n)*l!)
//
//	NOTES: 	i. !!! there is a better way to define this with the gamma
//		function
//
//
//========================================================================
//========================================================================
double dtwofac(int n)
{
  return (  dfactorial(2*n)/(   pow(2.0, n)*dfactorial(n)  )  );
}


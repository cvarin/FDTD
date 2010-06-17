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
// NAME: JSCIENCE SPECIAL FUNCTION SUBROUTINES
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_special.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
//	DESC:	Special subroutines are for calculating Bessel functions
//		of the first and second kinds and Hankel functions.
//
//
//	LIST OF SUBROUTINES:
//
//		BESSEL FUNCTIONS:
//
//			1) double bessel_1(int n, double arg) :: returns the  
//				Bessel function of the first kind of order n
//			2) complex<double> bessel_1_complex(int n, 
//				complex<double> arg) :: returns the Bessel 
//				function of the first kind of order n for 
//				complex arguments
//			3) double bessel_2(int n, double arg) :: returns the  
//				Bessel function of the second kind of order n
//			4) complex<double> bessel_2_complex(int n, 
//				complex<double> arg) :: returns the Bessel 
//				function of the second kind of order n for 
//				complex arguments
//			5) double h_s(int s) :: returns (1/1 + 1/2 + ... + 1/s)
//
//		BESSEL FUNCTION DERIVATIVES:
//
//			6) complex<double> bessel_1_deriv_complex(int n, 
//			 	complex<double> arg) :: returns the derivative 
//				of the first kind of order n for a complex
//				argument
//
//		HANKEL FUNCTIONS:
//
//			7) complex<double> hankel_1(int n, complex<double> arg) 
//				:: returns the Hankel function of the first kind 
//				of order n for complex arguments
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATH_SPECIAL_H
#define JMATH_SPECIAL_H



//************************************************************************
//			INCLUDE FILES
//************************************************************************
#include <cmath>
#include <complex>

#include "jconsts.h"
#include "jmath_basic.h"

using namespace std;	


//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			SUBROUTINES
//************************************************************************

// BESSEL FUNCTIONS
double bessel_1(int n, double arg);
complex<double> bessel_1_complex(int n, complex<double> arg);
double bessel_2(int n, double arg);
complex<double> bessel_2_complex(int n, complex<double> arg);
double h_s(int s);

// BESSEL FUNCTION DERIVATIVES
complex<double> bessel_1_deriv_complex(int n, complex<double> arg);

// LINEAR COMBINATIONS OF BESSEL FUNCTIONS
complex<double> hankel_1(int n, double arg);


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// BESSEL FUNCTIONS
void get_hn(int n, int m, double z, complex<double>& hn);

void get_g(double k0, double xpos, double xposp, complex<double>& g);
double digamma(int m);



#endif

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
// NAME: JSCIENCE CONSTANTS: STANDARD/PHYSICAL CONSTANTS FOR JPROGRAMS
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jconsts.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 7/22/08
//
// TODO:	i. !!! maybe separate standard, complex, and physical constants
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JCONSTS_H	
#define JCONSTS_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************

// STANDARD
// math
#include <cmath>
#include <complex>
// !!! temp
#include <iostream>

using namespace std;

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


extern const double PI, EULERC, NM;
extern const complex<double> COMPLEXZERO, COMPLEXONE, COMPLEXJ;
extern const double CSPEED, EPS0, EPS0R, MU0, MU0R, Z0, PLANCK, PLANCKEV, a0, HBAR, ECHARGE, EMASS;


//************************************************************************
//			SUBROUTINES
//************************************************************************

#endif	

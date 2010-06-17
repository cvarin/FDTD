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
// FILE: jconsts.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 7/22/08
//
//
// TODO:	i. !!! maybe separate standard, complex, and physical constants
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

using namespace std;

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "jconsts.h"

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


// STANDARD CONSTANTS
// calculate pi using Machin's formula (good to ~500 digits)
const double PI=16.0*atan(0.2)-4.0*atan(1.0/239.0);
const double EULERC = 0.57721566490153286060651209008;
const double NM=1.0e-9;

// COMPLEX CONSTANTS
const complex<double> COMPLEXZERO(0.0,0.0);
const complex<double> COMPLEXONE(1.0,0.0);
const complex<double> COMPLEXJ(0.0,1.0);

// PHYSICAL CONSTANTS
const double CSPEED=2.99792458e8;
const double MU0=PI*4.0e-7;
const double EPS0=(1.0e7)/(4.0*PI*CSPEED*CSPEED);
const double MU0R=1.0;
const double EPS0R=1.0;
const double Z0=sqrt(MU0/EPS0);
const double PLANCK=6.62606896e-34;
const double PLANCKEV=4.13566733e-15;

const double a0=5.2917720859e-11;
const double HBAR=1.054571629e-34;
const double ECHARGE=1.602176487e-19;
const double EMASS=9.10938215e-31;

//************************************************************************
//			SUBROUTINES
//************************************************************************


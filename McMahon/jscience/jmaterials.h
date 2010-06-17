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
// NAME: JSCIENCE DRUDE + 2 LORENTZ POLE MODEL FOR MATERIALS
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmaterials.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 9/13/08
//
//	DESC:
//
//	NOTES: 	i. remember my universal convention:
//
// 			i.a. types:
// 				99 == Free space
//				90-98 (0-8) == User defined (program specific)
// 				0 == Gold (D+2L)
//				1 == Silver (D+2L)
//
//	LIST OF SUBROUTINES:
//
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************


#ifndef JMATERIALS_H
#define JMATERIALS_H	

using namespace std;	

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include <cstdlib>
#include <cmath>
#include <complex>

#include <iostream>

#include "jconsts.h"		

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************

extern const int d2l_libsize;
extern int library_generated;
extern double d2l_inf[3], d2l_conduct[3], d2l_drude_omegap[3], d2l_drude_gammap[3], d2l_lorentz_omegap[3][2], d2l_lorentz_deltap[3][2], d2l_lorentz_depsr[3][2], d2l_lorentz_alphap[3][2], d2l_lorentz_zetap[3][2], d2l_lorentz_gammap[3][2];
extern double nl_dh_vf[1], nl_dh_inf[1], nl_dh_omegap[1], nl_dh_gammap[1];
extern complex<double> epsr_custom[9];


//************************************************************************
//			SUBROUTINES
//************************************************************************


void generate_epsrlib();

complex<double> get_epsr(int imaterial, double wavelength, int exp_convention);


#endif

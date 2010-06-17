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
// NAME: JFEM_CONFIG: CONFIGURATION FILE FOR JFEM3D
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jfem_config.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 12/20/08
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

#ifndef JPLANEWAVE_H	
#define JPLANEWAVE_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "jconsts.h"

#include "jmath_vector.h"

//************************************************************************
//			CLASSES / STRUCTS
//************************************************************************

void get_pw_E(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &Einc);
void get_pw_dE(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &dxEinc, VECTOR_CD &dyEinc, VECTOR_CD &dzEinc);
void get_pw_H(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &Hinc);
void get_pw_dH(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &dxHinc, VECTOR_CD &dyHinc, VECTOR_CD &dzHinc);

#endif	

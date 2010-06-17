/*
  Copyright (C) 2006 - 2008  Jeffrey M. McMahon

  This file is part of JFDTD2D.

  JFDTD2D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  JFDTD2D is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with JFDTD2D.  If not, see <http://www.gnu.org/licenses/>.
*/


//************************************************************************
//------------------------------------------------------------------------
//
// NAME: STRUCTURE SETUP: SUBROUTINE TO SET UP STRUCTURE
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jfdtd2d.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 1/7/09
//
// NOTES:	i. ! the flag for an NSOM probe is defined and used in jfdtd2d.cpp
//			
// TODO:	i. 
//
// CODES:	i.
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************


//************************************************************************
//			INCLUDE FILES
//************************************************************************	

// STANDARD	
// utilities
#include <cstdlib>
// math
#include <cmath>

// CUSTOM
// Drude + 2 Lorentz pole modes
#include "jmaterials.h"
// standard, complex, and physical constants
#include "jconsts.h"

using namespace std;


//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************

//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************

//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

//========================================================================
//========================================================================
//
// NAME: void get_structure(double xpos, double ypos, double zpos, int& mat_num, 
//          double& mat_eps, double& mat_mu)
//
// DESC: Allows the user to define a custom structure.
//
// INST: To define a structure the user must provide an algorithm here.  The
//          idea is, given a position (xpos, ypos, zpos) assign a material number (99
//          for a real dielectric or any other number <90 for a frequency
//          dependent material consistent with jmaterials.cpp in JSCIENCE).  The 
//          user must also specify the permittivity and permeability.  For a 
//          frequency dependent material the permittivity should be set to
//          EPS0*d2l_inf[mat_num].
//
//========================================================================
//========================================================================
void get_structure(double xpos, double ypos, double zpos, int& mat_num, double& mat_eps, double& mat_mu)
{
  // SET DEFAULT FREE SPACE VALUES
  mat_num = 99;
  mat_eps = EPS0;
  mat_mu = MU0;

  //=========================================================
  // DEFINE STRUCTURE BELOW
  //=========================================================


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================

  return;
}



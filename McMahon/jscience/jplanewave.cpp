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
// *last updated on 12/31/08
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

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "jplanewave.h"

//************************************************************************
//			CLASSES / STRUCTS
//************************************************************************


//************************************************************************
//			FULL CODE GLOBAL VARIABLES
//************************************************************************

//========================================================================
//========================================================================
//
//	NAME:	void getEinc(VECTOR_DOUBLE pos, VECTOR_CD &Einc)
//	DESC:	Get the initial field (Einc) at the point (pos)
//
//	NOTES:	i. the field is defined as, which was taken from pg. 451 in Jin's FEM book:
//			
//			E=(cos(alpha)*theta^ + sin(alpha)*psi^)*exp(-j*kinc*r)
//
//			ii.a. the direction is given by:
//
//				kinc=-k0*(sin(theta)*cos(psi)*x^ + sin(theta)*sin(psi)*y^ + cos(theta)*z^)
//
// *******************************************
//				!!! upon calculation of the power flow in (calc_scs2), it appears that:
//
//				kinc=k0*(sin(theta)*cos(psi)*x^ + sin(theta)*sin(psi)*y^ + cos(theta)*z^)
//
//				which would give:
//
//				E=(cos(alpha)*theta^ + sin(alpha)*psi^)*exp(j*kinc*r)
// *******************************************
//
//			ii.b. the polarization is given by (which depends on Kinc):
//
//				theta^=( cos(psi)*cos(theta), sin(psi)*cos(theta), -sin(theta) )
//				psi^=( -sin(psi), cos(psi), 0 )
//
//			ii.c. e.g. (theta=0, psi=pi/2, alpha=0)-->-z directed and y polarized
//
//		ii. the equations were translated from the mathematical notation at:
//			http://mathworld.wolfram.com/SphericalCoordinates.html
//
//========================================================================
//========================================================================
void get_pw_E(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &Einc)
{
  // GET THE PHASE
  complex<double> phase = exp( COMPLEXJ*k0*(sin(theta)*cos(psi)*pos.x + sin(theta)*sin(psi)*pos.y + cos(theta)*pos.z) );

  // GET THE POLARIZATION
  VECTOR_DOUBLE pol;
  pol.x = cos(alpha)*cos(psi)*cos(theta) + sin(alpha)*(-sin(psi));
  pol.y = cos(alpha)*sin(psi)*cos(theta) + sin(alpha)*cos(psi);
  pol.z = cos(alpha)*(-sin(theta));

  // GET THE INCIDENT FIELD
  Einc = pol*E0*phase;

  return;
}




//========================================================================
//========================================================================
//
//	NAME:	void getEinc(VECTOR_DOUBLE pos, VECTOR_CD &Einc)
//	DESC:	Get the initial field (Einc) at the point (pos)
//
//	NOTES:	i. the field is defined as, which was taken from pg. 451 in Jin's FEM book:
//			
//			E=(cos(alpha)*theta^ + sin(alpha)*psi^)*exp(-j*kinc*r)
//
//			ii.a. the direction is given by:
//
//				kinc=-k0*(sin(theta)*cos(psi)*x^ + sin(theta)*sin(psi)*y^ + cos(theta)*z^)
//
//			ii.b. the polarization is given by (which depends on Kinc):
//
//				theta^=( cos(psi)*cos(theta), sin(psi)*cos(theta), -sin(theta) )
//				psi^=( -sin(psi), cos(psi), 0 )
//
//			ii.c. e.g. (theta=0, psi=pi/2, alpha=0)-->-z directed and y polarized
//
//		ii. the equations were translated from the mathematical notation at:
//			http://mathworld.wolfram.com/SphericalCoordinates.html
//
//========================================================================
//========================================================================
void get_pw_dE(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &dxEinc, VECTOR_CD &dyEinc, VECTOR_CD &dzEinc)
{
  VECTOR_CD Einc;
  complex<double> dx, dy, dz;

  // GET THE DERIVATIVE FACTORS
  dx = COMPLEXJ*k0*sin(theta)*cos(psi);
  dy = COMPLEXJ*k0*sin(theta)*sin(psi);
  dz = COMPLEXJ*k0*cos(theta);

  // GET THE INCIDENT FIELD
  get_pw_E(pos, E0, k0, theta, psi, alpha, Einc);

  dxEinc.x = dx*Einc.x;
  dyEinc.x = dy*Einc.x;
  dzEinc.x = dz*Einc.x;

  dxEinc.y = dx*Einc.y;
  dyEinc.y = dy*Einc.y;
  dzEinc.y = dz*Einc.y;

  dxEinc.z = dx*Einc.z;
  dyEinc.z = dy*Einc.z;
  dzEinc.z = dz*Einc.z;

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void getEinc(VECTOR_DOUBLE pos, VECTOR_CD &Einc)
//	DESC:	Get the initial field (Einc) at the point (pos)
//
//	NOTES:	i. the field is defined as, which was taken from pg. 451 in Jin's FEM book:
//			
//			E=(cos(alpha)*theta^ + sin(alpha)*psi^)*exp(-j*kinc*r)
//
//			ii.a. the direction is given by:
//
//				kinc=-k0*(sin(theta)*cos(psi)*x^ + sin(theta)*sin(psi)*y^ + cos(theta)*z^)
//
//			ii.b. the polarization is given by (which depends on Kinc):
//
//				theta^=( cos(psi)*cos(theta), sin(psi)*cos(theta), -sin(theta) )
//				psi^=( -sin(psi), cos(psi), 0 )
//
//			ii.c. e.g. (theta=0, psi=pi/2, alpha=0)-->-z directed and y polarized
//
//		ii. the equations were translated from the mathematical notation at:
//			http://mathworld.wolfram.com/SphericalCoordinates.html
//
//========================================================================
//========================================================================
void get_pw_H(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &Hinc)
{

  VECTOR_CD dxEinc, dyEinc, dzEinc, curlEinc;

  // GET THE DERIVATIVES OF Einc
  get_pw_dE(pos, E0, k0, theta, psi, alpha, dxEinc, dyEinc, dzEinc);

  // GET curl(Einc)
  curlEinc.x = dyEinc.z - dzEinc.y;
  curlEinc.y = dzEinc.x - dxEinc.z;
  curlEinc.z = dxEinc.y - dyEinc.x;

  // NOW CALCULATE Hinc
  Hinc = curlEinc*( -1.0/(COMPLEXJ*MU0*k0*CSPEED) );

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void getEinc(VECTOR_DOUBLE pos, VECTOR_CD &Einc)
//	DESC:	Get the initial field (Einc) at the point (pos)
//
//	NOTES:	i. the field is defined as, which was taken from pg. 451 in Jin's FEM book:
//			
//			E=(cos(alpha)*theta^ + sin(alpha)*psi^)*exp(-j*kinc*r)
//
//			ii.a. the direction is given by:
//
//				kinc=-k0*(sin(theta)*cos(psi)*x^ + sin(theta)*sin(psi)*y^ + cos(theta)*z^)
//
//			ii.b. the polarization is given by (which depends on Kinc):
//
//				theta^=( cos(psi)*cos(theta), sin(psi)*cos(theta), -sin(theta) )
//				psi^=( -sin(psi), cos(psi), 0 )
//
//			ii.c. e.g. (theta=0, psi=pi/2, alpha=0)-->-z directed and y polarized
//
//		ii. the equations were translated from the mathematical notation at:
//			http://mathworld.wolfram.com/SphericalCoordinates.html
//
//========================================================================
//========================================================================
void get_pw_dH(VECTOR_DOUBLE pos, double E0, double k0, double theta, double psi, double alpha, VECTOR_CD &dxHinc, VECTOR_CD &dyHinc, VECTOR_CD &dzHinc)
{
  VECTOR_CD Hinc;
  complex<double> dx, dy, dz;

  // GET THE DERIVATIVE FACTORS
  dx = COMPLEXJ*k0*sin(theta)*cos(psi);
  dy = COMPLEXJ*k0*sin(theta)*sin(psi);
  dz = COMPLEXJ*k0*cos(theta);

  // GET THE INCIDENT FIELD
  get_pw_H(pos, E0, k0, theta, psi, alpha, Hinc);

  dxHinc.x = dx*Hinc.x;
  dyHinc.x = dy*Hinc.x;
  dzHinc.x = dz*Hinc.x;

  dxHinc.y = dx*Hinc.y;
  dyHinc.y = dy*Hinc.y;
  dzHinc.y = dz*Hinc.y;

  dxHinc.z = dx*Hinc.z;
  dyHinc.z = dy*Hinc.z;
  dzHinc.z = dz*Hinc.z;

  return;
}




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
// FILE: jmaterials.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 9/13/08
//
//	DESC:
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
#include "jmaterials.h"


//************************************************************************
//			SUBROUTINES
//************************************************************************


const int d2l_libsize = 3;

// DRUDE + 2 LORENTZ POLE MODEL
double d2l_inf[3];
double d2l_conduct[3];
double d2l_drude_omegap[3];
double d2l_drude_gammap[3];
double d2l_lorentz_omegap[3][2];
double d2l_lorentz_deltap[3][2];
double d2l_lorentz_depsr[3][2];
double d2l_lorentz_alphap[3][2];
double d2l_lorentz_zetap[3][2];
double d2l_lorentz_gammap[3][2];

// NON-LOCAL HYRDODYNAMIC DRUDE MODEL
double nl_dh_vf[1];
double nl_dh_inf[1];
double nl_dh_omegap[1];
double nl_dh_gammap[1];

// CUSTOM DIELECTRIC
complex<double> epsr_custom[9];

int library_generated = 0;


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
void generate_epsrlib()
{

  // i. Palik Al2O3 is approximately n=1.77 over 250-1000 nm  

/*
// *******************************************
  // METAL INFINITE FREQUENCY
  //d2l_inf[1] = 5.398334980;
  d2l_inf[1] = 1.0;

  // CONDUCTIVITY
  d2l_conduct[1] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[1] = 1.36e16;
  // damping constant
  d2l_drude_gammap[1] = 2.56e13 + 1.0*(1.39e6)/(3.0*NM - 2.5*NM);

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[1][0] = 0.0;
  d2l_lorentz_depsr[1][1] = 0.0;
  // plasmon frequencies
  d2l_lorentz_omegap[1][0] = 5.43497*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[1][1] = 2.0122*(2.0*PI/PLANCKEV);
  // plasmon frequencies
  d2l_lorentz_deltap[1][0] = 0.128755*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[1][1] = 0.56124*(2.0*PI/PLANCKEV);
// ********************************************
*/

  //=========================================================
  // Au
  //=========================================================

  //---------------------------------------------------------
  // Au - Lynch & Hunter (<1000nm)
  //---------------------------------------------------------
  // i. this is taken from JMM, et. al. Optics Express, 15, 18119 (2007)

  // METAL INFINITE FREQUENCY
  d2l_inf[0] = 5.398334980;

  // CONDUCTIVITY
  d2l_conduct[0] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[0] = 0.1397832433e17;
  // damping constant
  d2l_drude_gammap[0] = 0.1033372319e15;

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[0][0] = 2.541747093*0.2678716779;
  d2l_lorentz_depsr[0][1] = 2.541747093*(1.0 - 0.2678716779);
  // plasmon frequencies
  d2l_lorentz_omegap[0][0] = 0.4273918467e16;
  d2l_lorentz_omegap[0][1] = 0.5225397872e16;
  // plasmon frequencies 
  d2l_lorentz_deltap[0][0] = 0.4353347310e15;
  d2l_lorentz_deltap[0][1] = 0.6607732694e15;



  //=========================================================
  // Ag
  //=========================================================



  //---------------------------------------------------------
  // Ag METAL - Johnson & Christy fit (good for 300-800 nm)
  //---------------------------------------------------------

  // METAL INFINITE FREQUENCY
  d2l_inf[1] = 1.17152;

  // CONDUCTIVITY
  d2l_conduct[1] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[1] = sqrt(84.4362)*(2.0*PI/PLANCKEV);
  // damping constant
  d2l_drude_gammap[1] = (0.830175e-14)*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[1][0] = 2.23994;
  d2l_lorentz_depsr[1][1] = 0.222651;
  // plasmon frequencies
  d2l_lorentz_omegap[1][0] = 5.43497*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[1][1] = 2.0122*(2.0*PI/PLANCKEV);
  // plasmon frequencies
  d2l_lorentz_deltap[1][0] = 0.128755*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[1][1] = 0.56124*(2.0*PI/PLANCKEV);



/*
  //---------------------------------------------------------
  // Ag METAL - Palik (Lynch & Hunter) (good for ~250-1000 nm)
  //---------------------------------------------------------
  // i. this is taken from the TWL/SKG Optics Express, 13, 9652 (2005)

  // METAL INFINITE FREQUENCY
  d2l_inf[1] = 2.3646;

  // CONDUCTIVITY
  d2l_conduct[1] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[1] = 8.7377*(2.0*PI/PLANCKEV);
  // damping constant
  d2l_drude_gammap[1] = 0.07489*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[1][0] = 1.1831*0.2663;
  d2l_lorentz_depsr[1][1] = 1.1831*0.7337;
  // plasmon frequencies
  d2l_lorentz_omegap[1][0] = 4.3802*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[1][1] = 5.183*(2.0*PI/PLANCKEV);
  // plasmon frequencies
  d2l_lorentz_deltap[1][0] = 0.28*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[1][1] = 0.5482*(2.0*PI/PLANCKEV);
*/


  //=========================================================
  // Pd
  //=========================================================

  //---------------------------------------------------------
  // Pd - Johnson & Christy (300 - 1000 nm)
  //---------------------------------------------------------
  // i. Phys. Rev. B 9 (1974) 5056 - 5070

  // METAL INFINITE FREQUENCY
  d2l_inf[2] = 1.0000000000000018;
 
  // CONDUCTIVITY
  d2l_conduct[2] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[2] = 6.6813657338967589*(2.0*PI/PLANCKEV);
  // damping constant
  d2l_drude_gammap[2] = (4.47647945068469009E-017)*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[2][0] = 120.77863136686234*(2.26163049981847088E-002);
  d2l_lorentz_depsr[2][1] = 120.77863136686234*(1.0 - 2.26163049981847088E-002);
  // plasmon frequencies
  d2l_lorentz_omegap[2][0] = 6.1288546993501809*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[2][1] = 0.93595745548215437*(2.0*PI/PLANCKEV);
  // plasmon frequencies 
  d2l_lorentz_deltap[2][0] = 3.4185299135962706*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[2][1] = 1.4829360456231302*(2.0*PI/PLANCKEV);

/*
  //=========================================================
  // ITO
  //=========================================================


  //---------------------------------------------------------
  // ITO 1
  //---------------------------------------------------------
  // i.) this is taken from M. Losurdo paper J. Vac. Sci. Tech
  // on the parameterization of ITO

  // infinite frequency epsilon
  d2l_inf[1] = 1.0;

  // conductivity
  d2l_conduct[1] = 0.0;

  // DRUDE MODEL

  // plasmon frequency (eV conversion)
  d2l_drude_omegap[1] = 1.62*(2.0*PI/PLANCKEV);
  // damping constant (eV conversion)
  d2l_drude_gammap[1] = 0.05*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL

  // assign epsilon widths
  d2l_lorentz_depsr[1][0] = 0.47;
  d2l_lorentz_depsr[1][1] = 2.3;

  // plasmon frequencies (eV conversions)
  d2l_lorentz_omegap[1][0] = 4.82*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[1][1] = 7.1*(2.0*PI/PLANCKEV);

  // plasmon frequencies (eV conversion)
  d2l_lorentz_deltap[1][0] = (0.82*(2.0*PI/PLANCKEV))/2.0;
  d2l_lorentz_deltap[1][1] = (0.82*(2.0*PI/PLANCKEV))/2.0;



  //---------------------------------------------------------
  // ITO - SOPRA DATABASE
  //---------------------------------------------------------
 
  // infinite frequency epsilon
  d2l_inf[1] = 1.58198704500344;

  // conductivity
  d2l_conduct[1] = 0.0;

  // DRUDE MODEL

  // plasmon frequency (eV conversion)
  d2l_drude_omegap[1] = 1.61312991416747*(2.0*PI/PLANCKEV);
  // damping constant (eV conversion)
  d2l_drude_gammap[1] = 0.442862763442952*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL

  // assign epsilon widths
  d2l_lorentz_depsr[1][0] = 2.37383602797644*0.891550305353056;
  d2l_lorentz_depsr[1][1] = 2.37383602797644*(1.0-0.891550305353056);

  // plasmon frequencies (eV conversions)
  d2l_lorentz_omegap[1][0] = 6.87375733506245*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[1][1] = 4.78008554084326*(2.0*PI/PLANCKEV);

  // plasmon frequencies (eV conversion)
  d2l_lorentz_deltap[1][0] = (1.439454263791935E-004)*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[1][1] = 0.656146754926608*(2.0*PI/PLANCKEV);


  //=========================================================
  // SILICON
  //=========================================================

*/
  //=========================================================
  // CARBON
  //=========================================================

/*
  //---------------------------------------------------------
  // AMORPHEOUS CARBON - SOPRA DATABASE (300-840nm)
  //---------------------------------------------------------

  // METAL INFINITE FREQUENCY
  epsr_metal_inf[3] = 3.09622357327939;

  // CONDUCTIVITY
  metal_conductivity[3] = 0.0;

  // DRUDE MODEL
  // plasmon frequency (eV conversion)
  drude_omegap[3] = 1.21844814747305*(2.0*PI/PLANCKEV);
  // damping constant (eV conversion)
  drude_gammap[3] = (2.401419597583212E-005)*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  lorentz_depsr[3][0] = 0.759190185440699*(0.356086064908049);
  lorentz_depsr[3][1] = 0.759190185440699*(1.0-0.356086064908049);
  // plasmon frequencies (eV conversions)
  lorentz_omegap[3][0] = 3.52940394700369*(2.0*PI/PLANCKEV);
  lorentz_omegap[3][1] = 4.42518950947219*(2.0*PI/PLANCKEV);
  // plasmon frequencies (eV conversion)
  lorentz_deltap[3][0] = 0.498775653651519*(2.0*PI/PLANCKEV);
  lorentz_deltap[3][1] = 0.419181895479968*(2.0*PI/PLANCKEV);
*/

/*
  //---------------------------------------------------------
  // AMORPHEOUS CARBON
  //---------------------------------------------------------
  // i. this was taken from: PRB 31, 8097 (1985)
  // ii. this does not take into account any hydrogenation

  // METAL INFINITE FREQUENCY
  d2l_inf[2] = 2.554;

  // CONDUCTIVITY
  d2l_conduct[2] = 0.0;

  // DRUDE MODEL
  // plasmon frequency
  d2l_drude_omegap[2] = 1.415*(2.0*PI/PLANCKEV);
  // damping constant
  d2l_drude_gammap[2] = 6.237*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  d2l_lorentz_depsr[2][0] = 6.078*(0.587);
  d2l_lorentz_depsr[2][1] = 6.078*(1.0-0.587);
  // plasmon frequencies
  d2l_lorentz_omegap[2][0] = 4.489*(2.0*PI/PLANCKEV);
  d2l_lorentz_omegap[2][1] = 1.865*(2.0*PI/PLANCKEV);
  // plasmon frequencies
  d2l_lorentz_deltap[2][0] = 3.159*(2.0*PI/PLANCKEV);
  d2l_lorentz_deltap[2][1] = 1.202*(2.0*PI/PLANCKEV);
*/

/*
  //=========================================================
  // ALUMINUM
  //=========================================================

  //---------------------------------------------------------
  // Al (Pure) - SOPRA Database (300-1000nm)
  //---------------------------------------------------------

  // METAL INFINITE FREQUENCY
  epsr_metal_inf[4] = 1.81169007835483;

  // CONDUCTIVITY
  metal_conductivity[4] = 0.0;

  // DRUDE MODEL
  // plasmon frequency (eV conversion)
  drude_omegap[4] = 13.4014051138924*(2.0*PI/PLANCKEV);
  // damping constant (eV conversion)
  drude_gammap[4] = (6.920034702273238E-002)*(2.0*PI/PLANCKEV);

  // LORENTZ MODEL
  // assign epsilon widths
  lorentz_depsr[4][0] = 19.3928024262706*(0.440463385222862);
  lorentz_depsr[4][1] = 19.3928024262706*(1.0-0.440463385222862);
  // plasmon frequencies (eV conversions)
  lorentz_omegap[4][0] = 1.56781717564189*(2.0*PI/PLANCKEV);
  lorentz_omegap[4][1] = 1.96622898362407*(2.0*PI/PLANCKEV);
  // plasmon frequencies (eV conversion)
  lorentz_deltap[4][0] = 0.210537586110392*(2.0*PI/PLANCKEV);
  lorentz_deltap[4][1] = 0.881209383616846*(2.0*PI/PLANCKEV);
*/



  //=========================================================
  // NON-LOCAL HYDRODYNAMIC DRUDE MODEL
  //=========================================================

  // NON-LOCAL SILVER FOR A (2.5 nm, 3.0 nm) CORE-SHELL PARTICLE
  // i. these parameters were taken from Leung, PRB 73, 125438 (2006)
  // !!! I should probably put this in external files
  nl_dh_vf[0] = 1.39e6;
  nl_dh_inf[0] = 1.0;
  nl_dh_omegap[0] = 1.36e16;
  nl_dh_gammap[0] = 2.56e13 + 1.0*nl_dh_vf[0]/(3.0*NM - 2.5*NM);


  // LIBRARY HAS BEEN GENERATED
  library_generated = 1;

  return;
}




//========================================================================
//========================================================================
//
//	NAME:	complex<double> get_epsr(int imaterial, double wavelength)
//	DESC:	User-defined parameters
//
//	NOTES: 	i. the sign of the imaginary part of epsilon should be opposite to
// 		the sign in front of the time exponent
//
//	PARAMS: [in] imaterial == "universal" material identifier
//              [in] wavelength == wavelength (for wavelength-dependent dielectrics)
//              [in] exp_convention == +/-1 in front of time phasor factor exp(+/-i*omega*t)
//
//		[out] epsr == relative dielectric
//
//========================================================================
//========================================================================
complex<double> get_epsr(int imaterial, double wavelength, int exp_convention)
{

  // FIRST CHECK FOR FREE SPACE
  if(imaterial >= 90)
  {
    if(imaterial == 99)
    {
      return COMPLEXONE;
    }  
    else
    {
      return epsr_custom[imaterial-90];
    }
  }

  complex<double> epsr;

  int l;

  // GENERATE THE LIBRARY IF USER HAS NOT DONE SO
  if(library_generated == 0)
  {
    generate_epsrlib();
  }

  // GET THE WAVENUMBER
  double wavenumber = 2.0*PI*CSPEED/wavelength;

  // SET EPSR TO EPSILON INFINITY
  epsr = d2l_inf[imaterial];

  // ADD IN DRUDE TERMS
  epsr -= d2l_drude_omegap[imaterial]*d2l_drude_omegap[imaterial]/(wavenumber*wavenumber - COMPLEXJ*wavenumber*d2l_drude_gammap[imaterial]);

  // ADD IN LORENTZ TERMS
  for(l = 0; l < 2; ++l)
  {
    epsr += d2l_lorentz_depsr[imaterial][l]*d2l_lorentz_omegap[imaterial][l]*d2l_lorentz_omegap[imaterial][l]/(d2l_lorentz_omegap[imaterial][l]*d2l_lorentz_omegap[imaterial][l] + 2.0*COMPLEXJ*wavenumber*d2l_lorentz_deltap[imaterial][l] - wavenumber*wavenumber);
  }

  // RETURN EPSILON WITH THE APPROPRIATE SIGN
  // i. the sign of the imaginary part of epsilon should be opposite to
  // the sign in front of the time exponent
  if(exp_convention == 1)
  {
    return epsr;
  }
  else if(exp_convention == -1)
  {
    return complex<double>(real(epsr), fabs(imag(epsr)));
  }
  else
  {
    cout << "error in get_epsr(*), exp_convention " << exp_convention << " passed." << endl;
    exit(1);
  }
}

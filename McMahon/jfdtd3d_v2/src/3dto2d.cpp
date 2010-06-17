/*
  Copyright (c) 2006 - 2008  Jeffrey M. McMahon

  This file is part of JFDTD3D.

  JFDTD3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  JFDTD3D is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with JFDTD3D.  If not, see <http://www.gnu.org/licenses/>.
*/


//************************************************************************
//------------------------------------------------------------------------
//
// NAME: 3dto2d: CONVERTS A 3D PLOT TO A 2D PLOT
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: 3dto2d.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// DESC:	Converts a file with [x  y  z] data to [x  z] or [y  z] data.
//
// NOTES:	i. last updated on 7/04/08
//		ii. the file to read in needs to be in the form:
//			[x  y  z]
//		iv. the parameters need to be set at the top of the 
//		read_parameters() subroutine.
//
// TODO:	i. 
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************


using namespace std;	


//************************************************************************
//			INCLUDE FILES
//************************************************************************

// STANDARD
// utilities
#include <cstdlib>
// standard math
#include <cmath>
// screen I/O
#include <iostream>
// file I/O
#include <fstream> 

// CUSTOM

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			SUBROUTINES
//************************************************************************

void read_parameters();	


//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

int itype;

ifstream filein;
ofstream fileout;

double position, tolerance;


//========================================================================
//========================================================================
//
// 	NAME:	main() 
//	DESC:	Entryway to program
//
//	NOTES:
//
//========================================================================
//========================================================================
int main(int argc, char** argv)
{

  // local indices
  int i, l;


  //=========================================================
  // INITIALIZATION
  //=========================================================

  read_parameters();	

  //=========================================================
  // NOW READ IN
  //=========================================================

  double xpos, ypos, zpos;

  // while there is another xpos to read in ...
  while(filein >> xpos) 
  {
//cout << "here." << endl;
    // read in the ypos & zpos
    filein >> ypos;
    filein >> zpos;

    // ------------------------ xz ------------------------  
    if(itype==1)
    {
      if( (abs(ypos-position) <= tolerance) )
      {
        fileout << xpos << "     " << zpos << endl;
      }
    }
    // ------------------------ yz ------------------------  
    else if(itype==2)
    {
      if( (abs(xpos-position) <= tolerance) )
      {
        fileout << ypos << "     " << zpos << endl;
      }
    }

  }



  //=========================================================
  // CLEANUP
  //=========================================================
  filein.close();
  fileout.close();

  return 0;
}


//========================================================================
//========================================================================
//
//	NAME:	read_parameters()
//	DESC:	Read in initial parameters
//
//	NOTES: 	i.
//
//
//========================================================================
//========================================================================
void read_parameters()
{

  // FILENAME
  filein.open("e_field_stitched.3.10", ios::in);
  fileout.open("circ_hh.dat", ios::app);


  // TYPE OF PLOT
  // i. types:
  //	1 == xz
  //	2 == yz
  itype=1;

  // POSITION
  // i. the desired constant position needs to be set, as well
  // as a tolerance that dictates how close we need to be to consider
  // the correct position
  // ii. the tolerance better be set smaller than the distance between constant
  // values of two or more slices may be output
  position=1600.0e-9;
  tolerance=0.5e-9;

  return;

}



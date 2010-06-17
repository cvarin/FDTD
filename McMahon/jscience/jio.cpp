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
// NAME: JSCIENCE INPUT/OUTPUT: STANDARD SUBROUTINES FOR INPUT/OUTPUT
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jio.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 1/3/09
//
// TODO:
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************


//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "jio.h"

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************

//************************************************************************
//			SUBROUTINES
//************************************************************************

//========================================================================
//========================================================================
//
//	NAME:	void get_comments(ifstream &fin)
//	DESC:	Reads in comments from an input file.  The input files should
//		be set up the following way (! the # and % are important):
//
//		# comment
//		% variable_name =
//		value
//
//	NOTES: 	i.
//
//
//========================================================================
//========================================================================
void get_comments(ifstream &fin)
{
  string temp;
  fin >> temp;
  while(temp[0] == '#')
  {
    getline(fin, temp);
    fin >> temp;
  }
  getline(fin, temp);

  return;
}



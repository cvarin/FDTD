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
// NAME: JSCIENCE MATRIX CLASSES
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jmath_matrix.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 9/27/08
//
//	NAME: 	
//
//	DESC: 	
//
//	NOTES:	1) 
//
//	TODO:
//		1) 
//
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

#ifndef JMATRIX_H
#define JMATRIX_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************


#include <iostream>
#include <list>
#include <vector>

#include <cmath>
#include <complex>

using namespace std;	


//=========================================================
// SPARSE MATRIX (NON-ORDERED) STORAGE
//=========================================================

//---------------------------------------------------------
// DOUBLE
//---------------------------------------------------------
class VMATRIXSP_D
{
  public:
    vector< double > val;
    vector< int > col, row;    

    VMATRIXSP_D(int r, int c);
    ~VMATRIXSP_D();

    int nrows, ncols, nzeros;
};

//---------------------------------------------------------
// COMPLEX DOUBLE
//---------------------------------------------------------
class VMATRIXSP_CD
{
  public:
    vector< complex<double> > val;
    vector< int > col, row;    

    VMATRIXSP_CD(int r, int c);
    ~VMATRIXSP_CD();

    void reset();

    int nrows, ncols, nzeros;
};

//=========================================================
// COMPRESSED STORAGE (BY) ROW FORMAT
//=========================================================
// i. the range of the rowpntr and colpntr are extra big because they
// are used as ranges in loops and we need the last position [nrows+1], etc.

//---------------------------------------------------------
// COMPLEX DOUBLE
//---------------------------------------------------------

class MATRIXCSR_CD
{
  public:
    complex<double> *val;
    int *col, *rowpntr;

    MATRIXCSR_CD(int dim, int nel);
    ~MATRIXCSR_CD();

    int nrows, nzeros;
};


class VMATRIXCSR_CD
{
  public:
    vector< complex<double> > val;
    vector< int > col;    
    int *rowpntr;

    VMATRIXCSR_CD(int n);
    ~VMATRIXCSR_CD();

    void insert(complex<double> v, int r, int c);
    void remove(int k);

    int nrows, nzeros;
};



//=========================================================
// COMPRESSED STORAGE (BY) COLUMN FORMAT
//=========================================================

//---------------------------------------------------------
// ARRAY FORMAT
//---------------------------------------------------------

class MATRIXCSC_CD
{
  public:
    complex<double> *val;
    int *row, *colpntr;

    MATRIXCSC_CD(int dim, int nel);
    ~MATRIXCSC_CD();

    void insert(complex<double> v, int r, int c);
    void remove(int k);

    int ncols, nzeros;

  private:
    int nsize;
};


//---------------------------------------------------------
// VECTOR FORMAT
//---------------------------------------------------------

class VMATRIXCSC_CD
{
  public:
    vector< complex<double> > val;
    vector< int > row;    
    int *colpntr;

    VMATRIXCSC_CD(int n);
    ~VMATRIXCSC_CD();

    void insert(complex<double> v, int r, int c);
    void remove(int k);

    int ncols, nzeros;
};


//---------------------------------------------------------
// LIST FORMAT
//---------------------------------------------------------

class LMATRIXCSC_CD
{
  public:
    list< complex<double> > val;
    list< int > row;    
    int *colpntr;

    LMATRIXCSC_CD(int n);
    ~LMATRIXCSC_CD();

    //void get_column(int j, vector< complex<double> > &vect);
    void get_column(int j, int &nvals, vector< complex<double> > &vect, vector< int > &rows);
    void insert_column(int j, int nvals, vector< complex<double> > vect, vector< int > rows, double tol);

    int ncols, nzeros;
};



//************************************************************************
//			SUBROUTINES
//************************************************************************



void convert_matrix(VMATRIXSP_CD* smatrix, MATRIXCSR_CD*& csrmat);

//void multiply(MATRIXCSC_CD*& Amat, MATRIXCSC_CD*& Bmat, MATRIXCSR_CD*& Cmat);
//void get_transpose(MATRIXCSC_CD*& cscmat, MATRIXCSC_CD*& cscmatT);

//void convert_matrix(MATRIX_SPARSE_CD smatrix, MATRIXCSR_CD*& csrmat);
//void convert_matrix(MATRIX_SPARSE_CD smatrix, MATRIXCSC_CD*& cscmat);

//void convert_matrix(MATRIXCSC_CD*& origmat, MATRIXCSR_CD*& convertmat);
//void convert_matrix(MATRIXCSR_CD*& origmat, MATRIXCSC_CD*& convertmat);

//void convert_matrix(MATRIXCSC_CD*& origmat, VMATRIXSP_CD*& convertmat);


#endif


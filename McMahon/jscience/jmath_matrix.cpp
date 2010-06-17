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
// FILE: jmath_matrix.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 8/25/08
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

//************************************************************************
//			INCLUDE FILES
//************************************************************************	
	
#include "jmath_matrix.h"

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************


//************************************************************************
//			SUBROUTINES
//************************************************************************


//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

//=========================================================
// SPARSE MATRIX (NON-ORDERED) STORAGE
//=========================================================

//---------------------------------------------------------
// DOUBLE
//---------------------------------------------------------

VMATRIXSP_D::VMATRIXSP_D(int r, int c)
{
  nrows=r;
  ncols=c;
  nzeros=0;

  val.clear();
  row.clear();
  col.clear();
}

VMATRIXSP_D::~VMATRIXSP_D()
{

}


//---------------------------------------------------------
// COMPLEX DOUBLE
//---------------------------------------------------------

VMATRIXSP_CD::VMATRIXSP_CD(int r, int c)
{
  nrows=r;
  ncols=c;

  reset();
}

VMATRIXSP_CD::~VMATRIXSP_CD()
{

}


void VMATRIXSP_CD::reset()
{
  nzeros=0;

  val.clear();
  row.clear();
  col.clear();
}


//=========================================================
// COMPRESSED STORAGE (BY) ROW FORMAT
//=========================================================
// i. the range of the rowpntr and colpntr are extra big because they
// are used as ranges in loops and we need the last position [nrows+1], etc.



//---------------------------------------------------------
// ARRAY FORMAT
//---------------------------------------------------------

MATRIXCSR_CD::MATRIXCSR_CD(int dim, int nel)
{
  nrows=dim;
  nzeros=nel;

  val=new complex<double> [nzeros+1];
  col=new int [nzeros+1];
  rowpntr=new int [nrows+2];
}


MATRIXCSR_CD::~MATRIXCSR_CD()
{
  delete [] val;
  delete [] col;
  delete [] rowpntr;
}


//---------------------------------------------------------
// VECTOR FORMAT
//---------------------------------------------------------

VMATRIXCSR_CD::VMATRIXCSR_CD(int n)
{
  int i;

  nzeros=0;
  nrows=n;

  rowpntr=new int [nrows+2];

  for(i=1; i<=nrows+1; ++i)
  {
    rowpntr[i]=1;
  }
}

// Insert a value
void VMATRIXCSR_CD::insert(complex<double> v, int r, int c)
{

  int i, k;
 
  // LOOP OVER THE STORED VALUES FOR THIS COLUMN
  for(i=rowpntr[r]; i<rowpntr[r+1]; ++i)
  {
    // IF WE HAVE THE EXACT POSITION STORED
    if(col.at(i)==c)
    {
      val.at(i)+=v;
      return;
    }
    // ELSE THE CURRENT ROW IS HIGHER THAN THAT ADDING, SO ADD IN THE VALUE AT THIS i POSITION 
    else if(col[i]>c)
    {
      nzeros++;

      val.insert(val.begin() + i, v);
      col.insert(col.begin() + i, c);

      for(k=r+1;k<=nrows+1;++k)
      {
        rowpntr[k]++;
      }

      return;
      
    }
    // ELSE WE NEED TO ADD A VALUE AT THE END OF THE COLUMN
    else if(i==(rowpntr[r+1]-1))
    {

      nzeros++;

      val.insert(val.begin() + i, v);
      col.insert(col.begin() + i, c);

      for(k=r+1;k<=nrows+1;++k)
      {
        rowpntr[k]++;
      }

      return;
    }
  }

// j11
  // IF WE MAKE IN HERE THEN NO VALUES FOR THIS COLUMN HAVE BEEN STORED
  // SO ADD A NEW ROW
  nzeros++;

  for(k=r+1;k<=nrows+1;++k)
  {
    rowpntr[k]++;
  }

  val.insert(val.begin() + rowpntr[r], v);
  col.insert(col.begin() + rowpntr[r], c);

  //cout << "I don't think we should be here. " << endl;
}

// Remove a value
void VMATRIXCSR_CD::remove(int k)
{

  int i, j;

  for(i=1; i<=nrows; ++i)
  {
    if( (k>=rowpntr[i]) && (k < rowpntr[i+1]))
    {
      for(j=i+1; j<=nrows-1; ++j)
      {
        rowpntr[j]--;
      }
      break;
    }
  }

  for(i=k; i<nzeros; ++i)
  {
    val[i]=val[i+1];
    col[i]=col[i+1];
  }

  nzeros--;

}


VMATRIXCSR_CD::~VMATRIXCSR_CD()
{
  delete [] rowpntr;
}



//=========================================================
// COMPRESSED STORAGE (BY) COLUMN FORMAT
//=========================================================

//---------------------------------------------------------
// ARRAY FORMAT
//---------------------------------------------------------

MATRIXCSC_CD::MATRIXCSC_CD(int dim, int nel)
{
  nzeros=0;
  ncols=dim;
  nsize=nel;

  val=new complex<double> [nsize+1];
  row=new int [nsize+1];
  colpntr=new int [ncols+2];
}

// Insert a value
void MATRIXCSC_CD::insert(complex<double> v, int r, int c)
{

  int i, k;
 
  // LOOP OVER THE STORED VALUES FOR THIS COLUMN
  for(i=colpntr[c]; i<colpntr[c+1]; ++i)
  {
    // IF WE HAVE THE EXACT POSITION STORED
    if(row[i]==r)
    {
      val[i]+=v;
      return;
    }
    // ELSE THE CURRENT ROW IS HIGHER THAN THAT ADDING, SO ADD IN THE VALUE AT THIS i POSITION 
    else if(row[i]>r)
    {
      nzeros++;
      if(nzeros>nsize)
      {
        cout << "not enough space allocated in MATRIXCSC" << endl;
      }

      for(k=nzeros;k>=i+1;--k)
      {
        val[k]=val[k-1];
        row[k]=row[k-1];
      }

      val[i]=v;
      row[i]=r;

      for(k=c+1;k<=ncols+1;++k)
      {
        colpntr[k]++;
      }

      return;
      
    }
    // ELSE WE NEED TO ADD A VALUE AT THE END OF THE COLUMN
    else if(i==(colpntr[c+1]-1))
    {

      nzeros++;
      if(nzeros>nsize)
      {
        cout << "not enough space allocated in MATRIXCSC" << endl;
      }

      for(k=nzeros;k>=i+1;--k)
      {
        val[k]=val[k-1];
        row[k]=row[k-1];
      }

      val[i]=v;
      row[i]=r;

      for(k=c+1;k<=ncols+1;++k)
      {
        colpntr[k]++;
      }

      return;
    }
  }

  // IF WE MAKE IN HERE THEN NO VALUES FOR THIS COLUMN HAVE BEEN STORED
  // SO ADD A NEW COLUMN
  cout << "I don't think we should be here. " << endl;
}

// Remove a value
void MATRIXCSC_CD::remove(int k)
{

  int i, j;

  for(i=1; i<=ncols; ++i)
  {
    if( (k>=colpntr[i]) && (k < colpntr[i+1]))
    {
      for(j=i+1; j<=ncols-1; ++j)
      {
        colpntr[j]--;
      }
      break;
    }
  }

  for(i=k; i<nzeros; ++i)
  {
    val[i]=val[i+1];
    row[i]=row[i+1];
  }

  nzeros--;

}


MATRIXCSC_CD::~MATRIXCSC_CD()
{
  delete [] val;
  delete [] row;
  delete [] colpntr;
}




//---------------------------------------------------------
// VECTOR FORMAT
//---------------------------------------------------------


VMATRIXCSC_CD::VMATRIXCSC_CD(int n)
{
  int i;

  nzeros=0;
  ncols=n;

  colpntr=new int [ncols+2];

  for(i=1; i<=ncols+1; ++i)
  {
    colpntr[i]=1;
  }
}

// Insert a value
void VMATRIXCSC_CD::insert(complex<double> v, int r, int c)
{

  int i, k;
 
  // LOOP OVER THE STORED VALUES FOR THIS COLUMN
  for(i=colpntr[c]; i<colpntr[c+1]; ++i)
  {
    // IF WE HAVE THE EXACT POSITION STORED
    if(row.at(i)==r)
    {
      val.at(i)+=v;
      return;
    }
    // ELSE THE CURRENT ROW IS HIGHER THAN THAT ADDING, SO ADD IN THE VALUE AT THIS i POSITION 
    else if(row[i]>r)
    {
      nzeros++;

      val.insert(val.begin() + i, v);
      row.insert(row.begin() + i, r);

      for(k=c+1;k<=ncols+1;++k)
      {
        colpntr[k]++;
      }

      return;
      
    }
    // ELSE WE NEED TO ADD A VALUE AT THE END OF THE COLUMN
    else if(i==(colpntr[c+1]-1))
    {

      nzeros++;

      val.insert(val.begin() + i, v);
      row.insert(row.begin() + i, r);

      for(k=c+1;k<=ncols+1;++k)
      {
        colpntr[k]++;
      }

      return;
    }
  }

  // IF WE MAKE IN HERE THEN NO VALUES FOR THIS COLUMN HAVE BEEN STORED
  // SO ADD A NEW COLUMN
  nzeros++;

  for(k=c+1;k<=ncols+1;++k)
  {
    colpntr[k]++;
  }

  val.insert(val.begin() + colpntr[c], v);
  row.insert(row.begin() + colpntr[c], r);

  //cout << "I don't think we should be here. " << endl;
}

// Remove a value
void VMATRIXCSC_CD::remove(int k)
{

  int i, j;

  for(i=1; i<=ncols; ++i)
  {
    if( (k>=colpntr[i]) && (k < colpntr[i+1]))
    {
      for(j=i+1; j<=ncols-1; ++j)
      {
        colpntr[j]--;
      }
      break;
    }
  }

  for(i=k; i<nzeros; ++i)
  {
    val[i]=val[i+1];
    row[i]=row[i+1];
  }

  nzeros--;

}


VMATRIXCSC_CD::~VMATRIXCSC_CD()
{
  delete [] colpntr;
}


//---------------------------------------------------------
// LIST FORMAT
//---------------------------------------------------------

LMATRIXCSC_CD::LMATRIXCSC_CD(int n)
{
  int i;

  nzeros=0;
  val.clear();
  row.clear();
  ncols=n;

  colpntr=new int [ncols+2];

  for(i=1; i<=ncols+1; ++i)
  {
    colpntr[i]=1;
  }
}

/*
// Get a column
void LMATRIXCSC_CD::get_column(int j, vector< complex<double> > &vect)
{
  int i;

  vect.clear();
  vect.resize(ncols+1,0.0);

  list< complex<double> >::iterator itv=val.begin();
  list< int >::iterator itr=row.begin();

  advance(itv,colpntr[j]);
  advance(itr,colpntr[j]);

  for(i=1;i<=colpntr[j+1]-colpntr[j];++i)
  {
    vect.at(*itr)=*itv;
    itv++;
    itr++;
  }

}
*/


void LMATRIXCSC_CD::get_column(int j, int &nvals, vector< complex<double> > &vect, vector< int > &rows)
{
  int i;

  nvals = colpntr[j+1]-colpntr[j];

  vect.clear();
  rows.clear();

  list< complex<double> >::iterator itv=val.begin();
  list< int >::iterator itr=row.begin();

  advance(itv,colpntr[j]);
  advance(itr,colpntr[j]);

  for(i=0;i<nvals;++i)
  {
    vect.push_back(*itv);
    rows.push_back(*itr);
    itv++;
    itr++;
  }

}

/*
// Remove a value
void LMATRIXCSC_CD::insert_column(int j, vector< complex<double> > vect)
{

  double tol=1.0e-2;

  int i;

  //vector< complex<double> > colvect;

  //get_column(j, colvect);

  //for(i=1;i<=ncols;++i)
  //{ 
  // vect.at(i)+=colvect.at(i);
  //}

  list< complex<double> >::iterator itv1, itv2;
  list< int >::iterator itr1, itr2;
  itv1 = itv2 = val.begin();
  itr1 = itr2 = row.begin();
  advance(itv1,colpntr[j]);
  advance(itr1,colpntr[j]);
  advance(itv2,colpntr[j+1]);
  advance(itr2,colpntr[j+1]);

  // i. this will delete all elements between it1 and it2 including
  // it1 but not it2, and will return an iterator pointing to it2
  itv2 = val.erase(itv1, itv2);
  itr2 = row.erase(itr1, itr2);

  // i. insert the elements right before it2
  int nnz=0;
  for(i=1;i<=ncols;++i)
  { 
    if(abs(vect.at(i))>tol)
    {
      nnz++;
      val.insert(itv2, vect.at(i));
      row.insert(itr2, i);      
    }
  }


  int diff=nnz-(colpntr[j+1]-colpntr[j]);

  for(i=j+1; i<=ncols+1; ++i)
  {
    colpntr[i]+=diff;
  }


  nzeros+=nnz;


}
*/

// Remove a value
void LMATRIXCSC_CD::insert_column(int j, int nvals, vector< complex<double> > vect, vector< int > rows, double tol)
{
  int k;

  list< complex<double> >::iterator itv=val.begin();
  list< int >::iterator itr=row.begin();
/*
  advance(itv,colpntr[j+1]);
  advance(itr,colpntr[j+1]);

  //int nnz;
  for(k=0;k<nvals;++k)
  {
    //if(abs(vect.at(k))>0)
    //{
      //nnz++;
      val.insert(itv, vect.at(k));
      row.insert(itr, rows.at(k));      
    //}
  }
*/

  advance(itv,colpntr[j]);
  advance(itr,colpntr[j]);

  int nnz=0;
  int n=colpntr[j];

  for(k=0;k<nvals;++k)
  {
    if(n<colpntr[j+1])
    {
      // if we found an identical row ...
      if(*itr==rows.at(k))
      {
        *itv+=vect.at(k);
      }
      // else if the first row is higher than the one to add ...
      else if(*itr>rows.at(k))
      {
        if(abs(vect.at(k))>tol)
        {
          nnz++;
          val.insert(itv, vect.at(k));
          row.insert(itr, rows.at(k));      
        }
      }
      else
      {
        itr++;
        itv++;
        n++;
        k--;
      }
    }
    else
    {
      if(*itr==rows.at(k))
      {
        *itv+=vect.at(k);
      }
      else
      {
        if(abs(vect.at(k))>tol)
        {
          nnz++;
          val.insert(itv, vect.at(k));
          row.insert(itr, rows.at(k));      
        }
      }
    }

  } //++k

  for(k=j+1; k<=ncols+1; ++k)
  {
    colpntr[k]+=nnz;
  }

  nzeros+=nnz;

}


LMATRIXCSC_CD::~LMATRIXCSC_CD()
{
  delete [] colpntr;
}




/*
void get_transpose(MATRIXCSC_CD*& cscmat, MATRIXCSC_CD*& cscmatT)
{

  int k;

  MATRIXCSR_CD* csrmatT;
  csrmatT = new MATRIXCSR_CD(ninterpol_pts_assigned, cscmat->nzeros);
  csrmatT->nzeros=cscmat->nzeros;

  for(k=1;k<=cscmat->nzeros;++k)
  {
    csrmatT->val[k]=cscmat->val[k];
    csrmatT->col[k]=cscmat->row[k];
  }

  for(k=1;k<=ninterpol_pts_assigned+1;++k)
  {
    csrmatT->rowpntr[k]=cscmat->colpntr[k];
  }

  convert_matrix(csrmatT, cscmatT);

  delete csrmatT;

  return;
}

// j11
// i. remember there is no memory for the Cmat declares
void multiply(MATRIXCSC_CD*& Amat, MATRIXCSC_CD*& Bmat, MATRIXCSR_CD*& Cmat)
{

  int i, j, k;

  complex<double> *rowvals, COMPLEXZERO(0.0,0.0), sum;
  rowvals = new complex<double> [ninterpol_pts_assigned+1];

  // CONVER THE [A] MATRIX TO CSR FORMAT
  MATRIXCSR_CD* Amatr;
  convert_matrix(Amat, Amatr);

  // !!! make the Ctemp matrix to be totally full 
  //MATRIXCSC_CD* Ctemp;
  //Ctemp = new MATRIXCSC_CD(ninterpol_pts_assigned, ninterpol_pts_assigned*ninterpol_pts_assigned);
  //Cmat = new MATRIXCSR_CD(ninterpol_pts_assigned, ninterpol_pts_assigned*ninterpol_pts_assigned);
  //Cmat->nzeros=ninterpol_pts_assigned*ninterpol_pts_assigned;

  VMATRIXCSR_CD* vCmat;
  vCmat = new VMATRIXCSR_CD(ninterpol_pts_assigned);
  vCmat->val.push_back(COMPLEXZERO);
  vCmat->col.push_back(0);

  for(i=1;i<=ninterpol_pts_assigned;++i)
  {
    // GET rowvals[] 
    for(k=Amatr->rowpntr[i];k<Amatr->rowpntr[i+1];++k)
    { 
      rowvals[Amatr->col[k]]=Amatr->val[k];
    }

    for(j=1;j<=ninterpol_pts_assigned;++j)
    {

      // MULTIPLY
      sum=COMPLEXZERO;
      for(k=Bmat->colpntr[j];k<Bmat->colpntr[j+1];++k)
      { 
        sum+=rowvals[Bmat->row[k]]*Bmat->val[k];
      }

      if(sum!=COMPLEXZERO)
      {
        vCmat->insert(sum, i, j);
      }
    }

    // RESET rowvals[] 
    for(k=Amatr->rowpntr[i];k<Amatr->rowpntr[i+1];++k)
    { 
      rowvals[Amatr->col[k]]=COMPLEXZERO;
    }

  } // ++i

  // useful for loop ranges
  vCmat->rowpntr[ninterpol_pts_assigned+1]=vCmat->nzeros+1;


  Cmat = new MATRIXCSR_CD(ninterpol_pts_assigned, vCmat->nzeros);
  Cmat->nzeros=vCmat->nzeros;
  for(k=1;k<=vCmat->nzeros;++k)
  {
    Cmat->val[k]=vCmat->val.at(k);
    Cmat->col[k]=vCmat->col.at(k);
  }
  for(k=1;k<=ninterpol_pts_assigned+1;++k)
  {
    Cmat->rowpntr[k]=vCmat->rowpntr[k];
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] rowvals;

  delete Amatr;
  //delete Ctemp;
  delete vCmat;

  return;
}
*/

void convert_matrix(VMATRIXSP_CD* smatrix, MATRIXCSR_CD*& csrmat)
{

  int k,kk;

  complex<double> COMPLEXZERO(0.0,0.0);

  // SET UP A TEMP ARRAY TO HOLD ROW VALUES
  complex<double> *rowvals;
  rowvals=new complex<double> [smatrix->nrows+1];  

  // SET UP A TEMP CSR MATRIX TO MAXIMUM POSSIBLE SIZE
  // i. i.e. this size will be correct if there is no duplicate position values in smatrix
  MATRIXCSR_CD *csrmat_temp;
  csrmat_temp=new MATRIXCSR_CD(smatrix->nrows, smatrix->nzeros);
  csrmat_temp->nzeros=smatrix->nzeros;

  int nnz=0;


  // LOOP OVER ROWS
  for(k=1;k<=smatrix->nrows;++k)
  {
    // RESET ARRAY
    for(kk=1;kk<=smatrix->nrows;++kk)
    {
      rowvals[kk]=COMPLEXZERO;
    }

    // FILL ROWVALS ARRAY
    for(kk=0;kk<smatrix->nzeros;++kk)
    {
      if(smatrix->row.at(kk)==k)
      {
        rowvals[smatrix->col.at(kk)] += smatrix->val.at(kk);
      }
    }

    // SET ROWPNTR TO POINT TO 1 HIGHER THAN CURRENT NONZERO ELEMENTS
    // i. ! each row bettwe contain at least 1 nonzero element
    csrmat_temp->rowpntr[k]=nnz+1;

    // loop over "columns"
    for(kk=1;kk<=smatrix->ncols;++kk)
    {
      if(rowvals[kk]!=COMPLEXZERO)
      {
        // INCREMENT NONZERO ELEMENTS
        nnz++;

        // ADD VALUE
        csrmat_temp->val[nnz]=rowvals[kk];
        csrmat_temp->col[nnz]=kk;
      }
    }
    
  } // ++k

  // this is used in loops
  csrmat_temp->rowpntr[smatrix->nrows+1]=nnz+1;


  //=========================================================
  // COPY THE TEMP CSR MATRIX TO THE REAL CSR MATRIX
  //=========================================================
  csrmat=new MATRIXCSR_CD(smatrix->nrows, nnz);
  csrmat->nzeros=nnz;

  for(k=1;k<=nnz;++k)
  {
    csrmat->val[k]=csrmat_temp->val[k];
    csrmat->col[k]=csrmat_temp->col[k];
  }

  for(k=1;k<=smatrix->nrows+1;++k)
  {
    csrmat->rowpntr[k]=csrmat_temp->rowpntr[k];
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] rowvals;

  delete csrmat_temp;


  return;
}


/*


// i. this assumes a square matrix
void convert_matrix(MATRIXCSC_CD*& origmat, VMATRIXSP_CD*& convertmat)
{

  int j, k;

  complex<double> COMPLEXZERO(0.0,0.0);

  convertmat = new VMATRIXSP_CD(origmat->ncols, origmat->ncols);
  convertmat->nzeros = origmat->nzeros;

  // LOOP OVER COLUMNS
  for(j=1;j<=ninterpol_pts_assigned;++j)
  {
    for(k=origmat->colpntr[j];k<origmat->colpntr[j+1];++k)
    {
      convertmat->val.push_back(origmat->val[k]);
      convertmat->row.push_back(origmat->row[k]);
      convertmat->col.push_back(j);
    }
  }  
  
  return;
}

void convert_matrix(MATRIX_SPARSE_CD smatrix, MATRIX_CSR_CD csrmat)
{

  int k,kk;

  complex<double> *val;
  val=new complex<double> [smatrix.size+1];
  complex<double> *rowvals;
  rowvals=new complex<double> [ninterpol_pts_assigned+1];
  int *col,*rowpntr;
  col=new int [smatrix.size+1];
  rowpntr=new int [ninterpol_pts_assigned+1];

  int nnz=0;

  cout << "here." << endl;


  // loop over rows
  for(k=1;k<=ninterpol_pts_assigned;++k)
  {

    // reset array
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
      rowvals[kk]=complex_zero;

    // set the rowpointer to be 1 higher than the number of nonzero elements
    // i. here we need to be a little careful that each row in the matrix needs
    // at least 1 element
    rowpntr[k]=nnz+1;

    // fill in array with matrix
    for(kk=1;kk<=smatrix.size;++kk)
      if(smatrix.ia[kk]==k)
        rowvals[smatrix.ja[kk]]+=smatrix.a[kk];

    // loop over "columns"
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
      if(rowvals[kk].real!=0.0||rowvals[kk].imag!=0.0)
      {
        // the number of nonzero elements increases by 1
        nnz++;
        // add in the value
        val[nnz]=rowvals[kk];
        col[nnz]=kk;
      }
    
  }


  // COPY THE TEMP MATRICES TO THE FULL CSR MATRIX
  csrmat.init(nnz,ninterpol_pts_assigned);

  for(k=1;k<=ninterpol_pts_assigned;++k)
    csrmat.rowpntr[k]=rowpntr[k];
 
  for(k=1;k<=nnz;++k)
  {
    csrmat.val[k]=val[nnz];
    csrmat.col[k]=col[nnz];
  } 

  // CLEANUP
  delete [] val;
  delete [] rowvals;
  delete [] col;
  delete [] rowpntr;

  return;
}




// j7
void convert_matrix(MATRIX_SPARSE_CD smatrix, MATRIXCSR_CD*& csrmat)
{

  int k,kk;

  complex<double> COMPLEXZERO(0.0,0.0);

  // SET UP A TEMP ARRAY TO HOLD ROW VALUES
  complex<double> *rowvals;
  rowvals=new complex<double> [ninterpol_pts_assigned+1];  

  // SET UP A TEMP CSR MATRIX TO MAXIMUM POSSIBLE SIZE
  // i. i.e. this size will be correct if there is no duplicate position values in smatrix
  MATRIXCSR_CD *csrmat_temp;
  csrmat_temp=new MATRIXCSR_CD(ninterpol_pts_assigned, smatrix.size);
  csrmat_temp->nzeros=smatrix.size;

  int nnz=0;


  // LOOP OVER ROWS
  for(k=1;k<=ninterpol_pts_assigned;++k)
  {
    // RESET ARRAY
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
    {
      rowvals[kk]=COMPLEXZERO;
    }

    // FILL ROWVALS ARRAY
    for(kk=1;kk<=smatrix.size;++kk)
    {
      if(smatrix.ia[kk]==k)
      {
        rowvals[smatrix.ja[kk]]+= complex<double>(smatrix.a[kk].real, smatrix.a[kk].imag);
      }
    }

    // SET ROWPNTR TO POINT TO 1 HIGHER THAN CURRENT NONZERO ELEMENTS
    // i. ! each row bettwe contain at least 1 nonzero element
    csrmat_temp->rowpntr[k]=nnz+1;

    // loop over "columns"
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
    {
      if(rowvals[kk]!=COMPLEXZERO)
      {
        // INCREMENT NONZERO ELEMENTS
        nnz++;

        // ADD VALUE
        csrmat_temp->val[nnz]=rowvals[kk];
        csrmat_temp->col[nnz]=kk;
      }
    }
    
  } // ++k

  // this is used in loops
  csrmat_temp->rowpntr[ninterpol_pts_assigned+1]=nnz+1;


  //=========================================================
  // COPY THE TEMP CSR MATRIX TO THE REAL CSR MATRIX
  //=========================================================
  csrmat=new MATRIXCSR_CD(ninterpol_pts_assigned, nnz);
  csrmat->nzeros=nnz;

  for(k=1;k<=nnz;++k)
  {
    csrmat->val[k]=csrmat_temp->val[k];
    csrmat->col[k]=csrmat_temp->col[k];
  }

  for(k=1;k<=ninterpol_pts_assigned+1;++k)
  {
    csrmat->rowpntr[k]=csrmat_temp->rowpntr[k];
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] rowvals;

  delete csrmat_temp;


  return;
}



void convert_matrix(MATRIX_SPARSE_CD smatrix, MATRIXCSC_CD*& cscmat)
{

  int k,kk;

  complex<double> COMPLEXZERO(0.0,0.0);

  // SET UP A TEMP ARRAY TO HOLD COLUMN VALUES
  complex<double> *colvals;
  colvals=new complex<double> [ninterpol_pts_assigned+1];  

  // SET UP A TEMP CSC MATRIX TO MAXIMUM POSSIBLE SIZE
  // i. i.e. this size will be correct if there is no duplicate position values in smatrix
  MATRIXCSC_CD *cscmat_temp;
  cscmat_temp=new MATRIXCSC_CD(ninterpol_pts_assigned, smatrix.size);
  cscmat_temp->nzeros=smatrix.size;
  int nnz=0;


  // LOOP OVER COLUMNS
  for(k=1;k<=ninterpol_pts_assigned;++k)
  {
    // RESET ARRAY
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
    {
      colvals[kk]=COMPLEXZERO;
    }

    // FILL COLVALS ARRAY
    for(kk=1;kk<=smatrix.size;++kk)
    {
      if(smatrix.ja[kk]==k)
      {
        colvals[smatrix.ia[kk]]+= complex<double>(smatrix.a[kk].real, smatrix.a[kk].imag);
      }
    }

    // SET ROWPNTR TO POINT TO 1 HIGHER THAN CURRENT NONZERO ELEMENTS
    // i. ! each row bettwe contain at least 1 nonzero element
    cscmat_temp->colpntr[k]=nnz+1;

    // loop over "columns"
    for(kk=1;kk<=ninterpol_pts_assigned;++kk)
    {
      if(colvals[kk]!=COMPLEXZERO)
      {
        // INCREMENT NONZERO ELEMENTS
        nnz++;

        // ADD VALUE
        cscmat_temp->val[nnz]=colvals[kk];
        cscmat_temp->row[nnz]=kk;
      }
    }
    
  } // ++k

  // this is used in loops
  cscmat_temp->colpntr[ninterpol_pts_assigned+1]=nnz+1;

  //=========================================================
  // COPY THE TEMP CSC MATRIX TO THE REAL CSR MATRIX
  //=========================================================
  cscmat=new MATRIXCSC_CD(ninterpol_pts_assigned, nnz);
  cscmat->nzeros=nnz;
  for(k=1;k<=nnz;++k)
  {
    cscmat->val[k]=cscmat_temp->val[k];
    cscmat->row[k]=cscmat_temp->row[k];
  }

  for(k=1;k<=ninterpol_pts_assigned+1;++k)
  {
    cscmat->colpntr[k]=cscmat_temp->colpntr[k];
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] colvals;

  delete cscmat_temp;


  return;
}


// i. this assumes a square matrix
void convert_matrix(MATRIXCSC_CD*& origmat, MATRIXCSR_CD*& convertmat)
{

  int j, k, kk;


  convertmat = new MATRIXCSR_CD(origmat->ncols, origmat->nzeros);
  convertmat->nzeros=origmat->nzeros;

  complex<double> *colvals, COMPLEXZERO(0.0,0.0);
  colvals = new complex<double> [ninterpol_pts_assigned+1];

  int nnz=0;

  // LOOP OVER ROWS
  for(k=1;k<=ninterpol_pts_assigned;++k)
  {
    //if(origmat->row[origmat->colpntr[k]]==

    // RESET colvals[] 
    for(j=1;j<=ninterpol_pts_assigned;++j)
    {
      colvals[j]=COMPLEXZERO;
    }  


    for(kk=1;kk<=origmat->nzeros;++kk)
    {
      if(origmat->row[kk]==k)
      {
        // FIND OUT WHAT COLUMN WE ARE IN
        for(j=1;j<=ninterpol_pts_assigned;++j)
        {
          if((kk>=origmat->colpntr[j]) && (kk<origmat->colpntr[j+1]))
          {
            colvals[j]=origmat->val[kk];
            break;
          }
        }
      }
    } // ++kk

    // INCREMENT THE ROW POINTER FIRST
    convertmat->rowpntr[k]=nnz+1;

    // ADD THIS TO OUT TRANSPOSE MATRIX
    for(j=1;j<=ninterpol_pts_assigned;++j)
    {
      if(colvals[j]!=COMPLEXZERO)
      {
        nnz+=1;
        convertmat->val[nnz]=colvals[j];
        convertmat->col[nnz]=j;
      }
    }    

  }  // ++k


  // INCREMENT THE LAST OF THE ROW POINTER
  convertmat->rowpntr[ninterpol_pts_assigned+1]=nnz+1;


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] colvals;


  return;
}

void convert_matrix(MATRIXCSR_CD*& origmat, MATRIXCSC_CD*& convertmat)
{

  int j, k, kk;


  convertmat = new MATRIXCSC_CD(origmat->nrows, origmat->nzeros);
  convertmat->nzeros=origmat->nzeros;

  complex<double> *rowvals, COMPLEXZERO(0.0,0.0);
  rowvals = new complex<double> [ninterpol_pts_assigned+1];

  int nnz=0;

  // LOOP OVER COLUMNS
  for(k=1;k<=ninterpol_pts_assigned;++k)
  {

    // RESET rowvals[] 
    for(j=1;j<=ninterpol_pts_assigned;++j)
    {
      rowvals[j]=COMPLEXZERO;
    }  


    for(kk=1;kk<=origmat->nzeros;++kk)
    {
      if(origmat->col[kk]==k)
      {
        // FIND OUT WHAT COLUMN WE ARE IN
        for(j=1;j<=ninterpol_pts_assigned;++j)
        {
          if((kk>=origmat->rowpntr[j]) && (kk<origmat->rowpntr[j+1]))
          {
            rowvals[j]=origmat->val[kk];
            break;
          }
        }
      }
    } // ++kk




    // INCREMENT THE ROW POINTER FIRST
    convertmat->colpntr[k]=nnz+1;

    // ADD THIS TO OUT TRANSPOSE MATRIX
    for(j=1;j<=ninterpol_pts_assigned;++j)
    {
      if(rowvals[j]!=COMPLEXZERO)
      {
        nnz+=1;
        convertmat->val[nnz]=rowvals[j];
        convertmat->row[nnz]=j;
      }
    }    

  }  // ++k


  // INCREMENT THE LAST OF THE ROW POINTER
  convertmat->colpntr[ninterpol_pts_assigned+1]=nnz+1;


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] rowvals;


  return;
}

*/

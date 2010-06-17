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
// NAME: filecombine:: A FILE COMBINATION PROGRAM FOR JFDTD
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: filecombine.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@u.northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************



//************************************************************************
//			INCLUDE FILES
//************************************************************************	
	
// I/O
#include <cstring>
#include <cstdlib>
// for screen I/O
#include <iostream>
// for file I/O	
#include <fstream>


// EXTERNAL MATH ROUTINES
#include <cmath>

using namespace std;	



//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************



//************************************************************************
//			SUBROUTINES
//************************************************************************

// INPUT
void input();


// ARRAY SORTING
void sort_field_arrays();


// INPUT
void output(int ioutput_type, int iplane, int iunit);
void output_dx(int iplane, int iunit);


//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

int grid_nx, grid_ny, grid_nz;
double grid_dx, grid_dy, grid_dz;

int mpi_nxprocs, mpi_nyprocs, mpi_nzprocs;

int noutput, output_nplanes, *output_plane, *output_format;


double *field_pos1, *field_pos2, *field_val;

double *field_pos1_small, *field_pos2_small, *field_val_small;


// this will hold the size of the field arrays
int field_array_pt, field_array_pt_small;



//========================================================================
//========================================================================
//
// 	NAME:	int main(int argc, char** argv) 
//	DESC:	Entryway to program
//
//	NOTES:
//
//========================================================================
//========================================================================
int main(int argc, char** argv)
{

  // local indices
  int i, p, n, np1, np2;


  //=========================================================
  // READ IN PARAMETERS
  //=========================================================
 

  // INPUT THE FIELD STITCH INFORMATION
  input();


  //=========================================================
  // READ IN FIELDS
  //=========================================================

  // ALLOCATE MEMORY FOR FIELD POINTS

  // set maximum number of field points
  // i. a max size of 10000x10000 is probably plenty large for a plane 
  // from fdtd simulations
  int max_field_points = 10000*10000;

  field_pos1 = new double [max_field_points + 1]; 
  field_pos2 = new double [max_field_points + 1]; 
  field_val = new double [max_field_points + 1];


  field_pos1_small = new double [max_field_points + 1]; 
  field_pos2_small = new double [max_field_points + 1]; 
  field_val_small = new double [max_field_points + 1];

  int nprocs1, nprocs2;


  // standard filename
  char filename_1[] = "./output/e_field.";
	

  //=========================================================
  // LOOP OVER ALL PLANES
  //=========================================================
  // i. this will be the first number (i.e. e_field.#p...)
  // ii. this number ranges from 1 to number of planes we set to output
  // in the FDTD sim (output_nplanes)
  for(p = 1; p <= output_nplanes; ++p)
  {

    // OpenDX
    if(output_format[p] == 2)
    {
      for(n = 1; n <= noutput; ++n)
      {
        output_dx(p, n);
      }
      continue;
    }

    // SET UP PLANE PART OF FILENAME
    char filename_2[5];
    sprintf(filename_2,"%i",p); 

 
    //---------------------------------------------------------
    // GET NUMBER OF TOTAL OUTPUT FILES FOR PLANE
    //---------------------------------------------------------
    // i. this depends on the number of processors since each processor
    // outputs its own file

    // xy
    if(output_plane[p] == 1) 
    {
      nprocs1 = mpi_nxprocs;
      nprocs2 = mpi_nyprocs;
    }
    // xz
    else if(output_plane[p] == 2)
    {
      nprocs1 = mpi_nxprocs;
      nprocs2 = mpi_nzprocs;
    }
    // yz
    else if(output_plane[p] == 3)
    {
      nprocs1 = mpi_nyprocs;
      nprocs2 = mpi_nzprocs; 
    }

    cout << "nproc1: " << nprocs1 << endl;
    cout << "nproc2: " << nprocs2 << endl;


    //=========================================================
    // LOOP OVER ALL OUTPUT NUMBERS 
    //=========================================================
    // i. this will be the second number (i.e. e_field.#p.#n...)
    for(n = 1; n <= noutput; ++n)
    {


      // reset the point in our field arrays (total)
      field_array_pt = 0;


      // SET UP noutput PART OF FILENAME
      char filename_3[5];
      sprintf(filename_3,"%i",n);



      //=========================================================
      // LOOP OVER ALL OUTPUT FILES
      //=========================================================
      // i. we loop over all of the processors

      for(np1 = 0; np1 < nprocs1; ++np1)
      {

        // CONTINUE SETTING UP FILE
        // number proc1
        char filename_4[5];
        sprintf(filename_4,"%i",np1);


        // reset the point in our field arrays (small)
        field_array_pt_small = 0;


        for(np2 = 0; np2 < nprocs2; ++np2)
        {  

          // TOTALLY FINISH SETTING UP FILE

          // number proc2
          char filename_5[5];
          sprintf(filename_5,"%i",np2);

          // total filename
          char filename[31];

          // now copy all filenames together
          strcpy(filename, filename_1);
          strcat(filename, filename_2);
          strcat(filename, ".");
          strcat(filename, filename_3);
          strcat(filename, ".");
          strcat(filename, filename_4);
          strcat(filename, ".");
          strcat(filename, filename_5);


          //---------------------------------------------------------
          // OPEN FILE
          //---------------------------------------------------------

          // create file variable
          ifstream field_file;

          // open file
          field_file.open(filename, ios::in);


          //---------------------------------------------------------
          // READ FILE
          //---------------------------------------------------------
          // i. here we read in values until we hit the end of the file
          // ii. !!! I don't know which of the following ways is better of reading in a file

          double pos1;

          while(field_file >> pos1) 
          {

            // increase the point in our field array
            field_array_pt_small += 1;

            // read the pos1 (x) value
            field_pos1_small[field_array_pt_small] = pos1;

            // read the pos2 (y) value
            field_file >> field_pos2_small[field_array_pt_small];

            // read the value of the field here
            field_file >> field_val_small[field_array_pt_small];
          }


/*
          while(!field_file.eof()) 
          {

            // increase the point in our field array
            field_array_pt_small += 1;

            // read the pos1 (x) value
            field_file >> field_pos1_small[field_array_pt_small];

            // read the pos2 (y) value
            field_file >> field_pos2_small[field_array_pt_small];

            // read the value of the field here
            field_file >> field_val_small[field_array_pt_small];
          }
*/

          //---------------------------------------------------------
          // CLOSE FILE
          //---------------------------------------------------------
          field_file.close();


        } // ++np2 (number of procs 1)

        //---------------------------------------------------------
        // SORT FIELDS (SMALL ARRAY)
        //---------------------------------------------------------
        sort_field_arrays();

        //---------------------------------------------------------
        // ADD SMALL ARRAY TO TOTAL ARRAY
        //---------------------------------------------------------
        for(i = 1; i <= field_array_pt_small; ++i)
        {
          // increase the point in our field array
          field_array_pt += 1;

          // read the pos1 (x) value
          field_pos1[field_array_pt] = field_pos1_small[i];

          // read the pos2 (y) value
          field_pos2[field_array_pt] = field_pos2_small[i];

          // read the value of the field here
          field_val[field_array_pt] = field_val_small[i];

        }


      } // ++np1 (number of procs 1)


      //---------------------------------------------------------
      // OUTPUT TOTAL FIELD
      //---------------------------------------------------------
      output(1, p, n);


    } // ++n (next time output step)


  } // ++p (next plane)



  //=========================================================
  // CLEANUP
  //=========================================================

  delete [] field_pos1; 
  delete [] field_pos2; 
  delete [] field_val;


  delete [] field_pos1_small; 
  delete [] field_pos2_small; 
  delete [] field_val_small;



  return 0;
}


//========================================================================
//========================================================================
//
// 	NAME:	void input(int iinput_type)
//	DESC:	Reads in params
//
//	NOTES:
//
//========================================================================
//========================================================================
void input()
{

    //=========================================================
    // STITCH FILE
    //=========================================================

    //---------------------------------------------------------
    // OPEN FILE
    //---------------------------------------------------------
    // create file variable
    ifstream stitch_file;

    // open file
    stitch_file.open("./output/field_stitch.dat", ios::in);


    //---------------------------------------------------------
    // READ IN PARAMS
    //---------------------------------------------------------

    // read in number of processors
    stitch_file >> mpi_nxprocs;    
    stitch_file >> mpi_nyprocs; 
    stitch_file >> mpi_nzprocs; 


    // READ IN NUMBER OF GRID POINTS
    stitch_file >> grid_nx;    
    stitch_file >> grid_ny; 
    stitch_file >> grid_nz; 

    // READ IN GRID SIZE
    stitch_file >> grid_dx;    
    stitch_file >> grid_dy; 
    stitch_file >> grid_dz;  

    // read in number of output files per plane
    // !!! eventually this will be changed to not have ft every time step
    stitch_file >> noutput;

    // write number of planes
    stitch_file >> output_nplanes;


    // NOW ALLOCATE A LITTLE MEMORY
    // !!! bad mem
    output_format = new int [output_nplanes+1]; 
    output_plane = new int [output_nplanes+1]; 

    // now loop over the planes
    for(int p = 1; p <= output_nplanes; ++p)
    {
      // INPUT THE FORMAT 
      stitch_file >> output_format[p];

      // INPUT THE TYPE OF PLANE
      // i.) Types:
      //	1 == xy
      //	2 == xz
      //	3 == yz
      stitch_file >> output_plane[p];
    }


    //---------------------------------------------------------
    // CLOSE FILE
    //---------------------------------------------------------
    stitch_file.close();


  return;
}






//========================================================================
//========================================================================
//
// 	NAME:	void output(int ioutput_type, int iplane, int iunit)
//	DESC:	output total field
//
//	NOTES:
//
//========================================================================
//========================================================================
void output(int ioutput_type, int iplane, int iunit)
{

  int i;


  //=========================================================
  // SET-UP FILE (FILENAME, ETC)
  //=========================================================

  char filename_1[] = "./output/e_field_stitched.";
  char filename_2[5];
  char filename_3[5];

  char filename[31];	

  // SET UP iplane AND iunit
  sprintf(filename_2,"%i",iplane); 
  sprintf(filename_3,"%i",iunit);
   
  // NOW COPY ALL FILENAMES TOGETHER
  strcpy(filename, filename_1);
  strcat(filename, filename_2);
  strcat(filename, ".");
  strcat(filename, filename_3);

  // OpenDX specific  
  if(output_format[iplane] == 2)
  {
    strcat(filename, ".dx");
  }

  // OPEN FILE
  ofstream output_file(filename,ios::app);


  //=========================================================
  // WRITE FILE
  //=========================================================

  // loop over each value
  for(i = 1; i <= field_array_pt; ++i)
  {

    // if we are at a new "x" position leave an extra space (for pm3d in gnuplot)
    if( (field_pos1[i] != field_pos1[i-1]) && (i != 1) )
    {
      output_file << endl;
    }

    // output the field here
    output_file << field_pos1[i] << "     " << field_pos2[i] << "     " << field_val[i] << endl;
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================

  output_file.close();

  return;
}


void output_dx(int iplane, int iunit)
{

  int i, ii, j, jj, k, kk;
  double field;

  //=========================================================
  // SET-UP FILE (FILENAME, ETC)
  //=========================================================

  char filename_1[] = "./output/e_field_stitched.";
  char filename_2[5];
  char filename_3[5];

  char filename[31];	

  // SET UP iplane AND iunit
  sprintf(filename_2,"%i",iplane); 
  sprintf(filename_3,"%i",iunit);
   
  // NOW COPY ALL FILENAMES TOGETHER
  strcpy(filename, filename_1);
  strcat(filename, filename_2);
  strcat(filename, ".");
  strcat(filename, filename_3);

  // OpenDX specific  
  strcat(filename, ".dx");

  // OPEN FILE
  ofstream output_file(filename,ios::app);

  // standard filename
  int *minx, *maxx, *miny, *maxy, *minz, *maxz;
  minx = new int [mpi_nxprocs];
  maxx = new int [mpi_nxprocs];
  miny = new int [mpi_nxprocs];
  maxy = new int [mpi_nxprocs];
  minz = new int [mpi_nxprocs];
  maxz = new int [mpi_nxprocs];

  ifstream ***input_file;
  char filename_in1[] = "./output/e_field.";
  char filename_in2[5];
  char filename_in3[5];
  char filename_in4[5];

  input_file = new ifstream**[mpi_nxprocs];
  for(i = 0; i < mpi_nxprocs; ++i)
  {
    input_file[i] = new ifstream*[mpi_nyprocs];
    sprintf(filename_in2,"%i",i);
    for(j = 0; j < mpi_nyprocs; ++j)
    {
      input_file[i][j] = new ifstream [mpi_nzprocs];
      sprintf(filename_in3,"%i",j);
      for(k = 0; k < mpi_nzprocs; ++k)
      {
        sprintf(filename_in4,"%i",k);
      
        // NOW COPY ALL FILENAMES TOGETHER
        strcpy(filename, filename_in1);
        strcat(filename, filename_2);
        strcat(filename, ".");
        strcat(filename, filename_3);
        strcat(filename, ".");
        strcat(filename, filename_in2);
        strcat(filename, ".");
        strcat(filename, filename_in3);
        strcat(filename, ".");
        strcat(filename, filename_in4);

        // OPEN FILE
        input_file[i][j][k].open(filename,ios::app);

        
        input_file[i][j][k] >> minx[i];
        input_file[i][j][k] >> maxx[i];
        input_file[i][j][k] >> miny[j];
        input_file[i][j][k] >> maxy[j];
        input_file[i][j][k] >> minz[k];
        input_file[i][j][k] >> maxz[k];
      }
    }
  }

/*
  int *minx, *maxx, *miny, *maxy, *minz, *maxz;
  minx = new int [mpi_nxprocs];
  maxx = new int [mpi_nxprocs];
  miny = new int [mpi_nxprocs];
  maxy = new int [mpi_nxprocs];
  minz = new int [mpi_nxprocs];
  maxz = new int [mpi_nxprocs];


  int xnum, ynum, znum;

  // LOOP OVER ALL x POINTS
  for(i = 1; i <= grid_nx; ++i)
  {
    // GET THE x PROCESSOR NUMBER
    for(ii = 0; ii < mpi_nxprocs; ++ii)
    {
      if( (i >= minx[ii]) && (i <= maxx[ii]) )
      {
        xnum = ii;
        break;
      }
    }

    // LOOP OVER ALL y POINTS
    for(j = 1; j <= grid_ny; ++j)
    {

      // GET THE y PROCESSOR NUMBER
      for(jj = 0; jj < mpi_nyprocs; ++jj)
      {
        if( (j >= miny[jj]) && (j <= maxy[jj]) )
        {
          ynum = jj;
          break;
        }
      }

      // LOOP OVER ALL z POINTS
      for(k = 1; k <= grid_nz; ++k)
      {
    
        // GET THE z PROCESSOR NUMBER
        for(kk = 0; kk < mpi_nzprocs; ++kk)
        {
          if( (k >= minz[kk]) && (k <= maxz[kk]) )
          {
            znum = kk;
            break;
          }
        }

        


      } // ++k
    } // ++j
  } // ++i
  */


  //=========================================================
  // READ / WRITE FILES
  //=========================================================

  //---------------------------------------------------------
  // INITIAL OUTPUT
  //---------------------------------------------------------
  output_file << "object 1 class gridpositions counts " << grid_nx << " " << grid_ny << " " << grid_nz << endl;
  output_file << "origin 0.0 0.0 0.0" << endl;    
  output_file << "delta " << grid_dx << " 0.0 0.0" << endl;    
  output_file << "delta 0.0 " << grid_dy << " " << "0.0" << endl; 
  output_file << "delta 0.0 0.0 " << grid_dz << endl;
  output_file << "object 2 class gridconnections counts " << grid_nx << "   " << grid_ny << " " << grid_nz << endl; 
  output_file << "object 3 class array type float rank 0 items " << grid_nx*grid_ny*grid_nz << " data follows" << endl; 

  //---------------------------------------------------------
  // MAIN DATA
  //---------------------------------------------------------

  // LOOP OVER ALL x POINTS
  for(i = 0; i < mpi_nxprocs; ++i)
  {
    // GET THE x PROCESSOR NUMBER
    for(ii = minx[i]; ii <= maxx[i]; ++ii)
    {
      // LOOP OVER ALL y POINTS
      for(j = 0; j < mpi_nyprocs; ++j)
      {
        for(jj = miny[j]; jj <= maxy[j]; ++jj)
        {
          // LOOP OVER ALL y POINTS
          for(k = 0; k < mpi_nzprocs; ++k)
          {
            for(kk = minz[k]; kk <= maxz[k]; ++kk)
            {
              // READ IN THE FIELD
              input_file[i][j][k] >> field;

              // WRITE OUT THE FIELD
              output_file << field;
            } // ++kk
          } // ++k
        } // ++jj
      } // ++j
    } // ++ii
  } // ++i



  //---------------------------------------------------------
  // CLOSING OUTPUT
  //---------------------------------------------------------
  output_file << "object \"regular positions regular connections\" class field" << endl;
  output_file << " component \"positions\" value 1" << endl;    
  output_file << " component \"connections\" value 2" << endl; 
  output_file << " component \"data\" value 3" << endl; 
  output_file << "end" << endl; 


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================

  delete [] minx;
  delete [] maxx;
  delete [] miny;
  delete [] maxy;
  delete [] minz;
  delete [] maxz;

  for(i = 0; i < mpi_nxprocs; ++i)
  {
    for(j = 0; j < mpi_nyprocs; ++j)
    {
      for(k = 0; k < mpi_nzprocs; ++k)
      {
        input_file[i][j][k].close();
      }

      delete [] input_file[i][j];
    }

    delete [] input_file[i];
  }

  delete [] input_file;

  output_file.close();

  return;
}



//========================================================================
//========================================================================
//
// 	NAME:	void sort_field_arrays()
//	DESC:	output total field
//
//	NOTES:
//
//========================================================================
//========================================================================
void sort_field_arrays()
{

  int i, ii;

  int lowest_pos;

  // for swapping
  double temp_val;



  //====================================================
  // SORT THE pos1 VALUES
  //====================================================
  // i.) field_array_pt is the size

  // loop over each field point, except the last one
  for(i = 1; i <= (field_array_pt_small-1); ++i)
  {
    // set the lowest position to this i value
    lowest_pos = i;

    // now loop over all higher values looking for a smaller pos1 value
    for(ii = (i+1); ii <= field_array_pt_small; ++ii)
    {
      // if we found a lower value mark it
      if(field_pos1_small[ii] < field_pos1_small[lowest_pos])
      {
        lowest_pos = ii;
      }
    }

    // now swap the lowest position with the i position
    // i. this lowest position may actually be i, so we may swap i with
    // itself

    temp_val = field_pos1_small[i];
    field_pos1_small[i] = field_pos1_small[lowest_pos];
    field_pos1_small[lowest_pos] = temp_val;

    temp_val = field_pos2_small[i];
    field_pos2_small[i] = field_pos2_small[lowest_pos];
    field_pos2_small[lowest_pos] = temp_val;

    temp_val = field_val_small[i];
    field_val_small[i] = field_val_small[lowest_pos];
    field_val_small[lowest_pos] = temp_val;

  }

  //====================================================
  // NOW SORT THE pos2 VALUES
  //====================================================

  // loop over each field point, except the last one
  for(i = 1; i <= (field_array_pt_small-1); ++i)
  {
    // set the lowest position to this i value
    lowest_pos = i;

    // set ii to be 1 higher than i
    ii = i+1;

    // loop while pos1 values are the same
    while(field_pos1_small[ii] == field_pos1_small[i])
    {

      // if we found a lower value mark this value
      if(field_pos2_small[ii] < field_pos2_small[lowest_pos])
      {
        lowest_pos = ii;
      }

      // increment ii
      ii += 1;

    }
    
    // now swap the lowest position with the i position
    // i. this lowest position may actually be i, so we may swap i with
    // itself
    temp_val = field_pos1_small[i];
    field_pos1_small[i] = field_pos1_small[lowest_pos];
    field_pos1_small[lowest_pos] = temp_val;

    temp_val = field_pos2_small[i];
    field_pos2_small[i] = field_pos2_small[lowest_pos];
    field_pos2_small[lowest_pos] = temp_val;

    temp_val = field_val_small[i];
    field_val_small[i] = field_val_small[lowest_pos];
    field_val_small[lowest_pos] = temp_val;

  }



  return;
}



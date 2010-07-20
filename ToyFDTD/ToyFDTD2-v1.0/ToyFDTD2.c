// ToyFDTD2, version 1.0 
//      The if-I-can-do-it-you-can-do-it FDTD! 
// Copyright (C) 1998,1999 Laurie E. Miller, Paul Hayes, Matthew O'Keefe
// Copyright (C) 1999 John Schneider 

// This program is free software; you can redistribute it and/or 
//     modify it under the terms of the GNU General Public License 
//     as published by the Free Software Foundation; either version 2
//     of the License, or any later version, with the following conditions
//     attached in addition to any and all conditions of the GNU
//     General Public License:
//     When reporting or displaying any results or animations created
//     using this code or modification of this code, make the appropriate
//     citation referencing ToyFDTD2 by name and including the version
//     number.  
//
// This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty 
//     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//     See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  
//     02111-1307  USA

// Contacting the perpetratorss:
//
// Laurie E. Miller, Paul Hayes, Matthew O'Keefe
// Department of Electrical and Computer Engineering
//      200 Union Street S. E.
//      Minneapolis, MN 55455
//
// lemiller@borg.umn.edu
// 
// http://www.borg.umn.edu/toyfdtd/
// http://www.borg.umn.edu/toyfdtd/ToyFDTD2.html
// http://www.toyfdtd.org/

// This code is here for everyone, but not everyone will need something 
//      so simple, and not everyone will need to read all the comments.  
// This file is over 700 lines long, but less than 400 of that is actually
//      code. 

///////////////////////////////////////////////////////////////////////////////
// ToyFDTD2 builds into ToyFDTD1 a way of allocating memory so that each 
//      field verctor array is forced to be be contiguous, which the
//      ToyFDTD1 scheme does not.  In many cases this will be faster
//      than the ToyFDTD1 scheme, but in some cases it will be slower;
//      which will be better is highly platform-dependent.  Some optimizing
//      compilers simply force the ToyFDTD1 arrays to be contiguous.  It's
//      worth testing more than one option before undertaking a large run.
//      Both schemes are designed to be easy to read.  In ToyFDTD2, macros
//      are defined for accessing the arrays to make the notation easy to 
//      read: i.e. Ez(i,j,k).  The modifications were written by 
//      John Schneider.
// This ToyFDTD2 is a stripped-down, minimalist 3D FDTD code.  It 
//      illustrates the minimum factors that must be considered to 
//      create a simple FDTD simulation.  
///////////////////////////////////////////////////////////////////////////////

// This is a very simple Yee algorithm 3D FDTD code in C implementing
//      the free space form of Maxwell's equations on a Cartesian grid.
// There are no internal materials or geometry.  
// The code as delivered simulates an idealized rectangular waveguide 
//      by treating the interior of the mesh as free space/air and enforcing
//      PEC (Perfect Electric Conductor) conditions on the faces of the mesh.
// The problem is taken from Field and Wave Electromagnetics, 2nd ed., by 
//      David K. Cheng, pages 554-555.  It is a WG-16 waveguide 
//      useful for X-band applications, interior width = 2.29cm, 
//      interior height = 1.02cm.  The frequency (10 GHz) is chosen to be 
//      in the middle of the frequency range for TE10 operation.  
// Boundaries: PEC (Perfect Electric Conductor).
// Stimulus: A simplified sinusoidal plane wave emanates from x = 0 face.
// 3D output: The electric field intensity vector components in the direction 
//      of the height of the guide (ez) are output to file every 
//      PLOT_MODULUS timesteps (and for the last timestep), scaled to 
//      the range of integers from zero through 254.  (The colormap 
//      included in the tar file assigns rgb values to the range zero through
//      255.)   Scaling is performed on each timestep individually, since 
//      it is not known in advance what the maximum and minimum 
//      values will be for the entire simulation.  The integer value 127 is 
//      held to be equal to the data zero for every timestep.  This method 
//      of autoscaling every timestep can be very helpful in a simulation 
//      where the intensities are sometimes strong and sometimes faint, 
//      since it will highlight the presence and structure of faint signals 
//      when stronger signals have left the mesh.  
//      Each timestep has it's own output file.  This data output file format
//      can be used in several visualization tools, such as animabob and viz. 
// Other output: Notes on the progress of the simulation are written to standard
//      output as the program runs.  
//      A .viz file is output to feed parameters to viz, should viz later be
//      used to view the data files.

// Some terminology used here:
//
// This code implements a Cartesian mesh with space differentials 
//     of dx, dy, dz.
// This means that a point in the mesh has neighboring points dx meters 
//     away in the direction of the positive and negative x-axis,
//     and neighboring points dy meters away in the directions 
//     of the +- y-axis, and neighboring points dz meters away 
//     in the directions of the +- z-axis,
// The mesh has nx cells in the x direction, ny in the y direction, 
//     and nz in the z direction.
// ex, ey, and ez refer to the arrays of electric field intensity vectors 
//     -- for example, ex is a 3-dimensional array consisting of the 
//     x component of the E field intensity vector for every point in the 
//     mesh.  ex[i][j][k] refers to the x component of the E field intensity 
//     vector at point [i][j][k].  
// hx, hy, and hz refer to the arrays of magnetic field intensity vectors.
//
// dt is the time differential -- the length of each timestep in seconds.
//
// bob is a file format that stands for "brick of bytes", meaning a string 
//     of bytes that can be interpreted as a 3-dimensional array of byte 
//     values (integers from zero through 255).  animabob is a free 
//     visualization tool that displays and animates a sequence of bricks 
//     of bytes.  For more information on animabob or to download a copy, 
//     see the ToyFDTD website at http://www.borg.umn.edu/toyfdtd/
//
// viz is another free visualization tool that displays and animates 
//     brick-of-byte files. For more information on viz or to download a copy, 
//     see the ToyFDTD website at http://www.borg.umn.edu/toyfdtd/

#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

// program control constants
#define MAXIMUM_ITERATION 1000
          // total number of timesteps to be computed
#define PLOT_MODULUS 5
          // The program will output 3D data every PLOT_MODULUS timesteps,
          //     except for the last iteration computed, which is always
          //     output.  So if MAXIMUM_ITERATION is not an integer
          //     multiple of PLOT_MODULUS, the last timestep output will
          //     come after a shorter interval than that separating
          //     previous outputs.  
#define FREQUENCY 10.0e9     
          // frequency of the stimulus in Hertz
#define GUIDE_WIDTH 0.0229
          // meters
#define GUIDE_HEIGHT 0.0102
          // meters
#define LENGTH_IN_WAVELENGTHS 5.0
          // length of the waveguide in wavelengths of the stimulus wave
#define CELLS_PER_WAVELENGTH 25.0
          // minimum number of grid cells per wavelength in the x, y, and
          //     z directions

// physical constants
#define LIGHT_SPEED 299792458.0       
          // speed of light in a vacuum in meters/second
#define LIGHT_SPEED_SQUARED 89875517873681764.0        
          // m^2/s^2
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
          // permeability of free space in henry/meter
#define EPSILON_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
          // permittivity of free space in farad/meter

// The value used for pi is M_PI as found in /usr/include/math.h on SGI 
//     IRIX 6.2.  Other such constants used here are DBL_EPSILON and FLT_MAX

main()	
{// main

  /////////////////////////////////////////////////////////////////////////////
  // variable declarations
  int i,j,k;     
           // indices of the 3D array of cells
  int ioff_ex, ioff_ey, ioff_ez, ioff_hx, ioff_hy, ioff_hz;
           // offsets used in the macros below for calculating position
           // within the data arrays
  int nx, ny, nz;         
           // total number of cells along the x, y, and z axes, respectively
  int allocatedBytes = 0;          
           // a counter to track number of bytes allocated

  int iteration = 0;          
           // counter to track how many timesteps have been computed
  double stimulus = 0.0;       
           // value of the stimulus at a given timestep
  double currentSimulatedTime = 0.0;
           // time in simulated seconds that the simulation has progressed
  double totalSimulatedTime = 0.0;       
           // time in seconds that will be simulated by the program
  double omega;       
           // angular frequency in radians/second
  double lambda;       
           // wavelength of the stimulus in meters
  double dx, dy, dz; 
           // space differentials (or dimensions of a single cell) in meters
  double dt;
           // time differential (how much time between timesteps) in seconds
  double dtmudx, dtepsdx;       
           // physical constants used in the field update equations 
  double dtmudy, dtepsdy;       
           // physical constants used in the field update equations 
  double dtmudz, dtepsdz;       
           // physical constants used in the field update equations 

  double *ex, *ey, *ez; 
           // pointers to the arrays of ex, ey, and ez values
#define Ex(I,J,K) ex[ioff_ex*(I) + (nz+1)*(J) + (k)]
#define Ey(I,J,K) ey[ioff_ey*(I) + (nz+1)*(J) + (k)]
#define Ez(I,J,K) ez[ioff_ez*(I) + nz*(J) + (k)]
           // macros to make accessing Ex, Ey, and Ez convenient
           // i.e. Ez(i,j,k)
  double *hx, *hy, *hz;
           // pointers to the arrays of hx, hy, and hz values
#define Hx(I,J,K) hx[ioff_hx*(I) + nz*(J) + (k)]
#define Hy(I,J,K) hy[ioff_hy*(I) + nz*(J) + (k)]
#define Hz(I,J,K) hz[ioff_hz*(I) + (nz-1)*(J) + (k)]
           // macros to make accessing Hx, Hy, and Hz convenient
           // i.e. Hz(i,j,k)

  // bob output routine variables:
  char filename[1024];
           // filename variable for 3D bob files
  double simulationMin = FLT_MAX;
           // tracks minimum value output by the entire simulation
  double simulationMax = -FLT_MAX;
           // tracks maximum value output by the entire simulation
  double min, max;
           // these track minimum and maximum values output in one timestep
  double norm;
           // norm is set each iteration to be max or min, whichever is 
           //     greater in magnitude
  double scalingValue;
           // multiplier used in output scaling, calculated every timestep
  FILE *openFilePointer;
           // pointer to 3D bob output file
  FILE *vizFilePointer;
           // pointer to viz file


  /////////////////////////////////////////////////////////////////////////////
  // setting up the problem to be modeled
  //
  // David K. Cheng, Field and Wave Electromagnetics, 2nd ed., 
  //     pages 554-555. 
  // Rectangular waveguide, interior width = 2.29cm, interior height = 1.02cm.
  // This is a WG-16 waveguide useful for X-band applications.
  //
  // Choosing nx, ny, and nz:
  // There should be at least 20 cells per wavelength in each direction, 
  //     but we'll go with 25 so the animation will look prettier.   
  //     (CELLS_PER_WAVELENGTH was set to 25.0 in the global 
  //     constants at the beginning of the code.)
  // The number of cells along the width of the guide and the width of 
  //     those cells should fit the guide width exactly, so that ny*dy 
  //     = GUIDE_WIDTH meters.  
  //     The same should be true for nz*dz = GUIDE_HEIGHT meters.  
  // dx is chosen to be dy or dz -- whichever is smaller
  // nx is chosen to make the guide LENGTH_IN_WAVELENGTHS 
  //     wavelengths long.  
  // 
  // dt is chosen for Courant stability; the timestep must be kept small 
  //     enough so that the plane wave only travels one cell length 
  //     (one dx) in a single timestep.  Otherwise FDTD cannot keep up 
  //     with the signal propagation, since FDTD computes a cell only from 
  //     it's immediate neighbors.  

  // wavelength in meters:
  lambda = LIGHT_SPEED/FREQUENCY; 
  // angular frequency in radians/second:
  omega = 2.0*M_PI*FREQUENCY; 

  // set ny and dy:
  // start with a small ny:
  ny = 3;  
  // calculate dy from the guide width and ny:
  dy = GUIDE_WIDTH/ny;
  // until dy is less than a twenty-fifth of a wavelength,
  //     increment ny and recalculate dy:
  while(dy >= lambda/CELLS_PER_WAVELENGTH)
    {
      ny++;
      dy = GUIDE_WIDTH/ny;
    }

  // set nz and dz:
  // start with a small nz:
  nz = 3;  
  // calculate dz from the guide height and nz:
  dz = GUIDE_HEIGHT/nz;
  // until dz is less than a twenty-fifth of a wavelength,
  //     increment nz and recalculate dz:
  while(dz >= lambda/CELLS_PER_WAVELENGTH)
    {
      nz++;
      dz = GUIDE_HEIGHT/nz;
    }

  // set dx, nx, and dt:
  // set dx equal to dy or dz, whichever is smaller:
  dx = (dy < dz) ? dy : dz;
  // choose nx to make the guide LENGTH_IN_WAVELENGTHS 
  //     wavelengths long:  
  nx = (int)(LENGTH_IN_WAVELENGTHS*lambda/dx);
  // chose dt for Courant stability:
  dt = 1.0/(LIGHT_SPEED*sqrt(1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)));
  // time in seconds that will be simulated by the program:
  totalSimulatedTime = MAXIMUM_ITERATION*dt;

  // constants used in the field update equations:
  dtmudx = dt/(MU_0*dx);
  dtepsdx = dt/(EPSILON_0*dx);
  dtmudy = dt/(MU_0*dy);
  dtepsdy = dt/(EPSILON_0*dy);
  dtmudz = dt/(MU_0*dz);
  dtepsdz = dt/(EPSILON_0*dz);
  

  /////////////////////////////////////////////////////////////////////////////
  // memory allocation for the FDTD mesh:
  // There is a separate array for each of the six vector components, 
  //      ex, ey, ez, hx, hy, and hz.
  // The mesh is set up so that tangential E vectors form the outer faces of 
  //     the simulation volume.  There are nx*ny*nz cells in the mesh, but 
  //     there are nx*(ny+1)*(nz+1) ex component vectors in the mesh.
  //     There are (nx+1)*ny*(nz+1) ey component vectors in the mesh.
  //     There are (nx+1)*(ny+1)*nz ez component vectors in the mesh.
  // If you draw out a 2-dimensional slice of the mesh, you'll see why
  //     this is.  For example if you have a 3x3x3 cell mesh, and you 
  //     draw the E field components on the z=0 face, you'll see that
  //     the face has 12 ex component vectors, 3 in the x-direction
  //     and 4 in the y-direction.  That face also has 12 ey components,
  //     4 in the x-direction and 3 in the y-direction.  


  // Allocate memory for the E field arrays:

  // compute offset constants used in macros:
  ioff_ex = (ny+1)*(nz+1);
  ioff_ey = ny*(nz+1);
  ioff_ez = (ny+1)*nz;
  // allocate arrays:
  ex = calloc(nx*(ny+1)*(nz+1),sizeof(double));
  ey = calloc((nx+1)*ny*(nz+1),sizeof(double));
  ez = calloc((nx+1)*(ny+1)*nz,sizeof(double));
  // check if allocation successful
  if (ex == NULL || ey == NULL || ez == NULL) {
    fprintf(stderr,"E-field allocation failed.  Terminating...\n");
    exit(-1);
  }
  // keep rough track of allocated memory:
  allocatedBytes += ( (nx)*(ny+1)*(nz+1) * sizeof(double));
  allocatedBytes += ( (nx+1)*(ny)*(nz+1) * sizeof(double));
  allocatedBytes += ( (nx+1)*(ny+1)*(nz) * sizeof(double));

  // Allocate the H field arrays:
  // Since the H arrays are staggered half a step off 
  //     from the E arrays in every direction, the H 
  //     arrays are one cell smaller in the x, y, and z 
  //     directions than the corresponding E arrays. 
  // By this arrangement, the outer faces of the mesh
  //     consist of E components only, and the H 
  //     components lie only in the interior of the mesh.  

  // compute offset constants used in macros:
  ioff_hx = ny*nz;
  ioff_hy = (ny-1)*nz;
  ioff_hz = ny*(nz-1);
  // allocate arrays:
  hx = calloc((nx-1)*ny*nz,sizeof(double));
  hy = calloc(nx*(ny-1)*nz,sizeof(double));
  hz = calloc(nx*ny*(nz-1),sizeof(double));
  // check if allocation successful 
  if (hx == NULL || hy == NULL || hz == NULL) {
    fprintf(stderr,"E-field allocation failed.  Terminating...\n");
    exit(-1);
  }
  // keep rough track of allocated memory:
  allocatedBytes += ( (nx-1)*(ny)*(nz) * sizeof(double));
  allocatedBytes += ( (nx)*(ny-1)*(nz) * sizeof(double));
  allocatedBytes += ( (nx)*(ny)*(nz-1) * sizeof(double));


  /////////////////////////////////////////////////////////////////////////////
  // write some progress notes to standard output

  // print out some identifying information 
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "ToyFDTD2 version 1.0\n");
  fprintf(stdout, "Copyright (C) 1998,1999 Laurie E. Miller, Paul Hayes, ");
  fprintf(stdout, "Matthew O'Keefe\n");
  fprintf(stdout, "Copyright (C) 1999 John Schneider\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "ToyFDTD2 is free software published under the terms\n"); 
  fprintf(stdout, "of the GNU General Public License as published by the\n"); 
  fprintf(stdout, "Free Software Foundation.\n");  
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  // print out a bob command line,  
  //     including the dimensions of the output files: 
  fprintf(stdout, "bob -cmap chengGbry.cmap -s %dx%dx%d *.bob\n", nx+1, ny+1, nz); 
  fprintf(stdout, "\n"); 
  // print out a viz command line: 
  fprintf(stdout, "viz ToyFDTD2c.viz\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  // print out how much memory has been allocated:  
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes); 
  fprintf(stdout, "\n"); 
  // print out some simulation parameters: 
  fprintf(stdout, "PLOT_MODULUS = %d\n", PLOT_MODULUS); 
  fprintf(stdout, "rectangular waveguide\n"); 
  fprintf(stdout, "Stimulus = %lg Hertz ", FREQUENCY); 
  fprintf(stdout, "continuous plane wave\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Meshing parameters:\n"); 
  fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz); 
  fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz); 
  fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",  
	  GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Time simulated will be %lg seconds, %d timesteps\n",  
	  totalSimulatedTime, MAXIMUM_ITERATION); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "3D output scaling parameters:\n"); 
  fprintf(stdout, "Autoscaling every timestep\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  // print out some info on each timestep: 
  fprintf(stdout, "Following is the iteration number and current\n"); 
  fprintf(stdout, "simulated time for each timestep/iteration of\n"); 
  fprintf(stdout, "the simulation.  For each timestep that 3D data is\n"); 
  fprintf(stdout, "output to file, the maximum and minimum data\n"); 
  fprintf(stdout, "values are printed here with the maximum and\n"); 
  fprintf(stdout, "minimum scaled values in parentheses.\n"); 
  fprintf(stdout, "\n"); 

  /////////////////////////////////////////////////////////////////////////////
  // open and start writing the .viz file
  //    this file will be handy to feed parameters to viz if desired

  while ((vizFilePointer = fopen("ToyFDTD2c.viz", "w")) == NULL)
    {
      fprintf(stderr, "Difficulty opening ToyFDTD2c.viz");
      perror(" ");
    }
  fprintf(vizFilePointer, "#Viz V1.0\n");
  fprintf(vizFilePointer, "time: %lg %lg\n", currentSimulatedTime, dt);
  fprintf(vizFilePointer, "color: chengGbry.cmap\n");
  fprintf(vizFilePointer, "\n");


  /////////////////////////////////////////////////////////////////////////////
  // main loop:

  for(iteration = 0; iteration < MAXIMUM_ITERATION; iteration++)
    {// mainloop

      ///////////////////////////////////////////////////////////////////////// 
      // Output section: 

      // time in simulated seconds that the simulation has progressed:
      currentSimulatedTime = dt*(double)iteration;  
      // print to standard output the iteration number 
      //     and current simulated time:
      fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);

      // 3D data output every PLOT_MODULUS timesteps:
      //     The first time through the main loop all the data written to 
      //     file will be zeros.  If anything is nonzero, there's a bug.  :>
 
      if ( (iteration % PLOT_MODULUS) == 0)
	{// bob output section

	  // create the filename for this iteration, 
	  //     which includes the iteration number:
	  sprintf(filename, "c_%06d.bob", iteration);
	  // open a new data file for this iteration:
	  while ((openFilePointer = fopen(filename, "wb")) == NULL)
	    {// if the file fails to open, print an error message to 
	     //     standard output:
	      fprintf(stderr, "Difficulty opening c_%06d.bob", iteration);
	      perror(" ");
	    }

	  // find the max and min values to be output this timestep:  
	  min = FLT_MAX;
	  max = -FLT_MAX;
	  for(i=0; i<(nx+1); i++)
	    { 
	      for(j=0; j<(ny+1); j++)
		{
		  for(k=0; k<(nz); k++)
		    {
                      if (Ez(i,j,k) < min) min = Ez(i,j,k);
                      if (Ez(i,j,k) > max) max = Ez(i,j,k);
		    }
		}
	    }

	  // update the tracking variables for minimum and maximum 
	  //      values for the entire simulation:
	  if (min < simulationMin) simulationMin = min;
	  if (max > simulationMax) simulationMax = max;

	  // set norm to be max or min, whichever is greater in magnitude:
	  norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
	  if (norm == 0.0)
	    {// if everything is zero, give norm a tiny value 
	     //     to avoid division by zero:
	      norm = DBL_EPSILON;
	    }
	  scalingValue = 127.0/norm;
	  // write to standard output the minimum and maximum values 
	  //     from this iteration and the minimum and maximum values
	  //     that will be written to the bob file this iteration:
	  fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
		  min, (int)(127.0 + scalingValue*min),
		  max, (int)(127.0 + scalingValue*max));

	  // scale each ez value in the mesh to the range of integers 
	  //     from zero through 254 and write them to the output 
	  //     file for this iteration:
	  for(k=0; k<(nz); k++)
	    { 
	      for(j=0; j<(ny+1); j++)
		{
		  for(i=0; i<(nx+1); i++)
		    {
		      putc((int)(127.0 + 
				 scalingValue*Ez(i,j,k)), openFilePointer);
		    }
		}
	    }
	  
	  // close the output file for this iteration:
	  fclose(openFilePointer);

	  // write the dimensions and name of the output file for this 
	  //     iteration to the viz file:
	  fprintf(vizFilePointer, "%dx%dx%d %s\n", nx+1, ny+1, nz, filename);
	  // write identification of the corners of the mesh and the max
	  //     and min values for this iteration to the viz file:
	  fprintf(vizFilePointer, "bbox: 0.0 0.0 0.0 %lg %lg %lg %lg %lg\n",
		  dx*(double)nx, dy*(double)ny, dz*(double)nz, min, max);

	}// end bob output section

      /////////////////////////////////////////////////////////////////////////
      // Compute the stimulus: a plane wave emanates from the x=0 face:
      //     The length of the guide lies in the x-direction, the width of the 
      //     guide lies in the y-direction, and the height of the guide lies
      //     in the z-direction.  So the guide is sourced by all the ez 
      //     components on the stimulus face.  

      stimulus = sin(omega*currentSimulatedTime);
      for (i=0; i<(1); i++)
	{ 
	  for(j=0; j<(ny+1); j++)
	    {
	      for(k=0; k<nz; k++)
		{
                  Ez(i,j,k) = stimulus;
		}
	    }
	}

      /////////////////////////////////////////////////////////////////////////
      // Update the interior of the mesh:
      //    all vector components except those on the faces of the mesh.  
      //
      // Update all the H field vector components within the mesh:
      //     Since all H vectors are internal, all H values are updated here.
      //     Note that the normal H vectors on the faces of the mesh are not 
      //     computed here, and in fact were never allocated -- the normal 
      //     H components on the faces of the mesh are never used to update 
      //     any other value, so they are left out of the memory allocation  
      //     entirely.   
      
      // Update the hx values:
      for(i=0; i<(nx-1); i++)
	{  
	  for(j=0; j<(ny); j++)
	    {
	      for(k=0; k<(nz); k++)
		{                  
		  Hx(i,j,k) += (dtmudz*(Ey(i+1,j,k+1) - Ey(i+1,j,k)) - 
				dtmudy*(Ez(i+1,j+1,k) - Ez(i+1,j,k)));
		}
	    }
	}

      // Update the hy values:
      for(i=0; i<(nx); i++)
	{  
	  for(j=0; j<(ny-1); j++)
	    {
	      for(k=0; k<(nz); k++)
		{
                  Hy(i,j,k) += (dtmudx*(Ez(i+1,j+1,k) - Ez(i,j+1,k)) -
				dtmudz*(Ex(i,j+1,k+1) - Ex(i,j+1,k)));
		}
	    }
	}

      // Update the hz values:
      for(i=0; i<(nx); i++)
	{  
	  for(j=0; j<(ny); j++)
	    {
	      for(k=0; k<(nz-1); k++)
		{
                  Hz(i,j,k) += (dtmudy*(Ex(i,j+1,k+1) - Ex(i,j,k+1)) -
				dtmudx*(Ey(i+1,j,k+1) - Ey(i,j,k+1)));
		}
	    }
	}

      // Update the E field vector components.  
      // The values on the faces of the mesh are not updated here; they 
      //      are handled by the boundary condition computation 
      //      (and stimulus computation).  

      // Update the ex values:
      for(i=0; i<(nx); i++)
	{  
	  for(j=1; j<(ny); j++)
	    {
	      for(k=1; k<(nz); k++)
		{
                  Ex(i,j,k) += (dtepsdy*(Hz(i,j,k-1) - Hz(i,j-1,k-1)) -
				dtepsdz*(Hy(i,j-1,k) - Hy(i,j-1,k-1)));
		}
	    }
	}      

      // Update the ey values:
      for(i=1; i<(nx); i++)
	{  
	  for(j=0; j<(ny); j++)
	    {
	      for(k=1; k<(nz); k++)
		{
                  Ey(i,j,k) += (dtepsdz*(Hx(i-1,j,k) - Hx(i-1,j,k-1)) -
				dtepsdx*(Hz(i,j,k-1) - Hz(i-1,j,k-1)));
		}
	    }
	}

      // Update the ez values:
      for(i=1; i<(nx); i++)
	{  
	  for(j=1; j<(ny); j++)
	    {
	      for(k=0; k<(nz); k++)
		{
		  Ez(i,j,k) += (dtepsdx*(Hy(i,j-1,k) - Hy(i-1,j-1,k)) -
				dtepsdy*(Hx(i-1,j,k) - Hx(i-1,j-1,k)));
		}
	    }
	}
      
      fprintf(stdout, "\n");

      /////////////////////////////////////////////////////////////////////////
      // Compute the boundary conditions:

      // OK, so I'm yanking your chain on this one.  The PEC condition is 
      // enforced by setting the tangential E field components on all the
      // faces of the mesh to zero every timestep (except the stimulus 
      // face).  But the lazy/efficient way out is to initialize those 
      // vectors to zero and never compute them again, which is exactly 
      // what happens in this code.  
      
    }// end mainloop

  ///////////////////////////////////////////////////////////////////// 
  // Output section:  
  // The output routine is repeated one last time to write out  
  // the last data computed.   
  
  // time in simulated seconds that the simulation has progressed: 
  currentSimulatedTime = dt*(double)iteration;                           
  // print to standard output the iteration number 
  //     and current simulated time:
  fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);

  // 3D data output for the last timestep: 
  // create the filename for this iteration, 
  //     which includes the iteration number:
  sprintf(filename, "c_%06d.bob", iteration);
  // open a new data file for this iteration:
  while ((openFilePointer = fopen(filename, "wb")) == NULL)
    {// if the file fails to open, print an error message to 
      //     standard output:
      fprintf(stderr, "Difficulty opening c_%06d.bob", iteration);
      perror(" ");
    }
  
  // find the max and min values to be output this timestep:  
  min = FLT_MAX;
  max = -FLT_MAX;
  for(i=0; i<(nx+1); i++)
    { 
      for(j=0; j<(ny+1); j++)
        {
	  for(k=0; k<(nz); k++)
	    {
	      if (Ez(i,j,k) < min) min = Ez(i,j,k);
	      if (Ez(i,j,k) > max) max = Ez(i,j,k);
	    }
        }
    }
  
  // update the tracking variables for minimum and maximum 
  //      values for the entire simulation:
  if (min < simulationMin) simulationMin = min;
  if (max > simulationMax) simulationMax = max;
  
  // set norm to be max or min, whichever is greater in magnitude:	  
  norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
  if (norm == 0.0)
    {// if everything is zero, give norm a tiny value 
      //     to avoid division by zero:
      norm = DBL_EPSILON;
    }
  scalingValue = 127.0/norm;
  // write to standard output the minimum and maximum values 
  //     from this iteration and the minimum and maximum values
  //     that will be written to the bob file this iteration:
  fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
  	  min, (int)(127.0 + scalingValue*min),
  	  max, (int)(127.0 + scalingValue*max));
  
  // scale each ez value in the mesh to the range of integers 
  //     from zero through 254 and write them to the output 
  //     file for this iteration:
  for(k=0; k<(nz); k++)
    { 
      for(j=0; j<(ny+1); j++)
        {
	  for(i=0; i<(nx+1); i++)
	    {
	      putc((int)(127.0 + 
			 scalingValue*Ez(i,j,k)), openFilePointer);
	    }
        }
    }
  
  // close the output file for this iteration:
  fclose(openFilePointer);
  
  // write the dimensions and name of the output file for this 
  //     iteration to the viz file:
  fprintf(vizFilePointer, "%dx%dx%d %s\n", nx+1, ny+1, nz, filename);
  // write identification of the corners of the mesh and the max
  //     and min values for this iteration to the viz file:
  fprintf(vizFilePointer, "bbox: 0.0 0.0 0.0 %lg %lg %lg %lg %lg\n",
	  dx*(double)nx, dy*(double)ny, dz*(double)nz, min, max);

  /////////////////////////////////////////////////////////////////////////////
  // close the viz file for this simulation:
  fclose(vizFilePointer);

  /////////////////////////////////////////////////////////////////////////////
  // write some progress notes to standard output:

  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  // print out how much memory has been allocated:  
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes); 
  fprintf(stdout, "\n"); 
  // print out some simulation parameters: 
  fprintf(stdout, "PLOT_MODULUS = %d\n", PLOT_MODULUS); 
  fprintf(stdout, "rectangular waveguide\n"); 
  fprintf(stdout, "Stimulus = %lg Hertz ", FREQUENCY); 
  fprintf(stdout, "continuous plane wave\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Meshing parameters:\n"); 
  fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz); 
  fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz); 
  fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",  
	  GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "Time simulated was %lg seconds, %d timesteps\n",  
	  totalSimulatedTime, MAXIMUM_ITERATION); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "3D output scaling parameters:\n"); 
  fprintf(stdout, "Autoscaling every timestep\n"); 
  fprintf(stdout, "\n"); 
  // print out simulationMin and simulationMax: 
  fprintf(stdout, "Minimum output value was %lg \n", simulationMin); 
  fprintf(stdout, "Maximum output value was %lg \n", simulationMax); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  // print out a bob command line, including the dimensions  
  //      of the output files: 
  fprintf(stdout, "bob -cmap chengGbry.cmap -s %dx%dx%d *.bob\n",  
	  nx+1, ny+1, nz); 
  fprintf(stdout, "\n"); 
  // print out a viz command line: 
  fprintf(stdout, "viz ToyFDTD2c.viz\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 
  fprintf(stdout, "\n"); 

}// end main	

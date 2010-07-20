// ToyBoxFDTDbezhig, version 1.0
//      An if-I-can-do-it-you-can-do-it FDTD! 
// Copyright (C) 1998, 1999 Laurie E. Miller, Paul Hayes, Matthew O'Keefe 

// This program is free software; you can redistribute it and/or 
//     modify it under the terms of the GNU General Public License 
//     as published by the Free Software Foundation; either version 2
//     of the License, or any later version, with the following conditions
//     attached in addition to any and all conditions of the GNU
//     General Public License:
//     When reporting or displaying any results or animations created
//     using this code or modification of this code, make the appropriate
//     citation referencing ToyBoxFDTDbezhig by name and including the 
//     version number.  
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

// Contacting the authors:
//
// Laurie E. Miller, Paul Hayes, Matthew O'Keefe
// Department of Electrical and Computer Engineering
//      200 Union Street S. E.
//      Minneapolis, MN 55455
//
// lemiller@lcse.umn.edu
// info@cemtach.org
// 
// http://cemtach.org
// http://cemtach.org/software/ToyBoxFDTD/ToyBoxFDTD.html

///////////////////////////////////////////////////////////////////////////////
// This ToyBoxFDTDbezhig is an extremely simple 3D FDTD code. It demonstrates 
//      4 simple features by adding them to the original ToyFDTD 1.0 code:
//
//      1. Perfect Magnetic Conductor (PMC) boundary condition
//      2. Simplified sinusoidal pulse stimulus, rather than continuous 
//         plane wave
//      3. Output to a file tracking a field component at a single point 
//         in the mesh
//      4. No autoscaling of the output files each timestep: the files for 
//         each timestep can all be given uniform scaling by setting two
//         global constants to the maximum and minimum data values the
//         simulation will output
//         -- Note: this only works well if you supply appropriate numbers
//         for the maximum and minimum values.  If you're not sure what they
//         should be, stick with autoscaling.
//
// The original ToyFDTD is a stripped-down, minimalist 3D FDTD in C 
//      demonstrating the minimum considerations necessary in making a
//      simple 3D FDTD simulation. 
///////////////////////////////////////////////////////////////////////////////

// This is a very simple Yee algorithm 3D FDTD code in C implementing
//      the free space form of Maxwell's equations on a Cartesian grid.
//      This code is here for everyone, but not everyone will need something 
//      so simple, and not everyone will need to read all the comments.  
// There are no internal materials or geometry.  
// The code as delivered simulates an idealized parallel plate waveguide by 
//      treating the interior of the mesh as free space/air and enforcing PEC
//      (Perfect Electric Conductor) conditions on two faces of the mesh, and
//      PMC (Perfect Magnetic Conductor) conditions on two faces.  The
//      remaining two faces of the mesh form the ends of the waveguide; both
//      are PECs except at the beginning of the simulation when one emits a
//      sinusoidal pulse.
//          If the #define PMC statement near the beginning of the code is 
//      commented out, the sidewalls of the guide revert to PEC, and the
//      guide becomes a rectangular waveguide.
//          If the #define SINUSOIDAL_PULSE_STIMULUS line is commented out,
//      the input stimulus reverts to continuous-wave rather than the
//      pulse.
// The problem is modified from Field and Wave Electromagnetics, 2nd ed., by 
//      David K. Cheng, pages 554-555.  The original problem (simulated by 
//      ToyFDTD) was a WG-16 waveguide useful for X-band applications, 
//      interior width = 2.29cm, interior height = 1.02cm.  The frequency
//      (10 GHz) is chosen to be in the middle of the frequency range for
//      TE10 operation of the original rectangular waveguide.  
// Boundaries: PEC (Perfect Electric Conductor) and PMC (Perfect Magnetic 
//      Conductor).
// Stimulus: A simplified sinusoidal pulse emanates from the x = 0 face.
// 3D Output: The electric field intensity vector components in the
//      direction of the height of the guide (ez) are output to file every 
//      PLOT_MODULO timesteps, scaled to the range of integers from zero 
//      through 255 and with the integer value 127 held to be equal to data 
//      zero for every timestep.  Each timestep has it's own 3D output file.
//      This data output file format can be used in several visualization
//      tools, such as animabob and viz. There are several scaling options,
//      determined by #define statements:
//          1. Scaling is performed on each output timestep individually, 
//      since it is not known in advance what the maximum and minimum 
//      values will be for the entire simulation.  This method of auto-scaling
//      every output timestep can be very helpful in a simulation where the
//      intensities are sometimes strong and sometimes faint, since it will
//      highlight the  presence and structure of faint signals when stronger
//      signals have  left the mesh.  To use this option, comment out the
//      lines #define NO_AUTOSCALING, #define IDIOTPROOF_SCALING, and
//      #define NO_MIN_MAX_CALC.
//          2. Scaling uses the same scaling constant for the whole simulation.
//      This makes the resulting animation smoother in how it handles colors,
//      but less detail will show in timesteps where the data values are
//      comparatively small.  To use this option, uncomment the 
//      #define NO_AUTOSCALING line and set the global constants 
//      SCALING_MAX and SCALING_MIN to the maximum and minimum
//      values to be output for the entire simulation.  If you set these 
//      values too low in magnitude, you may get output values that
//      are outside the 0 to 255 range and various weirdnesses may result.
//      One way to find the correct values: run the simulation once, and look
//      at the end of the standard output.  There will be lines showing the
//      minimum and maximum data values sent to 3D output for the
//      simulation..  For the parallel-plate guide and the pulse source with
//      a PLOT_MODULUS of 5, the max/min values are about +-1.082, so
//      SCALING_MAX and SCALING_MIN are preset to these values.  
//          3. Scaling is global to all output timesteps as in option 2, but
//      a check is performed to make sure no output value falls outside the
//      0 to 255 range.  Scaled values greater than 255 are set to the
//      global constant OVERFLOW_DEFAULT and scaled values less than 0
//      are set to  UNDERFLOW_DEFAULT. To use this option, uncomment the 
//      #define NO_AUTOSCALING and #define IDIOTPROOF_SCALING lines,
//      and set set the global constants SCALING_MAX, SCALING_MIN, 
//      OVERFLOW_DEFAULT, and UNDERFLOW_DEFAULT.  
//          4. When autoscaling is turned off, as in option 2 and 3 above,
//      the code that tracks the maximum and minimum values for every
//      output timestep and for the whole simulation is no longer strictly
//      necessary.  It can be turned off by uncommenting the #define
//      NO_MIN_MAX_CALC.  Don't use this option if #define
//      NO_AUTOSCALING is uncommented; it will cause trouble. 
//
// Other Output: Notes on the progress of the simulation are written to 
//      standard output as the program runs.  
//
//      A .viz file is output to feed parameters to viz, should viz later be
//      used to view the data files.
//
//      The electric field intensity in the direction of the height of the guide
//      (ez) at one point near the center of the mesh is written out to a file
//      bezhigCenter.dat every timestep.  Viewing these data values as a 
//      function of time using a plotting tool like grace/xmgr or excel is 
//      fun because you can see the pulse bouncing back and forth through 
//      the mesh just by watching that one point.  First the positive pulse 
//      goes by, then it reflects back negative, then positive again...  You
//      can see the pulse start to get a little ragged over time.
//
//      The value of the stimulus function is also written out to a file, 
//      bezhigStimulus.dat, every timestep.  This is mainly for comparison with
//      the waveform in bezhigCenter.dat.  If the pulse initially seen by 
//      the center of the mesh doesn't look like the stimulus pulse, something
//      is wrong somewhere.

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
//     see the CEM TACH website at http://cemtach.org/software
//
// viz is another free visualization tool that displays and animates 
//     brick-of-byte files. For more information on viz or to download a copy, 
//     see the CEM TACH website at http://cemtach.org/software
//
// grace/xmgr is a free 2D plotting tool.  For more information on grace or
//     xmgr or to download a copy, see the CEM TACH website at
//     http://cemtach.org/software
//
// excel is an extremely non-free spreadsheet tool that can do plotting.
//     For more information on excel, see http://www.microsoft.com
//
// bezhig, as in ToyFDTDbezhig, is the Minnesota Ojibwe word for the
//     number one. If you want to know why it's called that, there's an
//     explanation of sorts on the web site somewhere.  Basically, it's 
//     about as descriptive without being misleading as anything else I
//     could think of.  


#include <math.h> 
#include <stdio.h>         
#include <float.h>

// constants that turn on/off the new features
#define PMC
          // If uncommented, this line includes the PMC part of the code at 
          //     compile time.  If this line is commented out, the open sides of
          //     the waveguide revert to PEC and the guide becomes rectangular
          //     rather than parallel-plate.
#define SINUSOIDAL_PULSE_STIMULUS
          // If uncommented, this line includes the sinusoidal pulse as the 
          //     stimulus for the simulation.  If this line is commented out,
          //     the stimulus will be a simplified continuous plane wave. 
#define POINT_OUTPUT
          // If uncommented, this line includes the routines that output
          //     single point values to bezhigStimulus.dat and bezhigCenter.dat
          //     If this line is commented out, those files will not be 
          //     created.  
//#define NO_AUTOSCALING
          // If uncommented, this line deactivates the feature that autoscales
          //     the 3D bob output files. The values for every output timestep
          //     are then scaled by the parameters SCALING_MAX and 
          //     SCALING_MIN set below.
//#define IDIOTPROOF_SCALING
          // If uncommented, this line includes code that makes sure the output
          //     values after scaling fall within the 0 to 255 range.  Values 
          //     over 255 are set to OVERFLOW_DEFAULT and values below 0
          //     are set to UNDERFLOW_DEFAULT as set below. 
//#define NO_MIN_MAX_CALC
          // If uncommented, this line excludes the code that finds the minimum
          //     and maximum value output each output timestep.  Generally,
          //     this option would be commented out.  If you uncomment this line,
          //     but leave NO_AUTOSCALING commented out, you're in for trouble.  

// program control constants
#define MAXIMUM_ITERATION 1000
          // total number of timesteps to be computed
          // should be equal to some_integer*PLOT_MODULO
          //     to make sure the last data computed is output
#define PLOT_MODULO 5
          // the program will output a 3D data file every PLOT_MODULO
          //     timesteps
#define FREQUENCY 10.0e9
          // frequency of the stimulus in Hertz
#define GUIDE_WIDTH 0.0229
          // meters -- section of the width of the guide to be simulated
          //    -- the parallel plate waveguide is actually infinite in width
#define GUIDE_HEIGHT 0.0102
          // interior height of the guide in meters
#define LENGTH_IN_WAVELENGTHS 5.0
          // length of the waveguide in wavelengths of the stimulus


// output scaling control constants
//   -- only used when NO_AUTOSCALING is defined
//   -- OVERFLOW_DEFAULT and UNDERFLOW_DEFAULT only used when 
//             IDIOTPROOF_SCALING also defined
#define SCALING_MAX 1.082
          // maximum raw data output value--defines top end of scaling range
#define SCALING_MIN -1.082
          // minimum raw data output value--defines low end of scaling range
          //     The code actually uses whichever of SCALING_MAX or SCALING_MIN
          //     has the greater absolute value to determine the scaling
          //     constant.  So if SCALING_MAX != -SCALING_MIN, the value with 
          //     the lesser absolute value will be ignored.  
#define OVERFLOW_DEFAULT 255
          // what scaled output values are set to if greater than SCALING_MAX
#define UNDERFLOW_DEFAULT 0
          // what scaled output values are set to if less than SCALING_MIN


// physical constants
#define LIGHT_SPEED     299792458.0       
          // speed of light in a vacuum in meters/second
#define LIGHT_SPEED_SQUARED 89875517873681764.0        
          // m^2/s^2
#define MU_0            (M_PI*4.0e-7)         
          // permeability of free space in henry/meter
#define EPSILON_0       (1.0/(MU_0*LIGHT_SPEED_SQUARED))      
          // permittivity of free space in farad/meter


main()	
{// main

  /////////////////////////////////////////////////////////////////////////////
  // variable declarations
  int i,j,k;     
           // indices of the 3D array of cells
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

  double ***ex, ***ey, ***ez;  
           // pointers to the arrays of ex, ey, and ez values
  double ***hx, ***hy, ***hz;
           // pointers to the arrays of hx, hy, and hz values

  // bob output routine variables:
  char filename[1024];
           // filename variable for 3D bob files
  double simulationMin = FLT_MAX;
           // tracks minimum value output by the entire simulation
  double simulationMax = -FLT_MAX;
           // tracks maximum value output by the entire simulation
  double min = 0.0;  
           // tracks minimum value output in one timestep
  double max = 0.0;
           // tracks maximum value output in one timestep
  double nextOut;
           // temp variable to hold next data value to be written
  double norm;
           // norm is set to be max or min, whichever is 
           //     greater in magnitude:
  double scalingValue;
           // multiplier used in output scaling
  FILE *openFilePointer;
           // pointer to 3D bob file
  FILE *vizFilePointer;
           // pointer to viz file

#ifdef POINT_OUTPUT
  // other output file variables:
  FILE *centerFilePointer;
           // pointer to center data output file
  FILE *stimulusFilePointer;
           // pointer to stimulus data output file
#endif

  /////////////////////////////////////////////////////////////////////////////
  // setting up the problem to be modeled
  //
  // Modified from David K. Cheng, Field and Wave Electromagnetics, 2nd ed., 
  //     pages 554-555. 
  // Parallel plate waveguide, interior height = 1.02cm, infinite width of which
  //     2.29 cm is simulated.
  //
  // Choosing nx, ny, and nz:
  // There should be at least 20 cells per wavelength in each direction, 
  //     but we'll go with 25 so the animation will look prettier.    
  // The number of cells along the width of the guide and the width of 
  //    those cells should fit the simulated guide width exactly, so that ny*dy 
  //    = GUIDE_WIDTH meters.  
  //    The same should be true for nz*dz = GUIDE_HEIGHT meters.  
  // dx is chosen to be dy or dz -- whichever is smaller
  // nx is chosen to make the guide LENGTH_IN_WAVELENGTHS 
  //     wavelengths long.  
  // 
  // dt is chosen for Courant stability; the time step must be kept small 
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
  while(dy >= lambda/25.0)
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
  while(dz >= lambda/25.0)
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
  dt = 1.0/(LIGHT_SPEED*sqrt(1/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz)));
  // time in seconds that will be simulated by the program:
  totalSimulatedTime = MAXIMUM_ITERATION*dt;

  // constants used in the field update equations:
  dtmudx = dt/(MU_0*dx);
  dtepsdx = dt/(EPSILON_0*dx);
  dtmudy = dt/(MU_0*dy);
  dtepsdy = dt/(EPSILON_0*dy);
  dtmudz = dt/(MU_0*dz);
  dtepsdz = dt/(EPSILON_0*dz);
  

#ifdef NO_AUTOSCALING
  /////////////////////////////////////////////////////////////////////////////
  // set the output scaling value for the main output files:
  //   -- This value is the same for every timestep when no autoscaling is 
  //      used, and can be set outside the main loop.  If autoscaling is used,
  //      scalingValue must be set within the main loop every output timestep.
  min = SCALING_MIN;
  max = SCALING_MAX;
  norm = (fabs(max) > fabs(min)) ? 
    fabs(max) : fabs(min);
    if (norm == 0.0)
      {// if everything is zero, give norm a tiny value 
       //     to avoid division by zero:
	norm = DBL_EPSILON;
      }
  scalingValue = 127.0/norm;
#endif

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

  // allocate the array of ex components:
  ex = (double ***)malloc((nx)*sizeof(double **));
  for(i=0; i<(nx); i++)
    {  
      ex[i] = (double **)malloc((ny+1)*sizeof(double *));
      for(j=0; j<(ny+1); j++)
	{
	  ex[i][j] = (double *)malloc((nz+1)*sizeof(double));
	  for(k=0; k<(nz+1); k++)
	    {
	      ex[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx)*(ny+1)*(nz+1) * sizeof(double));
  
  // allocate the array of ey components:
  ey = (double ***)malloc((nx+1)*sizeof(double **));
  for(i=0; i<(nx+1); i++)
    {  
      ey[i] = (double **)malloc((ny)*sizeof(double *));
      for(j=0; j<(ny); j++)
	{
	  ey[i][j] = (double *)malloc((nz+1)*sizeof(double));
	  for(k=0; k<(nz+1); k++)
	    {
	      ey[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx+1)*(ny)*(nz+1) * sizeof(double));

  // allocate the array of ez components:
  ez = (double ***)malloc((nx+1)*sizeof(double **));
  for(i=0; i<(nx+1); i++)
    {  
      ez[i] = (double **)malloc((ny+1)*sizeof(double *));
      for(j=0; j<(ny+1); j++)
	{
	  ez[i][j] = (double *)malloc((nz)*sizeof(double));
	  for(k=0; k<(nz); k++)
	    {
	      ez[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx+1)*(ny+1)*(nz) * sizeof(double));

  // Allocate the H field arrays:
  // Since the H arrays are staggered half a step off 
  //     from the E arrays in every direction, the H 
  //     arrays are one cell smaller in the x, y, and z 
  //     directions than the corresponding E arrays. 
  // By this arrangement, the outer faces of the mesh
  //     consist of E components only, and the H 
  //     components lie only in the interior of the mesh.  

  // allocate the array of hx components:
  hx = (double ***)malloc((nx-1)*sizeof(double **));
  for(i=0; i<(nx-1); i++)
    {  
      hx[i] = (double **)malloc((ny)*sizeof(double *));
      for(j=0; j<(ny); j++)
	{
	  hx[i][j] = (double *)malloc((nz)*sizeof(double));
	  for(k=0; k<(nz); k++)
	    {
	      hx[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx-1)*(ny)*(nz) * sizeof(double));
  
  // allocate the array of hy components:
  hy = (double ***)malloc((nx)*sizeof(double **));
  for(i=0; i<(nx); i++)
    {  
      hy[i] = (double **)malloc((ny-1)*sizeof(double *));
      for(j=0; j<(ny-1); j++)
	{
	  hy[i][j] = (double *)malloc((nz)*sizeof(double));
	  for(k=0; k<(nz); k++)
	    {
	      hy[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx)*(ny-1)*(nz) * sizeof(double));
  
  // allocate the array of hz components:
  hz = (double ***)malloc((nx)*sizeof(double **));
  for(i=0; i<(nx); i++)
    {  
      hz[i] = (double **)malloc((ny)*sizeof(double *));
      for(j=0; j<(ny); j++)
	{
	  hz[i][j] = (double *)malloc((nz-1)*sizeof(double));
	  for(k=0; k<(nz-1); k++)
	    {
	      hz[i][j][k] = 0.0;
	    }
	}
    }
  allocatedBytes += ( (nx)*(ny)*(nz-1) * sizeof(double));

  /////////////////////////////////////////////////////////////////////////////
  // write some progress notes to standard output

  // print out some identifying information 
  fprintf(stdout, "\n\nToyBoxFDTDbezhig version 1.0\n");
  fprintf(stdout, "Copyright (C) 1998,1999 Laurie E. Miller, Paul Hayes, ");
  fprintf(stdout, "Matthew O'Keefe\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "ToyBoxFDTDbezhig is free software published under the\n"); 
  fprintf(stdout, "terms of the GNU General Public License as published by\n"); 
  fprintf(stdout, "the Free Software Foundation.\n");  
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  // print out a bob command line, 
  //     including the dimensions of the output files:
  fprintf(stdout, "bob -cmap chengGbry2.cmap -s %dx%dx%d *.bob\n", 
	  nx+1, ny+1, nz);
  fprintf(stdout, "\n");
  // print out a viz command line:
  fprintf(stdout, "viz ToyBoxFDTDbezhig.viz\n");
  fprintf(stdout, "\n");
  // print out how much memory has been allocated: 
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes);
  fprintf(stdout, "\n");
  // print out some simulation parameters:
  fprintf(stdout, "PLOT_MODULO = %d\n", PLOT_MODULO);
#ifdef PMC
  fprintf(stdout, "parallel plate waveguide\n");
#else
  fprintf(stdout, "rectangular waveguide\n");
#endif
  fprintf(stdout, "Stimulus = %lg Hz ", FREQUENCY);
#ifdef SINUSOIDAL_PULSE_STIMULUS
  fprintf(stdout, "sinusoidal pulse\n");
#else
  fprintf(stdout, "continuous plane wave\n");
#endif
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
#ifndef NO_AUTOSCALING
  fprintf(stdout, "Autoscaling every timestep\n");
#endif
#ifdef NO_AUTOSCALING
  fprintf(stdout, "Global scaling over all timesteps\n");
  fprintf(stdout, "SCALING_MAX = %lg\n", SCALING_MAX);
  fprintf(stdout, "SCALING_MIN = %lg\n", SCALING_MIN);
#ifdef IDIOTPROOF_SCALING
  fprintf(stdout, "IDIOTPROOF_SCALING: test for overflow and underflow:\n");
  fprintf(stdout, "OVERFLOW_DEFAULT = %d\n", OVERFLOW_DEFAULT);
  fprintf(stdout, "UNDERFLOW_DEFAULT = %d\n", UNDERFLOW_DEFAULT);
#endif
#endif
  fprintf(stdout, "\n");
#ifdef POINT_OUTPUT
  fprintf(stdout, "Stimulus value for every timestep written to file\n");
  fprintf(stdout, "Ez at center of mesh for every timestep written to file\n");
#endif
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
  // open and start writing the .viz file:
  //    this file will be handy to feed parameters to viz if desired
  while ((vizFilePointer = fopen("ToyBoxFDTDbezhig.viz", "w")) == NULL)
    {
      fprintf(stderr, "Difficulty opening ToyBoxFDTDbezhig.viz");
      perror(" ");
    }
  fprintf(vizFilePointer, "#Viz V1.0\n");
  fprintf(vizFilePointer, "time: %lg %lg\n", currentSimulatedTime, dt);
  fprintf(vizFilePointer, "color: chengGbry2.cmap\n");
  fprintf(vizFilePointer, "\n");

#ifdef POINT_OUTPUT
  /////////////////////////////////////////////////////////////////////////////
  // open the stimulus data and center data output files:
  centerFilePointer = fopen("bezhigCenter.dat", "w");
  stimulusFilePointer = fopen("bezhigStimulus.dat", "w");
#endif

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

#ifdef POINT_OUTPUT
      // write values to the center data file:
      fprintf(centerFilePointer, "%lg %lg\n", currentSimulatedTime, 
	      ez[nx/2][ny/2][nz/2]);
      fflush(centerFilePointer);
#endif

      // 3D data output every PLOT_MODULO timesteps:
      //     The first time through the main loop all the data written to 
      //     file will be zeros.  If anything is nonzero, there's a bug.  :>

      // write out a bob file every PLOT_MODULO timesteps:
      if ( (iteration % PLOT_MODULO) == 0)
	{// bob output section

	  // print to standard output the filename for this iteration, 
	  //     which includes the iteration number:
	  sprintf(filename, "eh%05d.bob", iteration);
	  // open a new data file for this iteration:
	  while ((openFilePointer = fopen(filename, "wb")) == NULL)
	    {// if the file fails to open, print an error message to 
	     //     standard output:
	      fprintf(stderr, "Difficulty opening eh%05d.bob", iteration);
	      perror(" ");
	    }

#ifndef NO_MIN_MAX_CALC
	  // find the max and min values to be output this timestep:  
	  min = FLT_MAX;
	  max = -FLT_MAX;
	  for(i=0; i<(nx+1); i++)
	    { 
	      for(j=0; j<(ny+1); j++)
		{
		  for(k=0; k<(nz); k++)
		    {
		      if (ez[i][j][k] < min) min = ez[i][j][k];
		      if (ez[i][j][k] > max) max = ez[i][j][k];
		    }
		}
	    }

	  // update the tracking variables for minimum and maximum 
	  //      values for the entire simulation:
	  if (min < simulationMin) simulationMin = min;
	  if (max > simulationMax) simulationMax = max;
#endif

#ifndef NO_AUTOSCALING
	  // set the output scaling value for this timestep:
	  // set norm to be max or min, whichever is greater in magnitude:
	  norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
	  if (norm == 0.0)
	    {// if everything is zero, give norm a tiny value 
	     //     to avoid division by zero:
	      norm = DBL_EPSILON;
	    }
	  scalingValue = 127.0/norm;
#endif
	  // write to standard output the minimum and maximum values 
	  //     from this iteration and the minimum and maximum values
	  //     that will be written to the bob file this iteration:
	  fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
		  min, (int)(127.0 + scalingValue*min),
		  max, (int)(127.0 + scalingValue*max));

	  // scale each ez value in the mesh to the range of integers 
	  //     from zero through 255 and write them to the output 
	  //     file for this iteration:
	  for(k=0; k<(nz); k++)
	    { 
	      for(j=0; j<(ny+1); j++)
		{
		  for(i=0; i<(nx+1); i++)
		    {
#ifdef NO_AUTOSCALING
#ifdef IDIOTPROOF_SCALING
		      nextOut = (int)(127.0 + scalingValue*ez[i][j][k]);
		      if (nextOut < 0) nextOut = UNDERFLOW_DEFAULT;
		      else if (nextOut > 255) nextOut = OVERFLOW_DEFAULT;
		      putc(nextOut, openFilePointer);
#else 
		      putc((int)(127.0 + 
				 scalingValue*ez[i][j][k]), openFilePointer);
#endif
#endif

#ifndef NO_AUTOSCALING
		      putc((int)(127.0 + 
				 scalingValue*ez[i][j][k]), openFilePointer);
#endif
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
      // Compute the stimulus: a sinusoidal pulse emanates from the x=0 face:
      //     The length of the guide lies in the x-direction, the width of the 
      //     guide lies in the y-direction, and the height of the guide lies
      //     in the z-direction.  So the guide is sourced by all the ez 
      //     components on the stimulus face.  
      // The pulse is one complete wavelength of a shifted sinusoid that 
      //     ranges from zero to one in intensity.  So the stimulus varies 
      //     sinusoidally from zero to one to zero again and then terminates
      //     --the x=0 face becomes a PEC thereafter.  
      //
      // compute the current stimulus value:
#ifdef SINUSOIDAL_PULSE_STIMULUS
      if (currentSimulatedTime <= 1.0/FREQUENCY)
	{
	  stimulus = .5*(1.0 + sin(omega*currentSimulatedTime - M_PI*.5));
	}
      else
	{
	  stimulus = 0.0;
	}
#else
      stimulus = sin(omega*currentSimulatedTime);
#endif
      // set all vectors on the x=0 face to the value of stimulus:
      for (i=0; i<(1); i++)
	{ 
	  for(j=0; j<(ny+1); j++)
	    {
	      for(k=0; k<nz; k++)
		{
		  ez[i][j][k] = stimulus;
		}
	    }
	}
#ifdef POINT_OUTPUT
      // write the current stimulus value to bezhigStimulus.dat:
      fprintf(stimulusFilePointer, "%lg %lg\n",
	      currentSimulatedTime, stimulus);
      fflush(stimulusFilePointer);
#endif
	
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
		  hx[i][j][k] += (dtmudz*(ey[i+1][j][k+1] - ey[i+1][j][k]) - 
		                  dtmudy*(ez[i+1][j+1][k] - ez[i+1][j][k]));
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
		  hy[i][j][k] +=  (dtmudx*(ez[i+1][j+1][k] - ez[i][j+1][k]) -
      		                   dtmudz*(ex[i][j+1][k+1] - ex[i][j+1][k]));
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
		  hz[i][j][k] +=  (dtmudy*(ex[i][j+1][k+1] - ex[i][j][k+1]) -
      		                   dtmudx*(ey[i+1][j][k+1] - ey[i][j][k+1]));
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
		  ex[i][j][k] += (dtepsdy*(hz[i][j][k-1] - hz[i][j-1][k-1]) -
				  dtepsdz*(hy[i][j-1][k] - hy[i][j-1][k-1]));
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
		  ey[i][j][k] += (dtepsdz*(hx[i-1][j][k] - hx[i-1][j][k-1]) -
				  dtepsdx*(hz[i][j][k-1] - hz[i-1][j][k-1]));
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
		  ez[i][j][k] += (dtepsdx*(hy[i][j-1][k] - hy[i-1][j-1][k]) -
				  dtepsdy*(hx[i-1][j][k] - hx[i-1][j-1][k]));
		}
	    }
	}
      
      fprintf(stdout, "\n");

      /////////////////////////////////////////////////////////////////////////
      // Compute the boundary conditions:

#ifdef PMC
      // Compute the PMC symmetry boundaries:
      // j = 0 face, j = ny face
      for(i=0; i<nx; i++)
	{	
	  for(k=0; k<(nz+1); k++)
	    {
	      ex[i][0][k] = ex[i][ny-1][k];
	      ex[i][ny][k] = ex[i][1][k];
	    }
	}
      for(i=0; i<(nx+1); i++)
	{	
	  for(k=0; k<nz; k++)
	    {
	      ez[i][0][k] = ez[i][ny-1][k];
	      ez[i][ny][k] = ez[i][1][k];
	    }
	}
#endif

      // Compute the PEC boundaries:
      //
      // OK, so I'm yanking your chain on this one.  The PEC condition is 
      // enforced by setting the tangential E field components on the PEC
      // faces of the mesh to zero every timestep (except the stimulus 
      // face).  But the lazy/efficient way out is to initialize those 
      // vectors to zero and never compute them again, which is exactly 
      // what happens in this code.  

    }// end mainloop
  
  /////////////////////////////////////////////////////////////////////////
  // Output section: 
  // The output routine is repeated one last time to write out 
  // the last data computed.  
  
  // time in simulated seconds that the simulation has progressed:
  currentSimulatedTime = dt*(double)iteration;  
  // print to standard output the iteration number 
  //     and current simulated time:
  fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);
  
#ifdef POINT_OUTPUT
  // write values to the center data file:
  fprintf(centerFilePointer, "%lg %lg\n", currentSimulatedTime, 
	  ez[nx/2][ny/2][nz/2]);
  fflush(centerFilePointer);
#endif
  
  // 3D data output every PLOT_MODULO timesteps:
  
  // write out a bob file every PLOT_MODULO timesteps:
  if ( (iteration % PLOT_MODULO) == 0)
    {// bob output section
      
      // print to standard output the filename for this iteration, 
      //     which includes the iteration number:
      sprintf(filename, "eh%05d.bob", iteration);
      // open a new data file for this iteration:
      while ((openFilePointer = fopen(filename, "wb")) == NULL)
	{// if the file fails to open, print an error message to 
	  //     standard output:
	      fprintf(stderr, "Difficulty opening eh%05d.bob", iteration);
	      perror(" ");
	}
      
#ifndef NO_MIN_MAX_CALC
      // find the max and min values to be output this timestep:  
      min = FLT_MAX;
      max = -FLT_MAX;
      for(i=0; i<(nx+1); i++)
	{ 
	  for(j=0; j<(ny+1); j++)
	    {
	      for(k=0; k<(nz); k++)
		{
		  if (ez[i][j][k] < min) min = ez[i][j][k];
		  if (ez[i][j][k] > max) max = ez[i][j][k];
		}
	    }
	}

      // update the tracking variables for minimum and maximum 
      //      values for the entire simulation:
      if (min < simulationMin) simulationMin = min;
      if (max > simulationMax) simulationMax = max;
#endif

#ifndef NO_AUTOSCALING
      // set the output scaling value for this timestep:
      // set norm to be max or min, whichever is greater in magnitude:
      norm = (fabs(max) > fabs(min)) ? fabs(max) : fabs(min);
      if (norm == 0.0)
	{// if everything is zero, give norm a tiny value 
	  //     to avoid division by zero:
	  norm = DBL_EPSILON;
	}
      scalingValue = 127.0/norm;
#endif
      // write to standard output the minimum and maximum values 
      //     from this iteration and the minimum and maximum values
      //     that will be written to the bob file this iteration:
      fprintf(stdout, "\t%lg(%d) < ez BoB < %lg(%d)",
	      min, (int)(127.0 + scalingValue*min),
	      max, (int)(127.0 + scalingValue*max));
      
      // scale each ez value in the mesh to the range of integers 
      //     from zero through 255 and write them to the output 
      //     file for this iteration:
      for(k=0; k<(nz); k++)
	{ 
	  for(j=0; j<(ny+1); j++)
	    {
	      for(i=0; i<(nx+1); i++)
		{
#ifdef NO_AUTOSCALING
#ifdef IDIOTPROOF_SCALING
		  nextOut = (int)(127.0 + scalingValue*ez[i][j][k]);
		  if (nextOut < 0) nextOut = UNDERFLOW_DEFAULT;
		  else if (nextOut > 255) nextOut = OVERFLOW_DEFAULT;
		  putc(nextOut, openFilePointer);
#else 
		  putc((int)(127.0 + 
			     scalingValue*ez[i][j][k]), openFilePointer);
#endif
#endif

#ifndef NO_AUTOSCALING
		  putc((int)(127.0 + 
			     scalingValue*ez[i][j][k]), openFilePointer);
#endif
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

  /////////////////////////////////////////////////////////////////////////////
  // close the various files for this simulation:
  fclose(vizFilePointer);
#ifdef POINT_OUTPUT
  fclose(stimulusFilePointer);
  fclose(centerFilePointer);
#endif

  /////////////////////////////////////////////////////////////////////////////
  // write some progress notes to standard output:

  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  // print out a bob command line, 
  //     including the dimensions of the output files:
  fprintf(stdout, "bob -cmap chengGbry2.cmap -s %dx%dx%d *.bob\n", 
	  nx+1, ny+1, nz);
  fprintf(stdout, "\n");
  // print out a viz command line:
  fprintf(stdout, "viz ToyBoxFDTDbezhig.viz\n");
  fprintf(stdout, "\n");
  // print out how much memory has been allocated: 
  fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes);
  fprintf(stdout, "\n");
  // print out some simulation parameters:
  fprintf(stdout, "PLOT_MODULO = %d\n", PLOT_MODULO);
#ifdef PMC
  fprintf(stdout, "parallel plate waveguide\n");
#else
  fprintf(stdout, "rectangular waveguide\n");
#endif
  fprintf(stdout, "Stimulus = %lg Hz ", FREQUENCY);
#ifdef SINUSOIDAL_PULSE_STIMULUS
  fprintf(stdout, "sinusoidal pulse\n");
#else
  fprintf(stdout, "continuous plane wave\n");
#endif
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
#ifndef NO_AUTOSCALING
  fprintf(stdout, "Autoscaling every timestep\n");
#endif
#ifdef NO_AUTOSCALING
  fprintf(stdout, "Global scaling over all timesteps\n");
  fprintf(stdout, "SCALING_MAX = %lg\n", SCALING_MAX);
  fprintf(stdout, "SCALING_MIN = %lg\n", SCALING_MIN);
#ifdef IDIOTPROOF_SCALING
  fprintf(stdout, "IDIOTPROOF_SCALING: test for overflow and underflow:\n");
  fprintf(stdout, "OVERFLOW_DEFAULT = %d\n", OVERFLOW_DEFAULT);
  fprintf(stdout, "UNDERFLOW_DEFAULT = %d\n", UNDERFLOW_DEFAULT);
#endif
#endif
  fprintf(stdout, "\n");
#ifdef POINT_OUTPUT
  fprintf(stdout, "Stimulus value for every timestep written to file\n");
  fprintf(stdout, "Ez at center of mesh for every timestep written to file\n");
#endif
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  // print out simulationMin and simulationMax:
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Minimum output value was %lg \n", simulationMin);
  fprintf(stdout, "Maximum output value was %lg \n", simulationMax);
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  // print out a bob command line, including the dimensions 
  //      of the output files:
  fprintf(stdout, "bob -cmap chengGbry2.cmap -s %dx%dx%d *.bob\n", 
	  nx+1, ny+1, nz);
  fprintf(stdout, "\n");
  // print out a viz command line:
  fprintf(stdout, "viz ToyBoxFDTDbezhig.viz\n");
  fprintf(stdout, "\n");

}// end main	

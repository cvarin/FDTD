!! ToyFDTD1, version 1.03 (F90)
!!      The if-I-can-do-it-you-can-do-it FDTD! 
!! Copyright (C) 1998,1999 Laurie E. Miller, Paul Hayes, Matthew O'Keefe 
!!               1999 Max Smirnoff, Matt Rundquist (F90 translation)
!! This program is free software; you can redistribute it and/or 
!!     modify it under the terms of the GNU General Public License 
!!     as published by the Free Software Foundation; either version 2
!!     of the License, or any later version, with the following conditions
!!     attached in addition to any and all conditions of the GNU
!!     General Public License:
!!     When reporting or displaying any results or animations created
!!     using this code or modification of this code, make the appropriate
!!     citation referencing ToyFDTD1 by name and including the version
!!     number.  
!!
!! This program is distributed in the hope that it will be useful,
!!     but WITHOUT ANY WARRANTY; without even the implied warranty 
!!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!!     See the GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!!     along with this program; if not, write to the Free Software
!!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  
!!     02111-1307  USA

!! Contacting the authors:
!!
!! Laurie E. Miller, Paul Hayes, Matthew O'Keefe,
!! Max Smirnoff, Matt Rundquist
!! Department of Electrical and Computer Engineering
!!      200 Union Street S. E.
!!      Minneapolis, MN 55455
!! 
!! lemiller@borg.umn.edu
!! 
!! http://www.borg.umn.edu/toyfdtd/
!! http://www.borg.umn.edu/toyfdtd/ToyFDTD1.html
!! http://www.toyfdtd.org/

!! This code is here for everyone, but not everyone will need something 
!!      so simple, and not everyone will need to read all the comments.  
!! This file is over 700 lines long, but less than 400 of that is actually
!!      code. 

!! This ToyFDTD1 is a stripped-down, minimalist 3D FDTD code.  It 
!!      illustrates the minimum factors that must be considered to 
!!      create a simple FDTD simulation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Changes to version 1.03 from version 1.02: Updating contact & web info.
!! Changes to version 1.02 from version 1.0:  There is no 1.0 version of the 
!!      FORTRAN code.  For a listing of the changes to the C code, see 
!!      the changelog file.  
!! For some notes on the F90 translation, see the README file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This is a very simple Yee algorithm 3D FDTD code in F90 implementing
!!      the free space form of Maxwell's equations on a Cartesian grid.
!! There are no internal materials or geometry.  
!! The code simulates an idealized rectangular waveguide by treating the 
!!      interior of the mesh as free space/air and enforcing PEC (Perfect 
!!      Electric Conductor) conditions on the faces of the mesh.
!! The problem is taken from Field and Wave Electromagnetics, 2nd ed., by 
!!      David K. Cheng, pages 554-555.  It is a WG-16 waveguide 
!!      useful for X-band applications, interior width = 2.29cm, 
!!      interior height = 1.02cm.  The frequency (10 GHz) is chosen to be 
!!      in the middle of the frequency range for TE10 operation.  
!! Boundaries: PEC (Perfect Electric Conductor).
!! Stimulus: A simplified sinusoidal plane wave emanates from x = 0 face.
!! 3D output: The electric field intensity vector components in the direction 
!!      of the height of the guide (ez) are output to file every 
!!      PLOT_MODULUS timesteps (and for the last timestep), scaled to 
!!      the range of integers from zero through 254.  (The colormap 
!!      included in the tar file assigns rgb values to the range zero through
!!      255.)   Scaling is performed on each timestep individually, since 
!!      it is not known in advance what the maximum and minimum 
!!      values will be for the entire simulation.  The integer value 127 is 
!!      held to be equal to the data zero for every timestep.  This method 
!!      of autoscaling every timestep can be very helpful in a simulation 
!!      where the intensities are sometimes strong and sometimes faint, 
!!      since it will highlight the presence and structure of faint signals 
!!      when stronger signals have left the mesh.  
!!      Each timestep has it's own output file.  This data output file format
!!      can be used in several visualization tools, such as animabob and viz. 
!! Other output: Notes on the progress of the simulation are written to standard
!!      output as the program runs.  
!!      A .viz file is output to feed parameters to viz, should viz later be
!!      used to view the data files.

!! Some terminology used here:
!!
!! This code implements a Cartesian mesh with space differentials 
!!     of dx, dy, dz.
!! This means that a point in the mesh has neighboring points dx meters 
!!     away in the direction of the positive and negative x-axis,
!!     and neighboring points dy meters away in the directions 
!!     of the +- y-axis, and neighboring points dz meters away 
!!     in the directions of the +- z-axis,
!! The mesh has nx cells in the x direction, ny in the y direction, 
!!     and nz in the z direction.
!! ex, ey, and ez refer to the arrays of electric field intensity vectors 
!!     -- for example, ex is a 3-dimensional array consisting of the 
!!     x component of the E field intensity vector for every point in the 
!!     mesh.  ex(i,j,k) refers to the x component of the E field intensity 
!!     vector at point (i,j,k).  
!! hx, hy, and hz refer to the arrays of magnetic field intensity vectors.
!!
!! dt is the time differential -- the length of each timestep in seconds.
!!
!! bob is a file format that stands for "brick of bytes", meaning a string 
!!     of bytes that can be interpreted as a 3-dimensional array of byte 
!!     values (integers from zero through 255).  animabob is a free 
!!     visualization tool that displays and animates a sequence of bricks 
!!     of bytes.  For more information on animabob or to download a copy, 
!!     see the ToyFDTD homepage at http://www.borg.umn.edu/toyfdtd/
!!
!! viz is another free visualization tool that displays and animates 
!!     brick-of-byte files. For more information on viz or to download a copy, 
!!     see the ToyFDTD website at http://www.borg.umn.edu/toyfdtd/


!! program control constants
#define MAXIMUM_ITERATION 1000
          !! total number of timesteps to be computed
#define PLOT_MODULUS 5
          !! the program will output 3D data every PLOT_MODULUS timesteps
          !!     except for the last iteration computed, which is always output.
          !!     So if MAXIMUM_ITERATION is not an integer multiple of 
          !!     PLOT_MODULUS, the last timestep output will come after 
          !!     a shorter interval than that separating previous outputs.  
#define FREQUENCY 10.0d9
          !! frequency of the stimulus in Hertz
#define GUIDE_WIDTH 0.0229d0
          !! meters
#define GUIDE_HEIGHT 0.0102d0
          !! meters
#define LENGTH_IN_WAVELENGTHS 5.0d0
          !! length of the waveguide in wavelengths of the stimulus wave
#define CELLS_PER_WAVELENGTH 25.0d0
          !! minimum number of grid cells per wavelength in the x, y, and
          !!     z directions

!! physical constants
#define LIGHT_SPEED     299792458.0d0
          !! speed of light in a vacuum in meters/second
#define LIGHT_SPEED_SQUARED 89875517873681764.0d0    
          !! m^2/s^2
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271d-6
          !! permeability of free space in henry/meter
#define EPSILON_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932d-12
          !! permittivity of free space in farad/meter

!!
!! These parameters are allowed in C through include files in the compilation.
!!
#define M_PI 3.14159265358979323846d0
          !! pi as defined in /usr/include/math.h on SGI IRIX 6.2
#define FLT_MAX	3.40282347d+38
          !! max float as defined in /usr/include/float.h on SGI IRIX 6.2
#define	DBL_EPSILON	2.2204460492503131d-16
          !! DBL_EPSILON as defined in /usr/include/float.h on SGI IRIX 6.2
#define SIZEOF_DOUBLE 8
          !! bytes in a double-precision floating point number


program ToyFDTD1
implicit none !! All undeclared variables will be reported as errors;
              !! this prevents some major problems.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! variable declarations
  integer :: i,j,k
           !! indices of the 3D array cells
  integer :: nx,ny,nz
           !! total number of cells along the x,y,z axes, respectively
  integer :: allocatedBytes = 0;          
           !! a counter to track number of bytes allocated

  integer :: iteration = 0
           !! counter to track how many timesteps have been computed
  real*8 :: stimulus = 0.0
           !! value of the stimulus at current timestep
  real*8 :: currentSimulatedTime = 0.0
           !! time in simulated seconds that the simulation has progressed
  real*8 :: totalSimulatedTime = 0.0
           !! time in seconds that will be simulated by the program
  real*8 :: omega
           !! angular frequency in radians/second
  real*8 :: lambda
           !! wavelength of the stimulus in meters
  real*8 :: dx,dy,dz;
           !! space differentials (or dimensions of a single cell) in meters
  real*8 :: dt;
           !! time differential (how much time between timesteps) in seconds
  real*8 :: dtmudx,dtepsdx
           !! physical constants used in the field update equations  
  real*8 :: dtmudy,dtepsdy
           !! physical constants used in the field update equations  
  real*8 :: dtmudz,dtepsdz
           !! physical constants used in the field update equations  

  real*8,target,allocatable,dimension(:,:,:) :: ex,ey,ez
           !! pointers to the arrays of ex, ey, and ez values
  real*8,target,allocatable,dimension(:,:,:) :: hx,hy,hz
           !! pointers to the arrays of hx, hy, and hz values

  !! bob output routine variables:
  character(LEN=1024) :: filename 
           !! filename variable for 3D bob files
  real*8 :: simulationMin = FLT_MAX
           !! tracks minimum value output by the entire simulation
  real*8 :: simulationMax = -FLT_MAX
           !! tracks maximum value output by the entire simulation
  real*8 :: min, max
           !! these track minimum and maximum values output in one timestep
  real*8 :: norm
           !! norm is set each iteration to be max or min, whichever is 
           !!     greater in magnitude
  real*8 :: scalingValue
           !! multiplier used in output scaling, calculated every timestep

  character(LEN=10) :: numbers = "0123456789" 
           !! this is a handy string with all the decimal digits

  integer :: temp
  integer :: ios, record
           !! these are just dummy variables used during file writes 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! setting up the problem to be modeled
  !!
  !! David K. Cheng, Field and Wave Electromagnetics, 2nd ed., 
  !!     pages 554-555. 
  !! Rectangular waveguide, interior width = 2.29cm, interior height = 1.02cm.
  !! This is a WG-16 waveguide useful for X-band applications.
  !!
  !! Choosing nx, ny, and nz:
  !! There should be at least 20 cells per wavelength in each direction, 
  !!     but we'll go with 25 so the animation will look prettier.    
  !!     (CELLS_PER_WAVELENGTH was set to 25.0 in the global 
  !!     constants at the beginning of the code.)
  !! The number of cells along the width of the guide and the width of 
  !!    those cells should fit the guide width exactly, so that ny*dy 
  !!    = GUIDE_WIDTH meters.  
  !!    The same should be true for nz*dz = GUIDE_HEIGHT meters.  
  !! dx is chosen to be dy or dz -- whichever is smaller
  !! nx is chosen to make the guide LENGTH_IN_WAVELENGTHS 
  !!     wavelengths long.  
  !! 
  !! dt is chosen for Courant stability; the timestep must be kept small 
  !!     enough so that the plane wave only travels one cell length 
  !!     (one dx) in a single timestep.  Otherwise FDTD cannot keep up 
  !!     with the signal propagation, since FDTD computes a cell only from 
  !!     it's immediate neighbors.  


  !! wavelength in meters:
  lambda = LIGHT_SPEED/FREQUENCY

  !! angular frequency in radians/second:
  omega = 2.0d0*M_PI*FREQUENCY

  !! set ny and dy:
  !! start with small ny:
  ny = 3
  !! calculate dy from the guide width and ny:
  dy = GUIDE_WIDTH/ny
  !! until dy is less than a twenty-fifth of a wavelength,
  !!      increment ny and recalculate dy:
  do while (dy >= lambda/CELLS_PER_WAVELENGTH)
    ny = ny + 1
    dy = GUIDE_WIDTH/ny
  enddo

  !! start nz and dz:
  !! start with small nz:
  nz = 3
  !! calculate dz from the guide height and nz:
  dz = GUIDE_HEIGHT/nz
  !! until dz is less than a twenty-fifth of a wavelength,
  !!      increment nz and recalculate dz:
  do while (dz >= lambda/CELLS_PER_WAVELENGTH)
    nz = nz + 1
    dz = GUIDE_HEIGHT/nz
  enddo

  !! set dx, nx, and dt:
  !! set dx equal to dy or dz, whichever is smaller:
  if(dy < dz) then
    dx = dy
  else
    dx = dz
  end if
  !! choose nx to make the guide LENGTH_IN_WAVELENGTHS 
  !!     wavelengths long:  
  nx = int(LENGTH_IN_WAVELENGTHS*lambda/dx)
  !!  chose dt for Courant stability:
  dt = 1.0d0/(LIGHT_SPEED*sqrt(1.0d0/(dx*dx) + 1.0d0/(dy*dy) + 1.0d0/(dz*dz)))
  !!  time in seconds that will be simulated by the program:
  totalSimulatedTime = MAXIMUM_ITERATION*dt

  !!  constants used in the field update equations:
  dtmudx = dt/(MU_0*dx)
  dtepsdx = dt/(EPSILON_0*dx)
  dtmudy = dt/(MU_0*dy)
  dtepsdy = dt/(EPSILON_0*dy)
  dtmudz = dt/(MU_0*dz)
  dtepsdz = dt/(EPSILON_0*dz)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                  
  !! memory allocation for the FDTD mesh:
  !! There is a separate array for each of the six vector components, 
  !!      ex, ey, ez, hx, hy, and hz.
  !! There are nx*ny*nz cells in the mesh, but 
  !!     there are nx*(ny+1)*(nz+1) ex component vectors in the mesh.
  !!     There are (nx+1)*ny*(nz+1) ey component vectors in the mesh.
  !!     There are (nx+1)*(ny+1)*nz ez component vectors in the mesh.
  !! If you draw out a 2-dimensional slice of the mesh, you'll see why
  !!     this is.  For example if you have a 3x3x3 cell mesh, and you 
  !!     draw the E field components on the z=0 face, you'll see that
  !!     the face has 12 ex component vectors, 3 in the x-direction
  !!     and 4 in the y-direction.  That face also has 12 ey components,
  !!     4 in the x-direction and 3 in the y-direction.  


  !! Allocate memory for the E field arrays:

  !! allocate the array of ex components:
  allocate(ex(0:nx-1,0:ny,0:nz))
  do i=0,(nx-1),1
     do j=0,ny,1
        do k=0,nz,1
           ex(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx)*(ny+1)*(nz+1) * SIZEOF_DOUBLE)

  !! allocate the array of ey components:
  allocate(ey(0:nx,0:ny-1,0:nz))
  do i=0,nx,1
     do j=0,(ny-1),1
        do k=0,nz,1
           ey(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx+1)*(ny)*(nz+1) * SIZEOF_DOUBLE)

  !! allocate the array of ez components:
  allocate(ez(0:nx,0:ny,0:nz-1))
  do i=0,nx,1
     do j=0,ny,1
        do k=0,(nz-1),1
           ez(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx+1)*(ny+1)*(nz) * SIZEOF_DOUBLE)

  !! Allocate the H field arrays:
  !! Since the H arrays are staggered half a step off 
  !!     from the E arrays in every direction, the H 
  !!     arrays are one cell smaller in the x, y, and z 
  !!     directions than the corresponding E arrays. 
  !! By this arrangement, the outer faces of the mesh
  !!     consist of E components only, and the H 
  !!     components lie only in the interior of the mesh.  

  !! allocate the array of hx components:
  allocate(hx(0:nx-2,0:ny-1,0:nz-1))
  do i=0,(nx-2),1
     do j=0,(ny-1),1
        do k=0,(nz-1),1
           hx(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx-1)*(ny)*(nz) * SIZEOF_DOUBLE)

  !! allocate the array of hy components:
  allocate(hy(0:nx-1,0:ny-2,0:nz-1))
  do i=0,(nx-1),1
     do j=0,(ny-2),1
        do k=0,(nz-1),1
           hy(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx)*(ny-1)*(nz) * SIZEOF_DOUBLE)

  !! allocate the array of hz components:
  allocate(hz(0:nx-1,0:ny-1,0:nz-2))
  do i=0,(nx-1),1
     do j=0,(ny-1),1
        do k=0,(nz-2),1
           hz(i,j,k) = 0.0d0
        end do
     end do
  end do
  allocatedBytes = allocatedBytes + ( (nx)*(ny)*(nz-1) * SIZEOF_DOUBLE)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! write some progress notes to standard output

  !! print out some identifying information 
  print*," "
  print*," "
  print*," "
  print*," "
  print*,"ToyFDTD1 version 1.03 (F90)"
  print*,"Copyright (C) 1998, 1999 Laurie E. Miller, Paul Hayes, "
  print*,"Matthew O'Keefe, Max Smirnoff, Matt Rundquist"
  print*," "
  print*,"ToyFDTD1 is free software published under the terms" 
  print*,"of the GNU General Public License as published by the" 
  print*,"Free Software Foundation. " 
  print*," "
  print*," "
  !! print out a bob command line, 
  !! including the dimensions of the output files:
  print'(" bob -cmap chengGbry.cmap -s ",I5.5,"x",I5.5,"x",I5.5," *.bob")',nx+1,ny+1,nz
  print*," "
  !! print out a viz command line:
  print*,"viz ToyFDTD1f90.viz"
  print*," "
  print*," "
  !! print out how much memory has been allocated: 
  print*,"Dynamically allocated ", allocatedBytes, " bytes"
  print*," "
  !! print out some simulation parameters: 
  print*,"PLOT_MODULUS = ", PLOT_MODULUS
  print*,"rectangular waveguide"
  print*,"Stimulus = ", FREQUENCY, " Hertz continuous plane wave"
  print*," "
  print*,"Meshing parameters:"
  print '(I4,"x",I3,"x",I3," cells")', nx, ny, nz
  print '(" dx=",F11.8," dy=",F11.8," dz=",F11.8," meters")', dx, dy, dz
  print '(F11.8, " x ", F11.8, " x ", F11.8, " meter^3 simulation region")', &
       & GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda
  print*," "
  print '(" Time simulated will be ",E12.5E2," seconds, ",I5.5," timesteps")', &
       & totalSimulatedTime, MAXIMUM_ITERATION
  print*," "
  print*,"3D output scaling parameters:"
  print*,"Autoscaling every timestep"
  print*," "
  print*," "
  !! print out some info on each timestep: 
  print*,"Following is the iteration number and current "
  print*,"simulated time for each timestep/iteration of "
  print*,"the simulation.  For each timestep that data is "
  print*,"output to file, the maximum and minimum data "
  print*,"values are printed here with the maximum and "
  print*,"minimum scaled values in parentheses. "
  print*," "

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! open and start writing the .viz file
  !!    this file will be handy to feed parameters to viz if desired

  open(unit=33,file="ToyFDTD1f90.viz",iostat=ios)
  !! if the file fails to open, print an error message to 
  !!     standard output:
  if(ios /= 0) then
    print*,"Difficulty opening ToyFDTD1f90.viz"
    STOP
  end if
  write(unit=33, fmt='("#Viz V1.0")')
  write(unit=33, fmt='("time: ",F5.3,ES12.5E2)') currentSimulatedTime,dt
  write(unit=33, fmt='("color: chengGbry.cmap")')
  write(unit=33, fmt='("")')                                     

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! main loop:

  do iteration=0,(MAXIMUM_ITERATION-1),1
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Output section: 

      !! time in simulated seconds that the simulation has progressed:
      currentSimulatedTime = dt*iteration
      !! print to standard output the iteration number 
      !!     and current simulated time:
      write(*,372) iteration, currentSimulatedTime
372   format ("#",i5," ", ES14.5E2, "sec", $)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! 3D data output every PLOT_MODULUS timesteps:
      !!     The first time through the main loop all the data written to 
      !!     file will be zeros.  If anything is nonzero, there's a bug.  :>
 
      if(mod(iteration,PLOT_MODULUS) == 0) then

         !! Construct the filename:
         !! The numerical part is obtained by successively
         !!      dividing the current iteration by ten, and the
         !!      remainders of these divisions are the digits.
         !!      For example, 12345 would be converted as follows:
         !!      12345/10 = 1234+5/10 (5 is the first digit in the number)
         !!      1234/10 = 123 + 4/10 (4 is the next digit)
         !!      ....and so on..this is performed six times to create
         !!      a filename out of any 6-digit number.
         temp = iteration
         filename(1:4) = "f90_"
         filename(10:10) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(9:9) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(8:8) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(7:7) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(6:6) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(5:5) = numbers(mod(temp,10)+1:mod(temp,10)+1)
         temp = int(temp/10)
         filename(11:14) = ".bob"

         !! open a new data file for this iteration:
         open(unit=69,file=filename,access="direct",iostat=ios,recl=1)
         !! if the file fails to open, print an error message to 
	 !!     standard output:
         if(ios /= 0) then
           print*,"Difficulty opening ",filename
           STOP
         end if

         !! find the max and min values to be output this timestep:  
         min = FLT_MAX
         max = -FLT_MAX
         do i=0,(nx),1
            do j=0,(ny),1
               do k=0,(nz-1),1
                  if(ez(i,j,k) < min) then
                     min = ez(i,j,k)
                  end if
                  if(ez(i,j,k) > max) then
                     max = ez(i,j,k)
                  end if
               end do
            end do
         end do

         !!  update the tracking variables for minimum and maximum 
         !!  values for the entire simulation:
         if(min < simulationMin) then
            simulationMin = min
         end if
         if(max > simulationMax) then
            simulationMax = max
         end if
       
         !! set norm to be max or min, whichever is greater in magnitude:
         if(abs(max) > abs(min)) then
            norm = abs(max)
         else
            norm = abs(min)
         end if

         !! if everything is zero, give norm a tiny value 
	 !!     to avoid division by zero:
         if(norm == 0.0d0) then
            norm = DBL_EPSILON
         end if
         scalingValue = 127.0/norm
         !! write to standard output the minimum and maximum values 
         !!     from this iteration and the minimum and maximum values
         !!     that will be written to the bob file this iteration:
         write(*, 373) min, int(127.0d0 + scalingValue * min), max, int(127.0d0 + scalingValue * max)
373      format(" ", f9.5, " (", i3, ") < ez BoB < ", f9.5, " (", i3, ")", $)

         !! scale each ez value in the mesh to the range of integers 
         !!     from zero through 254 and write them to the output 
         !!     file for this iteration:
         record = 1
         do k=0,(nz-1),1
            do j=0,(ny),1
               do i=0,(nx),1
                  write(unit=69,rec=record) char(int(127.0d0+scalingValue*ez(i,j,k)))
                  record = record + 1
               end do
            end do
         end do
       
	 !! close the output file for this iteration:
         close(unit=69) 
	 !! write the dimensions and name of the output file for this 
	 !!     iteration to the viz file:
         write(unit=33,fmt=334) nx+1,ny+1,nz,filename
334      format(I4.4,"x",I4.4,"x",I4.4," ",A14)
	 !! write identification of the corners of the mesh and the max
	 !!     and min values for this iteration to the viz file:
         write(unit=33,fmt=335) dx*nx,dy*ny,dz*nz,min,max
335      format("bbox: 0.0 0.0 0.0 ",F7.4,F7.4,F7.4," ",F8.4," ",F8.4)

      end if !! end bob output section

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the stimulus: a plane wave emanates from the x=0 face:
      !!     The length of the guide lies in the x-direction, the width of the 
      !!     guide lies in the y-direction, and the height of the guide lies
      !!     in the z-direction.  So the guide is sourced by all the ez 
      !!     components on the stimulus face.  

      stimulus = sin(omega*currentSimulatedTime)
      do i=0,0,1
         do j=0,ny,1
            do k=0,(nz-1),1
               ez(i,j,k) = stimulus
            end do
         end do
      end do
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Update the interior of the mesh:
      !!    all vector components except those on the faces of the mesh.  
      !!
      !! Update all the H field vector components within the mesh:
      !!     Since all H vectors are internal, all H values are updated here.
      !!     Note that the normal H vectors on the faces of the mesh are not 
      !!     computed here, and in fact were never allocated -- the normal 
      !!     H components on the faces of the mesh are never used to update 
      !!     any other value, so they are left out of the memory allocation  
      !!     entirely.   
    
      !! Update the hx values:
      do i=0,(nx-2),1
         do j=0,(ny-1),1
            do k=0,(nz-1),1
               hx(i,j,k) = hx(i,j,k) + dtmudz*(ey(i+1,j,k+1) - ey(i+1,j,k)) &
                                   & - dtmudy*(ez(i+1,j+1,k) - ez(i+1,j,k))
            end do
         end do
      end do
    
      !! Update the hy values:
      do i=0,(nx-1),1
         do j=0,(ny-2),1
            do k=0,(nz-1),1
               hy(i,j,k) = hy(i,j,k) + dtmudx*(ez(i+1,j+1,k) - ez(i,j+1,k)) &
                                   & - dtmudz*(ex(i,j+1,k+1) - ex(i,j+1,k))
            end do
         end do
      end do
    
      !! Update the hz values:
      do i=0,(nx-1),1
         do j=0,(ny-1),1
            do k=0,(nz-2),1
               hz(i,j,k) = hz(i,j,k) + dtmudy*(ex(i,j+1,k+1) - ex(i,j,k+1)) &
                                   & - dtmudx*(ey(i+1,j,k+1) - ey(i,j,k+1))
            end do
         end do
      end do
    
      !! Update the E field vector components.  
      !! The values on the faces of the mesh are not updated here; they 
      !!      are handled by the boundary condition computation 
      !!      (and stimulus computation).  
    
      !! Update the ex values:
      do i=0,(nx-1),1
         do j=1,(ny-1),1
            do k=1,(nz-1),1
               ex(i,j,k) = ex(i,j,k) + dtepsdy*(hz(i,j,k-1) - hz(i,j-1,k-1)) &
                                   & - dtepsdz*(hy(i,j-1,k) - hy(i,j-1,k-1))
            end do
         end do
      end do
    
      !! Update the ey values:
      do i=1,(nx-1),1
         do j=0,(ny-1),1
            do k=1,(nz-1),1
               ey(i,j,k) = ey(i,j,k) + dtepsdz*(hx(i-1,j,k) - hx(i-1,j,k-1)) &
                                   & - dtepsdx*(hz(i,j,k-1) - hz(i-1,j,k-1))
            end do
         end do
      end do
    
      !! Update the ez values:
      do i=1,(nx-1),1
         do j=1,(ny-1),1
            do k=0,(nz-1),1
               ez(i,j,k) = ez(i,j,k) + dtepsdx*(hy(i,j-1,k) - hy(i-1,j-1,k)) &
                                   & - dtepsdy*(hx(i-1,j,k) - hx(i-1,j-1,k))
            end do
         end do
      end do
    
      print*," "
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the boundary conditions:
     
      !! OK, so I'm yanking your chain on this one.  The PEC condition is 
      !! enforced by setting the tangential E field components on all the
      !! faces of the mesh to zero every timestep (except the stimulus 
      !! face).  But the lazy/efficient way out is to initialize those 
      !! vectors to zero and never compute them again, which is exactly 
      !! what happens in this code.  

  end do !! end mainloop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Output section:  
  !! The output routine is repeated one last time to write out  
  !! the last data computed.   
  
  !! time in simulated seconds that the simulation has progressed: 
  currentSimulatedTime = dt*iteration                       
  !! print to standard output the iteration number 
  !!     and current simulated time:
  write(*,372) iteration, currentSimulatedTime

  !! 3D data output for the last timestep: 
  !! create the filename for this iteration, 
  !!     which includes the iteration number:

  !! Construct the filename:
  !! The numerical part is obtained by successively
  !!      dividing the current iteration by ten, and the
  !!      remainders of these divisions are the digits.
  !!      For example, 12345 would be converted as follows:
  !!      12345/10 = 1234+5/10 (5 is the first digit in the number)
  !!      1234/10 = 123 + 4/10 (4 is the next digit)
  !!      ....and so on..this is performed six times to create
  !!      a filename out of any 6-digit number.
  temp = iteration
  filename(1:4) = "f90_"
  filename(10:10) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(9:9) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(8:8) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(7:7) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(6:6) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(5:5) = numbers(mod(temp,10)+1:mod(temp,10)+1)
  temp = int(temp/10)
  filename(11:14) = ".bob"
  
  !! open a new data file for this iteration:
  open(unit=69,file=filename,access="direct",iostat=ios,recl=1)
  !! if the file fails to open, print an error message to 
  !!     standard output:
  if(ios /= 0) then
    print*,"Difficulty opening ",filename
    STOP
  end if
  
  !!! find the min and max values for this iteration
  min = FLT_MAX
  max = -FLT_MAX
  do i=0,(nx),1
     do j=0,(ny),1
        do k=0,(nz-1),1
           if(ez(i,j,k) < min) then
              min = ez(i,j,k)
           end if
           if(ez(i,j,k) > max) then
              max = ez(i,j,k)
           end if
        end do
     end do
  end do
  
  
  !!  update the tracking variables for minimum and maximum 
  !!  values for the entire simulation:
  if(min < simulationMin) then
     simulationMin = min
  end if
  if(max > simulationMax) then
     simulationMax = max
  end if
  
  !! set norm to be max or min, whichever is greater in magnitude:
  if(abs(max) > abs(min)) then
     norm = abs(max)
  else
     norm = abs(min)
  end if
  
  !! if everything is zero, give norm a tiny value 
  !!     to avoid division by zero:
  if(norm == 0.0) then
     norm = DBL_EPSILON
  end if
  scalingValue = 127.0d0/norm
  !! write to standard output the minimum and maximum values 
  !!     from this iteration and the minimum and maximum values
  !!     that will be written to the bob file this iteration:
  write(*, 373) min, int(127.0d0 + scalingValue * min), max, int(127.0d0 + scalingValue * max)
  
  !! scale each ez value in the mesh to the range of integers 
  !!     from zero through 254 and write them to the output 
  !!     file for this iteration:
  record = 1
  do k=0,(nz-1),1
     do j=0,(ny),1
        do i=0,(nx),1
           write(unit=69,rec=record) char(int(127.0d0+scalingValue*ez(i,j,k)))
           record = record + 1
        end do
     end do
  end do
  
  !! close the output file for this iteration:
  close(unit=69) 
  !! write the dimensions and name of the output file for this 
  !!     iteration to the viz file:
  write(unit=33,fmt=334) nx+1,ny+1,nz,filename
  !! write identification of the corners of the mesh and the max
  !!     and min values for this iteration to the viz file:
  write(unit=33,fmt=335) dx*nx,dy*ny,dz*nz,min,max

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! close the viz file for this simulation:
  close(unit=33)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! write some progress notes to standard output:

  print*," "
  print*," "
  !! print out how much memory has been allocated: 
  print*,"Dynamically allocated ", allocatedBytes, " bytes"
  !! print out some simulation parameters: 
  print*,"PLOT_MODULUS = ", PLOT_MODULUS
  print*,"rectangular waveguide"
  print*,"Stimulus = ", FREQUENCY, " Hertz continuous plane wave"
  print*," "
  print*,"Meshing parameters:"
  print '(I4,"x",I3,"x",I3," cells")', nx, ny, nz
  print '(" dx=",F11.8," dy=",F11.8," dz=",F11.8," meters")', dx, dy, dz
  print '(F11.8, " x ", F11.8, " x ", F11.8, " meter^3 simulation region")', &
       & GUIDE_WIDTH, GUIDE_HEIGHT, LENGTH_IN_WAVELENGTHS*lambda
  print*," "
  print '(" Time simulated was ",E12.5E2," seconds, ",I5.5," timesteps")', &
       & totalSimulatedTime, MAXIMUM_ITERATION
  print*," "
  print*,"3D output scaling parameters:"
  print*,"Autoscaling every timestep"
  print*," "
  print*," "
  !! print out simulationMin and simulationMax: 
  print'(" Minimum output value was : ", f10.5)', simulationMin
  print'(" Maximum output value was : ", f10.5)', simulationMax
  print*," "
  print*," "
  !! print out a bob command line, including the dimensions
  !!      of the output files:
  print'(" bob -cmap chengGbry.cmap -s ",I5.5,"x",I5.5,"x",I5.5," *.bob")',nx+1,ny+1,nz
  print*," "
  !! print out a viz command line:
  print*,"viz ToyFDTD1f90.viz"
  print*," "
  print*," "
  print*," "
  print*," "

end program ToyFDTD1 !! end main	



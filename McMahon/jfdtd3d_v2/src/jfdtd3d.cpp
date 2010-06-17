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
// NAME: JFDTD3D: A 3D FDTD CODE
//
// VERSION: 1.00
//		Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jfdtd3d.cpp
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// NOTES:	i. last updated on 7/14/08
//		ii. the code (format and actual coding) is very sloppy, but fully functional
//		iii.a. the (time) staggering is as follows:
//			E(t-1/2) --> E(t+1/2)
//			H(t) --> H(t+1)	
//
//			J(t-1/2) = J(t+1/2)
//			rho(t) --> rho(t+1)
//                      phi(t+1/2) = f( rho(t+1), rho(t) )
//			A(t) --> A(t+1)
//	
//		iii.b. the (spatial) staggering is:
//			Ex (i+1/2, j, k)
//			Ey (i, j+1/2, k)
//			Ez (i, j, k+1/2)
//
//			Hx (i, j+1/2, k+1/2)
//			Hy (i+1/2, j, k+1/2)
//			Hz (i+1/2, j+1/2, k)
// 		iv. it is set up such that mpi_master_id == 0 (for loops and communication
//		purposes) 
//		v. everything is set up in meters and seconds, etc except for the reading in
//		of particles information (but this is immediately converted) and calculations
//		of tetrahedra volumes
//		vi. there code is such that there is always periodicity, but PML can
//		be inserted to mimic non-periodicity
//		vii. the program essentially runs as follows:
//			vii.a. the code is copied to every process
//			vii.b. MPI is initialized at the beginning of main():
//				vii.b.1. the number of processes is found, and a processor ID is
//				assigned to each process
//			vii.c. mpi_initialize() is called which finishes setting up MPI:
//				vii.c.1. the MPI cartesian grid is setup, and the grid coordinates
//				of each process is assigned
//			vii.d. get_processor_responsibilities() finds out what exactly each process
//			is responsible for:
//				vii.d.1. find out the range of grid points the process handles
//				vii.d.2. find out if the process is responsible for TFSF or an output
//				plane
//				vii.d.3. find out what processes this process needs to send and receive
//				from during each time step
//			vii.e. perform_fdtd() is called which runs the actual FDTD algorithm:
//				vii.e.1. each process calculates it's responsible grid points and then
//				sends and receives values the processes set up in 
//				get_processor_responsibilities()
//			vii.f. free_memory() and MPI_Finalize() are called which cleans up and
//			finalizes MPI
//		viii. because the dielectric function evaluates to a NEGATIVE imaginary part
//		the Fourier transforms have to be such that they are done with exp(-i*omega*t),
//		and the time component has exp(i*omega*t)
//
//
// TODO:	i. !!! change loop iteration to count down rather than up
//		ii. !!! change loop check from >= to >
//		iii. !!! possibly implement the centering of grid components for output
//		iv. !!! delete arrays (pml psi's, scat arrays, fields and materials,
//		v. !!! NSOM not fully finished for a number of reasons and has not been tested
//		vi. !!! fix things so we don't need to calculate "local" coordinates
//		on each loop iteration
//		vii. !! in general it seems as if points for things lie on the boundary between 
//		two processors then it seems as if weird things can occur
//		viii. !!! each processor does not NEED to handle the entire auxillary grid
//		ix. !!! put in absorption and extinction cross section calcs
//		x. !!! maybe calculate all incident fields to add in TF/SF point up front
//
//		!!!!!!!!!!!!
//		ii. !!! there is a serious error in the PML, because setting izperiodic == 0 creates
// 		erroneous solutions
// 		iii. !!! so (see ii.) for the time being we will assume there is periodicity in all directions,
// 		except when the periodicity flags are set PML is placed along those sides.	
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

using namespace std;	

//************************************************************************
//			INCLUDE FILES
//************************************************************************	

// MPI
#include <mpi.h>

// STANDARD
// utilities
#include <cstdlib>
#include <ctime>
// math
#include <cmath>
#include <complex>
// I/O
#include <iostream>
#include <fstream> 

// CUSTOM
#include "global.h"

#include "jio.h"

// vector structures and operations
#include "jmath_vector.h"
// elemental structures, etc
#include "jmath_mesh.h"
// Drude + 2 Lorentz pole modes
#include "jmaterials.h"
// standard, complex, and physical constants
#include "jconsts.h"

//************************************************************************
//			PROGRAM CONSTANTS
//************************************************************************


//************************************************************************
//			CUSTOM STRUCTURES
//************************************************************************


//************************************************************************
//			SUBROUTINE DECLARATIONS
//************************************************************************

// INITIALIZATION
void read_parameters();	
void mpi_initialize(int argc, char** argv);
void simulation_initialize();

// STRUCTURE SET-UP, ETC
void init_pml();
void structure_initialize();
extern void get_structure(double xpos, double ypos, double zpos, int& mat_num, double& mat_eps, double& mat_mu);

// FDTD
void perform_fdtd();
void field_e_aux_update();
void field_h_aux_update();
void field_e_update();
void field_h_update();

// OUTPUT
void output(int ioutput_type, int iplane, int iunit);
void ft_fields();

// POST-PROCESSING
// reflection/transmission spectra
void ft_fields_spect();
void calc_transm();
// cross sections
void ft_fields_scat();
void calc_cs();
// NSOM
void ft_field_nsom();
void calc_transm_nsom();
void calc_transm_nf();

// CLEANUP, ETC
void free_memory();


//************************************************************************
//			GLOBAL VARIABLES
//************************************************************************

// *********************************
/*
// POINTS
int grid_nx, grid_ny, grid_nz;
int my_grid_nx, my_grid_ny, my_grid_nz, mpi_my_grid_nx1, mpi_my_grid_nx2, mpi_my_grid_ny1, mpi_my_grid_ny2, mpi_my_grid_nz1, mpi_my_grid_nz2;

// SIZES & SPACINGS
double grid_xsize, grid_ysize, grid_zsize, grid_dx, grid_dy, grid_dz;

// MATERIALS
int ***grid_material_ex, ***grid_material_ey, ***grid_material_ez;
double ***permittivity_ex, ***permittivity_ey, ***permittivity_ez, ***permeability_hx, ***permeability_hy, ***permeability_hz;
*/
// *********************************



//=========================================================
// MPI
//=========================================================

// NUMBER OF PROCESSES & IDs
int mpi_nprocessors, mpi_master_id, mpi_my_id;

// CARTESIAN GRID
// communicator
MPI_Comm mpi_grid_comm;
// sizes, periodicity flags, and IDs
int mpi_nxprocs, mpi_nyprocs, mpi_nzprocs, mpi_grid_dim_sizes[3], mpi_grid_iperiodicity[3], mpi_my_grid_coords[3];
// send/recieve info
int mpi_efield_isend, mpi_efield_irecv, mpi_efield_jsend, mpi_efield_jrecv, mpi_efield_ksend, mpi_efield_krecv, mpi_hfield_isend, mpi_hfield_irecv, mpi_hfield_jsend, mpi_hfield_jrecv, mpi_hfield_ksend, mpi_hfield_krecv;

// COMMUNICATION METHOD
int isendtype;


//=========================================================
// SYSTEM
//=========================================================

//---------------------------------------------------------
// TIMINGS
//---------------------------------------------------------
time_t systime_total, systime_start, systime_end;


//=========================================================
// COMPUTATIONAL DOMAIN & RELATED  
//=========================================================

//---------------------------------------------------------
// PERIODICITY
//---------------------------------------------------------
int ixperiodic, iyperiodic, izperiodic;

//---------------------------------------------------------
// COMPUTATIONAL GRID
//---------------------------------------------------------

// MPI
//int my_grid_nx, my_grid_ny, my_grid_nz, mpi_my_grid_nx1, mpi_my_grid_nx2, mpi_my_grid_ny1, mpi_my_grid_ny2, mpi_my_grid_nz1, mpi_my_grid_nz2;

//---------------------------------------------------------
// TIME & STEPS
//---------------------------------------------------------
double time_total, time_dt, time_courant_factor, time_current;
int nsteps, istep;

//---------------------------------------------------------
// SOURCE
//---------------------------------------------------------

// TF/SF
double tfsf_zpos;
int tfsf_depth;
// MPI 
int mpi_itfsf_k1, mpi_tfsf_k1;

// SOURCE
double source_intensity, source_wavelength, source_wavenumber, source_gauss_center, source_gauss_width;
// polarization
// !fix
int src_ieypol, src_iexpol, src_i45pol, src_icircpol;

//---------------------------------------------------------
// FIELDS & MATERIALS
//---------------------------------------------------------

// FIELDS
// i. the nm1 terms are for metal updates 
double ***field_ex, ***field_ey, ***field_ez, ***field_ex_nm1, ***field_ey_nm1, ***field_ez_nm1, ***field_hx, ***field_hy, ***field_hz;
// MPI send and receive planes
int ksurface_nvals, jsurface_nvals, isurface_nvals;
double *ksurface_x_send, *ksurface_y_send, *jsurface_x_send, *jsurface_z_send, *isurface_y_send, *isurface_z_send, *isurface_y_recv, *isurface_z_recv, *jsurface_x_recv, *jsurface_z_recv, *ksurface_x_recv, *ksurface_y_recv;

//---------------------------------------------------------
// CPML
//---------------------------------------------------------

int cpml_layers, cpml_m, cpml_ma;
double cpml_epsr, cpml_mur, cpml_kappamax, cpml_sigmamax_coeff, cpml_alphamax;
double ***psi_eyx, ***psi_ezx, ***psi_hyx, ***psi_hzx, ***psi_exy, ***psi_ezy, ***psi_hxy, ***psi_hzy, ***psi_exz, ***psi_eyz, ***psi_hxz, ***psi_hyz, *kedx, *khdx, *kedy, *khdy, *kedz, *khdz, *be_x, *ce_x, *bh_x, *ch_x, *be_y, *ce_y, *bh_y, *ch_y, *be_z, *ce_z, *bh_z, *ch_z;

//---------------------------------------------------------
// AUXILLARY GRID
//---------------------------------------------------------

// SIZE & IMPORTANT POINTS
int aux_npts, aux_tfsf_k1, aux_spect_transm_kpt, *ft_aux_coord;

// COEFFICIENTS & FIELDS
double aux_epsr, aux_mur, e_aux_coeff, h_aux_coeff, *ey_aux, *ex_aux, *hx_aux, *hy_aux;

//---------------------------------------------------------------
// OUTPUT					
//---------------------------------------------------------------

int noutput, output_nplanes, *output_format, *output_field, *output_plane, *output_plane_local_coord;
double *output_plane_midpos;

// FOURIER TRANSFORMED-PLANES
int nft_planes, *output_to_ft, *ft_to_output;
double *ft_wavelengths, *ft_wavenumbers;
complex<double> ***field_ft_ex, ***field_ft_ey, ***field_ft_ez, ***field_ft_hx, ***field_ft_hy, ***field_ft_hz, *aux_ftex, *aux_ftey, *aux_ftez, *aux_fthx, *aux_fthy, *aux_fthz; 
// MPI
int mpi_ift, *mpi_ioutput_plane;

// SIMFILE
int isimfile;
ofstream simfile;


//=========================================================
// SPECTRA & POST-PROCESSING CALCULATIONS
//=========================================================

// FOURIER TRANSNFORM "ROTATORS"
complex<double> **dft_multiplier_t, **dft_multiplier_ht;

//---------------------------------------------------------
// TRANSMISSION/REFLECTION SPECTRA
//---------------------------------------------------------
int itransm, itransm_nf;

int spect_npts, spect_transm_kpt, spect_refl_kpt;
double spect_transm_pos, spect_refl_pos, nf_aperture_xcenter, nf_aperture_ycenter, nf_aperture_size;
double spect_minwave, spect_maxwave, *spect_wavenumbers;
complex<double> ***transm_ftex, ***transm_ftey, ***transm_fthx, ***transm_fthy, ***refl_ftex, ***refl_ftey, ***refl_fthx, ***refl_fthy, *aux_ftex_transm, *aux_ftey_transm, *aux_ftez_transm, *aux_fthx_transm, *aux_fthy_transm, *aux_fthz_transm, *aux_ftex_refl, *aux_ftey_refl, *aux_ftez_refl, *aux_fthx_refl, *aux_fthy_refl, *aux_fthz_refl;
// MPI
int mpi_ft_spect_transm, mpi_ft_spect_refl;

//---------------------------------------------------------
// CROSS SECTIONS
//---------------------------------------------------------
// !q: are the scat_i1ft, ... MPI ???
int iscat_calc;
double scat_xcenter, scat_ycenter, scat_zcenter, scat_radius;
// FLAGS & RANGES
int scat_i1ft, scat_i2ft, scat_j1ft, scat_j2ft, scat_k1ft, scat_k2ft;
int scat_i1pt, scat_i2pt, scat_j1pt, scat_j2pt, scat_k1pt, scat_k2pt;
int scat_i1y1, scat_i1y2, scat_i1z1, scat_i1z2, scat_i2y1, scat_i2y2, scat_i2z1, scat_i2z2, scat_j1x1, scat_j1x2, scat_j1z1, scat_j1z2, scat_j2x1, scat_j2x2, scat_j2z1, scat_j2z2, scat_k1y1, scat_k1y2, scat_k1x1, scat_k1x2, scat_k2y1, scat_k2y2, scat_k2x1, scat_k2x2;
complex<double> ***scat_i1_ftez, ***scat_i1_ftey, ***scat_i1_fthz, ***scat_i1_fthy, ***scat_i2_ftez, ***scat_i2_ftey, ***scat_i2_fthz, ***scat_i2_fthy, ***scat_j1_ftex, ***scat_j1_ftez, ***scat_j1_fthx, ***scat_j1_fthz, ***scat_j2_ftex, ***scat_j2_ftez, ***scat_j2_fthx, ***scat_j2_fthz, ***scat_k1_ftex, ***scat_k1_ftey, ***scat_k1_fthx, ***scat_k1_fthy, ***scat_k2_ftex, ***scat_k2_ftey, ***scat_k2_fthx, ***scat_k2_fthy; 
complex<double> **scat_aux_ex, **scat_aux_ey, **scat_aux_hx, **scat_aux_hy;

//---------------------------------------------------------
// NSOM PROBE
//---------------------------------------------------------
// !q: are the nsom_ift, ... MPI ???
// !!! experimental

// MPI
//int nsom_ift;

// TRANSMISSION CHARACTERISTICS
//int transm_nsom_npts, transm_nsom_kpt, *transm_nsom_ptoi, *transm_nsom_ptoj, transm_nsom_ptok;
//double transm_nsom_zpos, nsom_wavelength, nsom_wavenumber;
complex<double> *transm_nsom_ex, *transm_nsom_ey, *transm_nsom_hx, *transm_nsom_hy; 
complex<double> aux_ftex_nsom, aux_ftey_nsom, aux_fthx_nsom, aux_fthy_nsom;


//=========================================================
// DRUDE + 2 LORENTZ POLE MATERIALS
//=========================================================

// CURRENTS & POLARIZATION TERMS
double *metal_c1, *metal_c2, *metal_c3, *drude_kp, *drude_betap, **lorentz_alphap, **lorentz_zetap, **lorentz_gammap, ***drude_current_x, ***drude_current_y, ***drude_current_z, ****lorentz_currentp_x, ****lorentz_currentp_y, ****lorentz_currentp_z,****lorentz_currentp_x_nm1, ****lorentz_currentp_y_nm1, ****lorentz_currentp_z_nm1;


//************************************************************************
//				SUBROUTINES
//************************************************************************


//========================================================================
//========================================================================
//
// 	NAME:	int main(int argc, char** argv)
//	DESC:	Entryway to program.
//
//	NOTES:
//
//========================================================================
//========================================================================
int main(int argc, char** argv)
{

  // INITIALIZE MPI
  // i. some initial MPI setup is done here first thing, since some odd behavior
  // was observed in JFEM3D when read_parameters() was before MPI_Init
  MPI_Init(&argc, &argv);

  //=========================================================
  // INITIALIZATION
  //=========================================================

  // RIGHT OFF THE BAT GET THE SYSTEM TIME
  // i. this placement will not catch the time taken to initialize MPI
  systime_start = time(NULL);


  // FIRST READ IN PARAMETERS BECAUSE WE NEED SOME OF IT FOR SETTING UP MPI
  read_parameters();

  // INITIALIZE AND SET UP MPI
  mpi_initialize(argc, argv);

  // PERFORM SIMULATION SET-UP
  simulation_initialize();

  // SET UP PML
  init_pml();

  // GET STRUCTURAL INFO
  structure_initialize();


  //=========================================================
  // FDTD
  //=========================================================
  perform_fdtd();


  //=========================================================
  // POST-PROCESSING
  //=========================================================

  cout << "here 1" << endl;

  // CALCULATE THE TRANSMISSION
  if(itransm == 1)
  {
    calc_transm();
  }

  cout << "here 2" << endl;

  // CALCULATE THE NEAR-FIELD TRANSMISSION
  if(itransm_nf == 1)
  {
    calc_transm_nf();
  }

  cout << "here 3" << endl;

  // CALCULATE CROSS SECTIONS
  if(iscat_calc == 1)
  {
    calc_cs();
  }

  cout << "here 4" << endl;

  // CALCULATE NSOM POWER FLOW
  if(insom == 1)
  {
    calc_transm_nsom();
  }


  cout << "here 5" << endl;

  // GET THE SYSTEM TIME
  // i. this neglects the final output and final cleanup times
  systime_total = time(NULL) - systime_start;


  // OUTPUT FINAL DATA FOR SIMULATION INFO
  if(mpi_my_id == mpi_master_id)
  {
    output(3, 0, 0);
  }

  cout << "here 6" << endl;

  //=========================================================
  // FINALIZE SIMULATION
  //=========================================================

  // FREE DYNAMIC MEMORY
  //free_memory();

  // SHUT DOWN MPI
  // i. the last thing for MPI simulatons is to shut it down
  MPI_Finalize();


  return 0;

}



//========================================================================
//========================================================================
//
//	NAME:	void mpi_initialize(int argc, char** argv)
//
//	DESC:	Initialized MPI settings.
//
//	NOTES: 	i. 
//
//
//========================================================================
//========================================================================
void mpi_initialize(int argc, char** argv)
{

  // local indices
  int i, ii;

  //---------------------------------------------------------
  // SET-UP BASIC MPI
  //---------------------------------------------------------

  // GET THE NUMBER OF PROCESSES
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocessors);

  // GET THE ID OF THE PROCESSOR RUNNING THIS PROGRAM
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_id);


  // CHECK THAT THE DESIRED NUMBER OF PROCESSORS HAS BEEN ALLOCATED
  // i. in reality we do not HAVE to exit here, but setting up the responisibility
  // grid later may be difficult if we have done something wrong
  if( (mpi_nxprocs*mpi_nyprocs*mpi_nzprocs) != mpi_nprocessors)
  {
    cout << endl << " !! NUMBER OF REQUESTED PROCESSORS IS NOT EQUAL " << endl;
    cout << " THE NUMBER OF PROCESSORS AVAILABLE !!" << endl;

    cout << "requested: " << (mpi_nxprocs*mpi_nyprocs*mpi_nzprocs) << endl;
    cout << "available: " << mpi_nprocessors << endl;

    exit(1);
  }


  // SPECIFY MASTER ID
  mpi_master_id = 0;


  //---------------------------------------------------------
  // SET-UP MPI CARTESIAN TOPOLOGY
  //---------------------------------------------------------

  int mpi_my_grid_id;

  // SET NUMBER OF DIMENSIONS
  int mpi_grid_ndims = 3;

  // SPECIFY REORDORING OF RANKS IS ALLOWED
  int mpi_grid_ireorder = 1;


  // ASSIGN THE NUMBER OF GRID BLOCKS IN EACH DIRECTION
  // i. we only do this step because an array is needed for MPI_Cart_create
  mpi_grid_dim_sizes[0] = mpi_nxprocs;
  mpi_grid_dim_sizes[1] = mpi_nyprocs;
  mpi_grid_dim_sizes[2] = mpi_nzprocs;

  // SET WRAP AROUND IN ALL DIMENSIONS
  mpi_grid_iperiodicity[0] = mpi_grid_iperiodicity[1] = mpi_grid_iperiodicity[2] = 1;

  // NOW CREATE THE CARTESIAN TOPOLOGY
  MPI_Cart_create(MPI_COMM_WORLD, mpi_grid_ndims, mpi_grid_dim_sizes, mpi_grid_iperiodicity, mpi_grid_ireorder, &mpi_grid_comm);

  // GET THE PROCESS(OR) ID WITHIN OUR NEW COMMUNICATOR
  // i.) see the note with the last Comm_rank call
  MPI_Comm_rank(mpi_grid_comm, &mpi_my_grid_id);

  // NOW GET THE CARTESIAN TOPOLOGY COORDINATES FOR THIS PROCESS(OR)
  // i.) note the topology goes from 0, 1, ..., mpi_grid_ndims-1
  MPI_Cart_coords(mpi_grid_comm, mpi_my_grid_id, mpi_grid_ndims, mpi_my_grid_coords);


  //---------------------------------------------------------
  // DETERMINE COMMUNICATION INFO
  //---------------------------------------------------------

  int mpi_coord_send[3], mpi_coord_recv[3], mpi_proc_send[3], mpi_proc_recv[3];

  // E-FIELD
  // i. for the e-field we send to "below" and receive from "above"

  for(i = 0; i <= (3-1); ++i)
  {
    for(ii = 0; ii <= (3-1); ++ii)
    {
      mpi_coord_send[ii] = mpi_my_grid_coords[ii];
      mpi_coord_recv[ii] = mpi_my_grid_coords[ii];
    }

    // SEND BELOW
    // make sure we wrap around ...
    if(mpi_coord_send[i] == 0)
    {
      mpi_coord_send[i] = mpi_grid_dim_sizes[i] - 1;
    }
    // ... otherwise send to the proc right below us
    else
    {
      mpi_coord_send[i] -= 1;
    }
  
    // RECEIVE FROM ABOVE
    // make sure we wrap around ...
    if(mpi_coord_recv[i] == (mpi_grid_dim_sizes[i] - 1))
    {
      mpi_coord_recv[i] = 0;
    }
    // ... otherwise receieve from the proc right above us
    else
    {
      mpi_coord_recv[i] += 1;
    }

    // NOW GET THE RANKS FROM THE COORDINATE INFORMATION
    MPI_Cart_rank(mpi_grid_comm, mpi_coord_send, &mpi_proc_send[i]);
    MPI_Cart_rank(mpi_grid_comm, mpi_coord_recv, &mpi_proc_recv[i]);
  } // ++i

  // ASSIGN
  mpi_efield_isend = mpi_proc_send[0];
  mpi_efield_irecv = mpi_proc_recv[0];

  mpi_efield_jsend = mpi_proc_send[1];
  mpi_efield_jrecv = mpi_proc_recv[1];

  mpi_efield_ksend = mpi_proc_send[2];
  mpi_efield_krecv = mpi_proc_recv[2];


  // H-FIELD
  // i. for the h-field we send "above" and receive from "below"

  for(i = 0; i <= (3-1); ++i)
  {
    for(ii = 0; ii <= (3-1); ++ii)
    {
      mpi_coord_send[ii] = mpi_my_grid_coords[ii];
      mpi_coord_recv[ii] = mpi_my_grid_coords[ii];
    }

    // SEND ABOVE
    // make sure we wrap around ...
    if(mpi_coord_send[i] == (mpi_grid_dim_sizes[i]-1))
    {
      mpi_coord_send[i] = 0;
    }
    // ... otherwise send to the proc right above us
    else
    {
      mpi_coord_send[i] += 1;
    }
  
    // RECEIVE FROM BELOW
    // make sure we wrap around ... 
    if(mpi_coord_recv[i] == 0)
    {
      mpi_coord_recv[i] = mpi_grid_dim_sizes[i] - 1;
    }
    // ... otherwise receieve from the proc right below us
    else
    {
      mpi_coord_recv[i] -= 1;
    }

    // NOW GET THE RANKS FROM THE COORDINATE INFORMATION
    MPI_Cart_rank(mpi_grid_comm, mpi_coord_send, &mpi_proc_send[i]);
    MPI_Cart_rank(mpi_grid_comm, mpi_coord_recv, &mpi_proc_recv[i]);
  } // i++

  // ASSIGN
  mpi_hfield_isend = mpi_proc_send[0];
  mpi_hfield_irecv = mpi_proc_recv[0];

  mpi_hfield_jsend = mpi_proc_send[1];
  mpi_hfield_jrecv = mpi_proc_recv[1];

  mpi_hfield_ksend = mpi_proc_send[2];
  mpi_hfield_krecv = mpi_proc_recv[2];


  //=========================================================
  // SIMULATION FILE
  //=========================================================
  // i. for debugging purposes it would be ideal to create the simulation file 
  // at the very beginning, but since MPI information is needed for this we
  // will do this here
  // !fix

  if(isimfile==1)
  {
    // SET-UP FILENAME
    char simfilename_1[] = "./output/sim_file.";
    char simfilename_2[5];
    char simfilename[15];

    sprintf(simfilename_2,"%i",mpi_my_id); 
    strcpy(simfilename, simfilename_1);
    strcat(simfilename, simfilename_2);

    // OPEN FILE
    simfile.open(simfilename);

    // WRITE INITIAL INFORMATION
    simfile << "In COMM_WORLD I am processor numbered: " << mpi_my_id << endl;
    simfile << "In COMM_GRID I am numbered: " << endl;
    simfile << " -mpi_my_grid_coords[0] (x): " << mpi_my_grid_coords[0] << endl;
    simfile << " -mpi_my_grid_coords[1] (y): " << mpi_my_grid_coords[1] << endl;
    simfile << " -mpi_my_grid_coords[2] (z): " << mpi_my_grid_coords[2] << endl;
    simfile << "***********************************" << endl << endl;
  }


  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void simulation_initialize()
//	DESC:	Initializes simulation settings and allocate memory.
//
//	NOTES: 	i. this subroutine is distinctly different than the 2D code
//
//
//========================================================================
//========================================================================
void simulation_initialize()
{

  // local indices
  int i, j, k, l, m, n, p;


  //=========================================================
  // SET-UP GRID
  //=========================================================

  // GET NUMBER OF GRID POINTS
  grid_nx = static_cast<int>(grid_xsize/grid_dx);
  grid_ny = static_cast<int>(grid_ysize/grid_dy);
  grid_nz = static_cast<int>(grid_zsize/grid_dz);

  // ADJUST GRID SIZES 
  grid_xsize = (static_cast<double>(grid_nx))*grid_dx;
  grid_ysize = (static_cast<double>(grid_ny))*grid_dy;
  grid_zsize = (static_cast<double>(grid_nz))*grid_dz;

  //---------------------------------------------------------
  // ADJUST FOR NSOM
  //---------------------------------------------------------
  // i. this is here so we can resize the grid so to truncate the NSOM probe 
  // (if there is one) with PML
  // !fix
  if(insom == 1)
  {
    // RESIZE THE GRID
    grid_zsize = nsom_zcenter + nsom_tip_height + nsom_cylinder_height;
    // i. we use -1 to make sure to truncate the nsom probe with PML
    grid_nz = static_cast<int>(grid_zsize/grid_dz) - 1;
 
    // ALLOCATE MEMORY & INITIALIZE
    // i. !!! this allocation of memory is WAY more than enough because not every
    // process needs this, but will have it, and the size of arrays is for a square
    // with sides as long as the diameter of.  In addition it would make more sense to
    // ? put this somewhere else
    int max_nsom_pts = (static_cast<int>((nsom_cylinder_width+nsom_cladding_width)/grid_dy))*(static_cast<int>((nsom_cylinder_width+nsom_cladding_width)/grid_dx));
    transm_nsom_ptoi = new int [max_nsom_pts+1];
    transm_nsom_ptoj = new int [max_nsom_pts+1];
    transm_nsom_ex = new complex<double> [max_nsom_pts+1];
    transm_nsom_ey = new complex<double> [max_nsom_pts+1];
    transm_nsom_hx = new complex<double> [max_nsom_pts+1];
    transm_nsom_hy = new complex<double> [max_nsom_pts+1];

    for(i = 0; i <= max_nsom_pts; ++i)
    {
      transm_nsom_ex[i] = COMPLEXZERO;
      transm_nsom_ey[i] = COMPLEXZERO;
      transm_nsom_hx[i] = COMPLEXZERO;
      transm_nsom_hy[i] = COMPLEXZERO;
    }

    aux_ftex_nsom = COMPLEXZERO;
    aux_ftey_nsom = COMPLEXZERO;
    aux_fthx_nsom = COMPLEXZERO;
    aux_fthy_nsom = COMPLEXZERO;


    // SET UP THE NSOM WAVENUMBER 
    nsom_wavenumber = 2.0*PI*CSPEED/nsom_wavelength;


    // GET THE POSITION TO CALCULATE THE FIELD ACROSS THE PROBE

    // -- right at the tip --
    // we add 1 to make sure that we end up INSIDE the tip
    //transm_nsom_kpt = static_cast<int>(nsom_zcenter/grid_dz) + 1;

    // -- right above the opening in the cylinder --
    // I set it to be only 1/10 into the cylinder and then add 1
    transm_nsom_zpos = nsom_zcenter + nsom_tip_height + nsom_cylinder_height/10.0;
    transm_nsom_kpt = static_cast<int>(transm_nsom_zpos/grid_dz) + 1;

    // SET NOW TO NOT DO ANY FT ON THE NSOM
    // i. this will be changed when we initialize the structure
    nsom_ift = 0;
    transm_nsom_npts = 0;

  } // end if insom==1


  //---------------------------------------------------------
  // GET PROCESSOR GRID 
  //---------------------------------------------------------

  // GET THE (ROUNDED) NUMBER OF GRID POINTS TO HANDLE
  int dnx, dny, dnz;
  dnx = grid_nx/mpi_grid_dim_sizes[0];
  dny = grid_ny/mpi_grid_dim_sizes[1];
  dnz = grid_nz/mpi_grid_dim_sizes[2];

  // GET THE REMAINED FROM THE ROUNDED GRID BLOCKS
  int dnx_remainder, dny_remainder, dnz_remainder;
  dnx_remainder = (grid_nx % mpi_grid_dim_sizes[0]);
  dny_remainder = (grid_ny % mpi_grid_dim_sizes[1]);
  dnz_remainder = (grid_nz % mpi_grid_dim_sizes[2]);

  // ASSIGN THE BLOCKS TO THIS PROCESS(OR) 

  //--------------------------- x ---------------------------
  if(dnx_remainder!=0)
  {
    if(mpi_my_grid_coords[0] < dnx_remainder)
    {
      mpi_my_grid_nx1 = mpi_my_grid_coords[0]*(dnx+1) + 1;
      mpi_my_grid_nx2 = (mpi_my_grid_coords[0]+1)*(dnx+1);
    }
    else
    {
      mpi_my_grid_nx1 = mpi_my_grid_coords[0]*dnx + 1 + dnx_remainder;
      mpi_my_grid_nx2 = (mpi_my_grid_coords[0]+1)*dnx + dnx_remainder;
    }
  }
  else
  {
    mpi_my_grid_nx1 = mpi_my_grid_coords[0]*dnx + 1;
    mpi_my_grid_nx2 = (mpi_my_grid_coords[0]+1)*dnx;
  }

  //--------------------------- y ---------------------------
  if(dny_remainder!=0)
  {
    if(mpi_my_grid_coords[1] < dny_remainder)
    {
      mpi_my_grid_ny1 = mpi_my_grid_coords[1]*(dny+1) + 1;
      mpi_my_grid_ny2 = (mpi_my_grid_coords[1]+1)*(dny+1);
    }
    else
    {
      mpi_my_grid_ny1 = mpi_my_grid_coords[1]*dny + 1 + dny_remainder;
      mpi_my_grid_ny2 = (mpi_my_grid_coords[1]+1)*dny + dny_remainder;
    }
  }
  else
  {
    mpi_my_grid_ny1 = mpi_my_grid_coords[1]*dny + 1;
    mpi_my_grid_ny2 = (mpi_my_grid_coords[1]+1)*dny;
  }

  //--------------------------- z ---------------------------
  if(dnz_remainder!=0)
  {
    if(mpi_my_grid_coords[2] < dnz_remainder)
    {
      mpi_my_grid_nz1 = mpi_my_grid_coords[2]*(dnz+1) + 1;
      mpi_my_grid_nz2 = (mpi_my_grid_coords[2]+1)*(dnz+1);
    }
    else
    {
      mpi_my_grid_nz1 = mpi_my_grid_coords[2]*dnz + 1 + dnz_remainder;
      mpi_my_grid_nz2 = (mpi_my_grid_coords[2]+1)*dnz + dnz_remainder;
    }
  }
  else
  {
    mpi_my_grid_nz1 = mpi_my_grid_coords[2]*dnz + 1;
    mpi_my_grid_nz2 = (mpi_my_grid_coords[2]+1)*dnz;
  }

  // MAKE SURE THAT THE LOWEST PROCESSORS HANDLE POINT 0 AND TOP PROCESSOR HANDLES grid_npts+1
  // i. only if there is no periodicity
  // !!! this should probably be done

  // NOW GET TOTAL AMOUNT OF GRID POINTS
  my_grid_nx = mpi_my_grid_nx2 - mpi_my_grid_nx1 + 1;
  my_grid_ny = mpi_my_grid_ny2 - mpi_my_grid_ny1 + 1;
  my_grid_nz = mpi_my_grid_nz2 - mpi_my_grid_nz1 + 1;


  //---------------------------------------------------------
  // ALLOCATE SURFACES TO SEND/RECEIVE TO/FROM OTHER PROCESSORS
  //---------------------------------------------------------

  // SIZES
  ksurface_nvals = (mpi_my_grid_nx2 - mpi_my_grid_nx1 + 1)*(mpi_my_grid_ny2 - mpi_my_grid_ny1 + 1);
  jsurface_nvals = (mpi_my_grid_nx2 - mpi_my_grid_nx1 + 1)*(mpi_my_grid_nz2 - mpi_my_grid_nz1 + 1);
  isurface_nvals = (mpi_my_grid_ny2 - mpi_my_grid_ny1 + 1)*(mpi_my_grid_nz2 - mpi_my_grid_nz1 + 1);

  // ALLOCATE
  ksurface_x_send = new double [ksurface_nvals];
  ksurface_y_send = new double [ksurface_nvals];
  ksurface_x_recv = new double [ksurface_nvals];
  ksurface_y_recv = new double [ksurface_nvals];

  jsurface_x_send = new double [jsurface_nvals];
  jsurface_z_send = new double [jsurface_nvals];
  jsurface_x_recv = new double [jsurface_nvals];
  jsurface_z_recv = new double [jsurface_nvals];

  isurface_y_send = new double [isurface_nvals];
  isurface_z_send = new double [isurface_nvals];
  isurface_y_recv = new double [isurface_nvals];
  isurface_z_recv = new double [isurface_nvals];


  //=========================================================
  // SET-UP TIME AND STEPS
  //=========================================================

  time_current = 0.0;
 
  // GET TIME STEP
  // i. this is a user-specified fraction of the Courant stability limit
  time_dt = time_courant_factor/(CSPEED*sqrt(1.0/(grid_dx*grid_dx) + 1.0/(grid_dy*grid_dy) + 1.0/(grid_dz*grid_dz)));

  // GET NUMBER OF STEPS
  // i. add 1 to make sure the simulation is at least as long as the specified time
  nsteps = static_cast<int>(time_total/time_dt) + 1;
  
  // ADJUST TOTAL TIME
  time_total = (static_cast<double>(nsteps))*time_dt;


  //=========================================================
  // SOURCE SET-UP
  //=========================================================

  // SET WAVENUMBER	
  source_wavenumber = 2.0*PI*CSPEED/source_wavelength;


  //=========================================================
  // GET TF/SF INFO
  //=========================================================
  // i. because of normal incidence we only worry about the lower TF/SF plane

  // GET TF/SF POINT
  mpi_tfsf_k1 = static_cast<int>(tfsf_zpos/grid_dz) + 1;

  // CHECK RESPONSIBILITY
  mpi_itfsf_k1 = 0;
  if((mpi_tfsf_k1 >= mpi_my_grid_nz1) && (mpi_tfsf_k1 <= mpi_my_grid_nz2))
  {
    mpi_itfsf_k1 = 1;
  }


  //=========================================================
  // MATERIALS & FIELDS
  //=========================================================
  // i. the values come from the loops that look like [ for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk) 
  // k = kk - mpi_my_grid_nz1 + 1 ], but we will need 1 additional value over for received values from other
  // processors
  // ii. xxx??? I think the size of the nm1 arrays, permittivities, permeabilities, and lorentz arrays can be reduced by 1

  //---------------------------------------------------------
  // ALOCATE MEMORY & INITIALIZE
  //---------------------------------------------------------
  field_ex=new double**[my_grid_nx+1+1];
  field_ey=new double**[my_grid_nx+1+1];
  field_ez=new double**[my_grid_nx+1+1];
  field_ex_nm1=new double**[my_grid_nx+1+1];
  field_ey_nm1=new double**[my_grid_nx+1+1];
  field_ez_nm1=new double**[my_grid_nx+1+1];
  field_hx=new double**[my_grid_nx+1+1];
  field_hy=new double**[my_grid_nx+1+1];
  field_hz=new double**[my_grid_nx+1+1];

  grid_material_ex=new int**[my_grid_nx+1+1];
  grid_material_ey=new int**[my_grid_nx+1+1];
  grid_material_ez=new int**[my_grid_nx+1+1];

  permittivity_ex=new double**[my_grid_nx+1+1];
  permittivity_ey=new double**[my_grid_nx+1+1];
  permittivity_ez=new double**[my_grid_nx+1+1];
  permeability_hx=new double**[my_grid_nx+1+1];
  permeability_hy=new double**[my_grid_nx+1+1];
  permeability_hz=new double**[my_grid_nx+1+1];

  for(i=0;i<=my_grid_nx+1;++i)
  {
    field_ex[i]=new double*[my_grid_ny+1+1];
    field_ey[i]=new double*[my_grid_ny+1+1];
    field_ez[i]=new double*[my_grid_ny+1+1];
    field_ex_nm1[i]=new double*[my_grid_ny+1+1];
    field_ey_nm1[i]=new double*[my_grid_ny+1+1];
    field_ez_nm1[i]=new double*[my_grid_ny+1+1];
    field_hx[i]=new double*[my_grid_ny+1+1];
    field_hy[i]=new double*[my_grid_ny+1+1];
    field_hz[i]=new double*[my_grid_ny+1+1];

    grid_material_ex[i]=new int*[my_grid_ny+1+1];
    grid_material_ey[i]=new int*[my_grid_ny+1+1];
    grid_material_ez[i]=new int*[my_grid_ny+1+1];

    permittivity_ex[i]=new double*[my_grid_ny+1+1];
    permittivity_ey[i]=new double*[my_grid_ny+1+1];
    permittivity_ez[i]=new double*[my_grid_ny+1+1];
    permeability_hx[i]=new double*[my_grid_ny+1+1];
    permeability_hy[i]=new double*[my_grid_ny+1+1];
    permeability_hz[i]=new double*[my_grid_ny+1+1];

    for(j=0;j<=my_grid_ny+1;++j)
    {
      field_ex[i][j]=new double [my_grid_nz+1+1];
      field_ey[i][j]=new double [my_grid_nz+1+1];
      field_ez[i][j]=new double [my_grid_nz+1+1];
      field_ex_nm1[i][j]=new double [my_grid_nz+1+1];
      field_ey_nm1[i][j]=new double [my_grid_nz+1+1];
      field_ez_nm1[i][j]=new double [my_grid_nz+1+1];
      field_hx[i][j]=new double [my_grid_nz+1+1];
      field_hy[i][j]=new double [my_grid_nz+1+1];
      field_hz[i][j]=new double [my_grid_nz+1+1];

      grid_material_ex[i][j]=new int [my_grid_nz+1+1];
      grid_material_ey[i][j]=new int [my_grid_nz+1+1];
      grid_material_ez[i][j]=new int [my_grid_nz+1+1];

      permittivity_ex[i][j]=new double [my_grid_nz+1+1];
      permittivity_ey[i][j]=new double [my_grid_nz+1+1];
      permittivity_ez[i][j]=new double [my_grid_nz+1+1];
      permeability_hx[i][j]=new double [my_grid_nz+1+1];
      permeability_hy[i][j]=new double [my_grid_nz+1+1];
      permeability_hz[i][j]=new double [my_grid_nz+1+1];

      // zero memory
      for(k=0;k<=my_grid_nz+1;++k)
      {
        field_ex[i][j][k]=0.0;
        field_ey[i][j][k]=0.0;
        field_ez[i][j][k]=0.0;
        field_ex_nm1[i][j][k]=0.0;
        field_ey_nm1[i][j][k]=0.0;
        field_ez_nm1[i][j][k]=0.0;
        field_hx[i][j][k]=0.0;
        field_hy[i][j][k]=0.0;
        field_hz[i][j][k]=0.0;
      }
    }
  }


  //=========================================================
  // SET-UP AUXILLARY GRID
  //=========================================================
  // i. a 1D grid is set up to be a size of 2*nsteps + 20 (which is a little large)
  // but allows the wave to propagate from the center of the grid to the end

  // GET AUXILLARY GRID SIZE
  aux_npts = 2*nsteps + 20;

  
  // ALLOCATE MEMORE & INITIALIZE GRID
  ex_aux = new double [aux_npts+1];
  ey_aux = new double [aux_npts+1];
  hx_aux = new double [aux_npts+1];
  hy_aux = new double [aux_npts+1];

  for(i = 0; i <= aux_npts; ++i)
  {
    ex_aux[i] = 0.0;
    ey_aux[i] = 0.0;
    hx_aux[i] = 0.0;
    hy_aux[i] = 0.0;
  }


  // SET THE AUX TF/SF POSITION TO BE MIDDLE OF GRID
  aux_tfsf_k1 = aux_npts/2 + 1;

  // SET THE MAXWELL EQUATION COEFFICIENTS
  // i. time_dt must be set before this
  e_aux_coeff = time_dt/(EPS0*aux_epsr);
  h_aux_coeff = time_dt/MU0;

  // GET THE TRANSMISSION POINT IN THE AUXILLARY GRID
  aux_spect_transm_kpt = aux_tfsf_k1 + spect_transm_kpt - mpi_tfsf_k1;

  // !fix
  // maybe allocate auxillary ft fields, etc here

  //=========================================================
  // SET UP LOSSY MATERIALS
  //=========================================================
  // i. time_dt must be set before this is done

  // GENERATE MATERIAL LIBRARY
  generate_epsrlib();

  // ALLOCATE SOME MEMORY
  lorentz_alphap = new double*[d2l_libsize];
  lorentz_zetap = new double*[d2l_libsize];
  lorentz_gammap = new double*[d2l_libsize];

  for(i = 0; i < d2l_libsize; ++i)
  {
    lorentz_alphap[i] = new double [2];
    lorentz_zetap[i] = new double [2];
    lorentz_gammap[i] = new double [2];
  }


  //---------------------------------------------------------
  // ALOCATE MEMORY & INITIALIZE
  //---------------------------------------------------------

  // DRUDE MODEL
  drude_current_x=new double**[my_grid_nx+1+1];
  drude_current_y=new double**[my_grid_nx+1+1];
  drude_current_z=new double**[my_grid_nx+1+1];

  for(i=0;i<=my_grid_nx+1;++i)
  {
    drude_current_x[i]=new double*[my_grid_ny+1+1];
    drude_current_y[i]=new double*[my_grid_ny+1+1];
    drude_current_z[i]=new double*[my_grid_ny+1+1];

    for(j=0;j<=my_grid_ny+1;++j)
    {
      drude_current_x[i][j]=new double [my_grid_nz+1+1];
      drude_current_y[i][j]=new double [my_grid_nz+1+1];
      drude_current_z[i][j]=new double [my_grid_nz+1+1];

      // zero memory
      for(k=0;k<=my_grid_nz+1;++k)
      {
        drude_current_x[i][j][k]=0.0;
        drude_current_y[i][j][k]=0.0;
        drude_current_z[i][j][k]=0.0;
      }
    }
  }

  // LORENTZ MODEL
  lorentz_currentp_x=new double***[2];
  lorentz_currentp_y=new double***[2];
  lorentz_currentp_z=new double***[2];
  lorentz_currentp_x_nm1=new double***[2];
  lorentz_currentp_y_nm1=new double***[2];
  lorentz_currentp_z_nm1=new double***[2];

  for(l=0;l<2;++l)
  {
    lorentz_currentp_x[l]=new double**[my_grid_nx+1+1];
    lorentz_currentp_y[l]=new double**[my_grid_nx+1+1];
    lorentz_currentp_z[l]=new double**[my_grid_nx+1+1];
    lorentz_currentp_x_nm1[l]=new double**[my_grid_nx+1+1];
    lorentz_currentp_y_nm1[l]=new double**[my_grid_nx+1+1];
    lorentz_currentp_z_nm1[l]=new double**[my_grid_nx+1+1];

    for(i=0;i<=my_grid_nx+1;++i)
    {
      lorentz_currentp_x[l][i]=new double*[my_grid_ny+1+1];
      lorentz_currentp_y[l][i]=new double*[my_grid_ny+1+1];
      lorentz_currentp_z[l][i]=new double*[my_grid_ny+1+1];
      lorentz_currentp_x_nm1[l][i]=new double*[my_grid_ny+1+1];
      lorentz_currentp_y_nm1[l][i]=new double*[my_grid_ny+1+1];
      lorentz_currentp_z_nm1[l][i]=new double*[my_grid_ny+1+1];

      for(j=0;j<=my_grid_ny+1;++j)
      {
        lorentz_currentp_x[l][i][j]=new double [my_grid_nz+1+1];
        lorentz_currentp_y[l][i][j]=new double [my_grid_nz+1+1];
        lorentz_currentp_z[l][i][j]=new double [my_grid_nz+1+1];
        lorentz_currentp_x_nm1[l][i][j]=new double [my_grid_nz+1+1];
        lorentz_currentp_y_nm1[l][i][j]=new double [my_grid_nz+1+1];
        lorentz_currentp_z_nm1[l][i][j]=new double [my_grid_nz+1+1];

        // zero memory
        for(k=0;k<=my_grid_nz+1;++k)
        {
          lorentz_currentp_x[l][i][j][k]=0.0;
          lorentz_currentp_y[l][i][j][k]=0.0;
          lorentz_currentp_z[l][i][j][k]=0.0;
          lorentz_currentp_x_nm1[l][i][j][k]=0.0;
          lorentz_currentp_y_nm1[l][i][j][k]=0.0;
          lorentz_currentp_z_nm1[l][i][j][k]=0.0;
        }
      }
    }
  }


  // FIELD (METAL) UPDATE COEFFICIENTS
  drude_kp = new double [d2l_libsize];
  drude_betap = new double [d2l_libsize];

  metal_c1 = new double [d2l_libsize];
  metal_c2 = new double [d2l_libsize];
  metal_c3 = new double [d2l_libsize];


  //---------------------------------------------------------
  // CALCULATE
  //---------------------------------------------------------

  // DRUDE COEFFICIENTS
  for(m = 0; m < d2l_libsize; ++m)
  {
    drude_kp[m] = (1.0 - d2l_drude_gammap[m]*time_dt/2.0)/(1.0 + d2l_drude_gammap[m]*time_dt/2.0);
    drude_betap[m] = (d2l_drude_omegap[m]*d2l_drude_omegap[m]*EPS0*time_dt/2.0)/(1.0 + d2l_drude_gammap[m]*time_dt/2.0);
  }

  // LORENTZ COEFFICIENTS
  for(m=0;m<d2l_libsize;++m)
  {
    for(l=0;l<2;++l)
    {
      lorentz_alphap[m][l] = (2.0 - d2l_lorentz_omegap[m][l]*d2l_lorentz_omegap[m][l]*time_dt*time_dt)/(1.0 + d2l_lorentz_deltap[m][l]*time_dt);
      lorentz_zetap[m][l] = (d2l_lorentz_deltap[m][l]*time_dt - 1.0)/(1.0 + d2l_lorentz_deltap[m][l]*time_dt);
      lorentz_gammap[m][l] = (EPS0*d2l_lorentz_depsr[m][l]*d2l_lorentz_omegap[m][l]*d2l_lorentz_omegap[m][l]*time_dt*time_dt)/(1.0 + d2l_lorentz_deltap[m][l]*time_dt);
    }
  }


  // CALCULATE COEFFICIENTS
  for(m = 0; m < d2l_libsize; ++m)
  {
    // C1
    metal_c1[m] = EPS0*d2l_inf[m]/time_dt + d2l_conduct[m]/2.0 + 0.5*drude_betap[m];
    for(l=0;l<2;++l)
    {
      metal_c1[m] += (1.0/(4.0*time_dt))*lorentz_gammap[m][l];
    } 

    // C2
    metal_c2[m] = EPS0*d2l_inf[m]/time_dt - d2l_conduct[m]/2.0 - 0.5*drude_betap[m];

    // C3
    metal_c3[m] = 0.0;
    for(l=0;l<2;++l)
    {
      metal_c3[m] += (1.0/(4.0*time_dt))*lorentz_gammap[m][l];
    } 
  } // ++m


  //=========================================================
  // DETERMINE OUTPUT RESPONSIBILITIES
  //=========================================================

  //---------------------------------------------------------
  // OUTPUT DATA FOR STITCHING FILES TOGETHER
  //---------------------------------------------------------
  if(mpi_my_id == mpi_master_id)
  {
    output(2, 0, 0);
  }


  //---------------------------------------------------------
  // FT PLANE INITIALIZATION
  //---------------------------------------------------------

  // set having to ft fields to 0
  mpi_ift = 0;
  // set the number of planes we must ft to 0
  nft_planes = 0;



  //---------------------------------------------------------
  // FIGURE OUT PLANE OUTPUT RESPONSIBILITIES
  //---------------------------------------------------------

  // ALLOCATE MEMORY
  // ! these arrays are allocated to maximal size assuming that every output is FT, most of these
  // should really be set to nft_planes at the end of this assignment.  But the additional memory
  // neede here is minimal
  // !!! bad mem
  mpi_ioutput_plane = new int [output_nplanes+1];

  output_plane_local_coord = new int [output_nplanes+1];

  ft_wavenumbers = new double [output_nplanes+1];

  ft_aux_coord = new int [output_nplanes + 1];

  output_to_ft = new int [output_nplanes + 1];
  ft_to_output = new int [output_nplanes + 1];


  // max/min positions this processor handles
  double xpos_min, xpos_max, ypos_min, ypos_max, zpos_min, zpos_max;

  xpos_min = (static_cast<double>(mpi_my_grid_nx1))*grid_dx;
  xpos_max = (static_cast<double>(mpi_my_grid_nx2))*grid_dx;

  ypos_min = (static_cast<double>(mpi_my_grid_ny1))*grid_dy;
  ypos_max = (static_cast<double>(mpi_my_grid_ny2))*grid_dy;

  zpos_min = (static_cast<double>(mpi_my_grid_nz1))*grid_dz;
  zpos_max = (static_cast<double>(mpi_my_grid_nz2))*grid_dz;

  // ADJUST THE MIN/MAX POSITIONS THIS PROCESSOR HANDLES
  // i. this is necessary because the user may have set a plane, etc to occur
  // between two grid points, so doing this should catch most of these cases
  xpos_min -= grid_dx/2.0;
  xpos_max += grid_dx/1.99;

  ypos_min -= grid_dy/2.0;
  ypos_max += grid_dy/1.99;

  zpos_min -= grid_dy/2.0;
  zpos_max += grid_dy/1.99;


  // loop over the planes
  // i. types:
  //	1 == xy
  //	2 == xz
  //	3 == yz
  // ii. the local coordinate that is calculated to use for output could be off by 1
  // due to round off errors
  // ii. it should not matter where we do FT in the auxillary grid because we are 
  // taking E0*conj(E0), which will not depend on the position, so I will just set 
  // this to a point above where the TF/SF is.
  for(j = 1; j <= output_nplanes; ++j)
  {

    // ---------------------------------- xy ---------------------------------- 
    if(output_plane[j] == 1)
    {
      // if the mid-position of this plane falls within this processors responsibilities ...
      if((output_plane_midpos[j] >= zpos_min) && (output_plane_midpos[j] <= zpos_max))
      {
        // set this processor to output the plane
        mpi_ioutput_plane[j] = 1;

        // GET THE LOCAL COORDINATE TO OUTPUT AT
        output_plane_local_coord[j] = static_cast<int>((output_plane_midpos[j] - zpos_min)/grid_dz);
        // if it is the last grid point, shift local coord down
        if(output_plane_local_coord[j] > mpi_my_grid_nz2)
        {
          output_plane_local_coord[j]--;
        }

        // if this is an FT plane ...
        if(output_field[j] >= 9)
        {
          // SET FT FLAG
          mpi_ift = 1;

          // CONVERSION ARRAYS FOR OUTPUT-TO-FT AND VICE-VERSA
          output_to_ft[j] = nft_planes;
          ft_to_output[nft_planes] = j;

          // GET THE FT WAVENUMBER
          ft_wavenumbers[nft_planes] = 2.0*PI*CSPEED/ft_wavelengths[j];

          // GET THE FT POINT IN THE AUXILLARY GRID 
          ft_aux_coord[nft_planes] = aux_tfsf_k1 + 20;

          // INCREMENT NUMBER OF FT PLANES
          // i. we do this down here do our arrays start at 0
          nft_planes++;
        }
    
      }
      else
      {
        mpi_ioutput_plane[j] = 0;
      }
    }
    // ---------------------------------- xz ---------------------------------- 
    else if(output_plane[j] == 2)
    {

      // now check to see if the mid-position of this plane falls
      // within this processors responsibilities
      if((output_plane_midpos[j] >= ypos_min) && (output_plane_midpos[j] <= ypos_max))
      {
        // set this processor to output the plane
        mpi_ioutput_plane[j] = 1;

        // GET THE LOCAL COORDINATE TO OUTPUT AT
        output_plane_local_coord[j] = static_cast<int>((output_plane_midpos[j] - ypos_min)/grid_dy);
        // if it is the last grid point, shift local coord down
        if(output_plane_local_coord[j] > mpi_my_grid_ny2)
        {
          output_plane_local_coord[j]--;
        }


        // if this is an FT plane ...
        if(output_field[j] >= 9)
        {
          // SET FT FLAG
          mpi_ift = 1;

          // CONVERSION ARRAYS FOR OUTPUT-TO-FT AND VICE-VERSA
          output_to_ft[j] = nft_planes;
          ft_to_output[nft_planes] = j;

          // GET THE FT WAVENUMBER
          ft_wavenumbers[nft_planes] = 2.0*PI*CSPEED/ft_wavelengths[j];

          // GET THE FT POINT IN THE AUXILLARY GRID 
          ft_aux_coord[nft_planes] = aux_tfsf_k1 + 20;

          // INCREMENT NUMBER OF FT PLANES
          // i. we do this down here do our arrays start at 0
          nft_planes++;
        }

      }
      else
      {
        mpi_ioutput_plane[j] = 0;
      }
    }  
    // ---------------------------------- yz ---------------------------------- 
    // i.) note that we use an else if and not else in case we used a bad type value
    else if(output_plane[j] == 3)
    {
      // if the mid-position of this plane falls within this processors responsibilities ...
      if((output_plane_midpos[j] >= xpos_min) && (output_plane_midpos[j] <= xpos_max))
      {
        // set this processor to output the plane
        mpi_ioutput_plane[j] = 1;

        // GET THE LOCAL COORDINATE TO OUTPUT AT
        output_plane_local_coord[j] = static_cast<int>((output_plane_midpos[j] - xpos_min)/grid_dx);
        // if it is the last grid point, shift local coord down
        if(output_plane_local_coord[j] > mpi_my_grid_nx2)
        {
          output_plane_local_coord[j]--;
        }

        // if this is an FT plane ...
        if(output_field[j] >= 9)
        {
          // SET FT FLAG
          mpi_ift = 1;

          // CONVERSION ARRAYS FOR OUTPUT-TO-FT AND VICE-VERSA
          output_to_ft[j] = nft_planes;
          ft_to_output[nft_planes] = j;

          // GET THE FT WAVENUMBER
          ft_wavenumbers[nft_planes] = 2.0*PI*CSPEED/ft_wavelengths[j];

          // GET THE FT POINT IN THE AUXILLARY GRID 
          ft_aux_coord[nft_planes] = aux_tfsf_k1 + 20;

          // INCREMENT NUMBER OF FT PLANES
          // i. we do this down here do our arrays start at 0
          nft_planes++;
        }

      }
      else
      {
        mpi_ioutput_plane[j] = 0;
      }
    }

  } // ++j


  //---------------------------------------------------------
  // ALLOCATE MEMORY
  //---------------------------------------------------------
  if(nft_planes!=0)
  {
    aux_ftex = new complex<double> [nft_planes];
    aux_ftey = new complex<double> [nft_planes];
    aux_ftez = new complex<double> [nft_planes]; 

    aux_fthx = new complex<double> [nft_planes];
    aux_fthy = new complex<double> [nft_planes];
    aux_fthz = new complex<double> [nft_planes]; 

    field_ft_ex=new complex<double>**[nft_planes];
    field_ft_ey=new complex<double>**[nft_planes];
    field_ft_ez=new complex<double>**[nft_planes];
    field_ft_hx=new complex<double>**[nft_planes];
    field_ft_hy=new complex<double>**[nft_planes];
    field_ft_hz=new complex<double>**[nft_planes];
  }


  for(p=0;p<nft_planes;++p)
  {

    aux_ftex[p] = COMPLEXZERO;
    aux_ftey[p] = COMPLEXZERO;
    aux_ftez[p] = COMPLEXZERO;

    aux_fthx[p] = COMPLEXZERO;
    aux_fthy[p] = COMPLEXZERO;
    aux_fthz[p] = COMPLEXZERO;

    // ---------------------------------- xy ---------------------------------- 
    if(output_plane[ft_to_output[p]] == 1)
    {
      field_ft_ex[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_ey[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_ez[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hx[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hy[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hz[p]=new complex<double>*[my_grid_nx+1+1];

      for(i=0;i<=my_grid_nx+1; ++i)
      {
        field_ft_ex[p][i]=new complex<double> [my_grid_ny+1+1];
        field_ft_ey[p][i]=new complex<double> [my_grid_ny+1+1];
        field_ft_ez[p][i]=new complex<double> [my_grid_ny+1+1];
        field_ft_hx[p][i]=new complex<double> [my_grid_ny+1+1];
        field_ft_hy[p][i]=new complex<double> [my_grid_ny+1+1];
        field_ft_hz[p][i]=new complex<double> [my_grid_ny+1+1];

        // INITIALIZE MEMORY
        for(j=0;j<=my_grid_ny+1; ++j)
        {
          field_ft_ex[p][i][j]=COMPLEXZERO;
          field_ft_ey[p][i][j]=COMPLEXZERO;
          field_ft_ez[p][i][j]=COMPLEXZERO;
          field_ft_hx[p][i][j]=COMPLEXZERO;
          field_ft_hy[p][i][j]=COMPLEXZERO;
          field_ft_hz[p][i][j]=COMPLEXZERO;
        }
      }
    }
    // ---------------------------------- xz ---------------------------------- 
    else if(output_plane[ft_to_output[p]] == 2)
    {
      field_ft_ex[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_ey[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_ez[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hx[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hy[p]=new complex<double>*[my_grid_nx+1+1];
      field_ft_hz[p]=new complex<double>*[my_grid_nx+1+1];

      for(i=0;i<=my_grid_nx+1; ++i)
      {
        field_ft_ex[p][i]=new complex<double> [my_grid_nz+1+1];
        field_ft_ey[p][i]=new complex<double> [my_grid_nz+1+1];
        field_ft_ez[p][i]=new complex<double> [my_grid_nz+1+1];
        field_ft_hx[p][i]=new complex<double> [my_grid_nz+1+1];
        field_ft_hy[p][i]=new complex<double> [my_grid_nz+1+1];
        field_ft_hz[p][i]=new complex<double> [my_grid_nz+1+1];

        // INITIALIZE MEMORY
        for(k=0;k<=my_grid_nz+1; ++k)
        {
          field_ft_ex[p][i][k]=COMPLEXZERO;
          field_ft_ey[p][i][k]=COMPLEXZERO;
          field_ft_ez[p][i][k]=COMPLEXZERO;
          field_ft_hx[p][i][k]=COMPLEXZERO;
          field_ft_hy[p][i][k]=COMPLEXZERO;
          field_ft_hz[p][i][k]=COMPLEXZERO;
        }
      }
    }
    // ---------------------------------- yz ---------------------------------- 
    else if(output_plane[ft_to_output[p]] == 3)
    {
      field_ft_ex[p]=new complex<double>*[my_grid_ny+1+1];
      field_ft_ey[p]=new complex<double>*[my_grid_ny+1+1];
      field_ft_ez[p]=new complex<double>*[my_grid_ny+1+1];
      field_ft_hx[p]=new complex<double>*[my_grid_ny+1+1];
      field_ft_hy[p]=new complex<double>*[my_grid_ny+1+1];
      field_ft_hz[p]=new complex<double>*[my_grid_ny+1+1];

      for(j=0;j<=my_grid_ny+1; ++j)
      {
        field_ft_ex[p][j]=new complex<double> [my_grid_nz+1+1];
        field_ft_ey[p][j]=new complex<double> [my_grid_nz+1+1];
        field_ft_ez[p][j]=new complex<double> [my_grid_nz+1+1];
        field_ft_hx[p][j]=new complex<double> [my_grid_nz+1+1];
        field_ft_hy[p][j]=new complex<double> [my_grid_nz+1+1];
        field_ft_hz[p][j]=new complex<double> [my_grid_nz+1+1];

        // INITIALIZE MEMORY
        for(k=0;k<=my_grid_nz+1; ++k)
        {
          field_ft_ex[p][j][k]=COMPLEXZERO;
          field_ft_ey[p][j][k]=COMPLEXZERO;
          field_ft_ez[p][j][k]=COMPLEXZERO;
          field_ft_hx[p][j][k]=COMPLEXZERO;
          field_ft_hy[p][j][k]=COMPLEXZERO;
          field_ft_hz[p][j][k]=COMPLEXZERO;
        }
      }
    }

  } // ++p




  //=========================================================
  // TRANSMISSION/REFLECTION - CROSS SECTION CALCULATIONS
  //=========================================================
  // i. we only set this up if the user has set to calculate the spectra
  // so we can conserve memory

  // i. these transmission arrays are needed for the cross section calculations
  // to get the incident energy flux
  if((itransm==1) || (iscat_calc==1))
  {


    //=========================================================
    // SET UP FOURIER TRANSFORM ARRAYS (AND DFT MULTIPLIER ARRAY)
    //=========================================================
    // i. the reason we set up these arrays is that we use these at almost every time step
    // (to do any fourier transform), and these can become very expensive to calculate, but
    // we can access the memory rather quickly
    // ii. we better have set up nsteps, spect_wavenumbers[], and time_dt before doing this
  
    // ALLOCATE MEMORY
    spect_wavenumbers = new double [spect_npts+1];

    // GET THE FT WAVENUMBERS
    double spect_ft_wavelength;

    // zero all arrays
    for(n = 1; n <= spect_npts; n++)
    {
      // get the wavelength for this point
      spect_ft_wavelength = (static_cast<double>(n-1))*(spect_maxwave - spect_minwave)/(static_cast<double>(spect_npts - 1)) + spect_minwave;

      // calculate the wavevector for this wavelength
      spect_wavenumbers[n] = 2.0*PI*CSPEED/spect_ft_wavelength;
    }


    double time_temp=0.0;

    dft_multiplier_t = new complex<double>*[spect_npts+1];
    dft_multiplier_ht = new complex<double>*[spect_npts+1]; 

    for(m = 0; m <= spect_npts; ++m)
    {
      dft_multiplier_t[m] = new complex<double> [nsteps+1];
      dft_multiplier_ht[m] = new complex<double> [nsteps+1];
    }

    for(m = 1; m <= spect_npts; ++m)
    {
      for(n = 1; n <= nsteps; ++n)
      {
        // i. the first pass through the time will be 0.0, but these arrays are used after the
        // fields are updated, so we add time_dt or time_dt/2.0
        dft_multiplier_t[m][n] = exp(0.0 - COMPLEXJ*spect_wavenumbers[m]*(time_temp+time_dt));
        dft_multiplier_ht[m][n] = exp(0.0 - COMPLEXJ*spect_wavenumbers[m]*(time_temp+time_dt/2.0));

        // UPDATE THE TIME
        time_temp+=time_dt; 
      }
    }


    aux_ftex_transm = new complex<double> [spect_npts+1];
    aux_ftey_transm = new complex<double> [spect_npts+1];
    aux_ftez_transm = new complex<double> [spect_npts+1];
    aux_fthx_transm = new complex<double> [spect_npts+1];
    aux_fthy_transm = new complex<double> [spect_npts+1];
    aux_fthz_transm = new complex<double> [spect_npts+1];

    for(m=0;m<=spect_npts;++m)
    {
      aux_ftex_transm[m] = COMPLEXZERO;
      aux_ftey_transm[m] = COMPLEXZERO;
      aux_ftez_transm[m] = COMPLEXZERO;
      aux_fthx_transm[m] = COMPLEXZERO;
      aux_fthy_transm[m] = COMPLEXZERO;
      aux_fthz_transm[m] = COMPLEXZERO;
    }
  }


  if(itransm==1)
  {
    transm_ftex=new complex<double>**[spect_npts+1];
    transm_ftey=new complex<double>**[spect_npts+1];
    transm_fthx=new complex<double>**[spect_npts+1];
    transm_fthy=new complex<double>**[spect_npts+1];
    refl_ftex=new complex<double>**[spect_npts+1];
    refl_ftey=new complex<double>**[spect_npts+1];
    refl_fthx=new complex<double>**[spect_npts+1];
    refl_fthy=new complex<double>**[spect_npts+1];

    //aux_ftex_transm = new complex<double> [spect_npts+1];
    //aux_ftey_transm = new complex<double> [spect_npts+1];
    //aux_ftez_transm = new complex<double> [spect_npts+1];
    //aux_fthx_transm = new complex<double> [spect_npts+1];
    //aux_fthy_transm = new complex<double> [spect_npts+1];
    //aux_fthz_transm = new complex<double> [spect_npts+1];
    aux_ftex_refl = new complex<double> [spect_npts+1];
    aux_ftey_refl = new complex<double> [spect_npts+1];
    aux_ftez_refl = new complex<double> [spect_npts+1];
    aux_fthx_refl = new complex<double> [spect_npts+1];
    aux_fthy_refl = new complex<double> [spect_npts+1];
    aux_fthz_refl = new complex<double> [spect_npts+1];

    for(m=0;m<=spect_npts;++m)
    {
      //aux_ftex_transm[m] = COMPLEXZERO;
      //aux_ftey_transm[m] = COMPLEXZERO;
      //aux_ftez_transm[m] = COMPLEXZERO;
      //aux_fthx_transm[m] = COMPLEXZERO;
      //aux_fthy_transm[m] = COMPLEXZERO;
      //aux_fthz_transm[m] = COMPLEXZERO;
      aux_ftex_refl[m] = COMPLEXZERO;
      aux_ftey_refl[m] = COMPLEXZERO;
      aux_ftez_refl[m] = COMPLEXZERO;
      aux_fthx_refl[m] = COMPLEXZERO;
      aux_fthy_refl[m] = COMPLEXZERO;
      aux_fthz_refl[m] = COMPLEXZERO;

      transm_ftex[m]=new complex<double>*[my_grid_nx+1+1];
      transm_ftey[m]=new complex<double>*[my_grid_nx+1+1];
      transm_fthx[m]=new complex<double>*[my_grid_nx+1+1];
      transm_fthy[m]=new complex<double>*[my_grid_nx+1+1];
      refl_ftex[m]=new complex<double>*[my_grid_nx+1+1];
      refl_ftey[m]=new complex<double>*[my_grid_nx+1+1];
      refl_fthx[m]=new complex<double>*[my_grid_nx+1+1];
      refl_fthy[m]=new complex<double>*[my_grid_nx+1+1];

      for(i=0;i<=my_grid_nx+1;++i)
      {
        transm_ftex[m][i]=new complex<double> [my_grid_ny+1+1];
        transm_ftey[m][i]=new complex<double> [my_grid_ny+1+1];
        transm_fthx[m][i]=new complex<double> [my_grid_ny+1+1];
        transm_fthy[m][i]=new complex<double> [my_grid_ny+1+1];
        refl_ftex[m][i]=new complex<double> [my_grid_ny+1+1];
        refl_ftey[m][i]=new complex<double> [my_grid_ny+1+1];
        refl_fthx[m][i]=new complex<double> [my_grid_ny+1+1];
        refl_fthy[m][i]=new complex<double> [my_grid_ny+1+1];

        // zero memory
        for(j=0;j<=my_grid_ny+1;++j)
        {
          transm_ftex[m][i][j]=COMPLEXZERO;
          transm_ftey[m][i][j]=COMPLEXZERO;
          transm_fthx[m][i][j]=COMPLEXZERO;
          transm_fthy[m][i][j]=COMPLEXZERO;
          refl_ftex[m][i][j]=COMPLEXZERO;
          refl_ftey[m][i][j]=COMPLEXZERO;
          refl_fthx[m][i][j]=COMPLEXZERO;
          refl_fthy[m][i][j]=COMPLEXZERO;
        }
      }
    }


    // FIRST GET THE k POINTS FOR SPECTRA CALCULATIONS
    spect_transm_kpt = static_cast<int>(spect_transm_pos/grid_dz) + 1;
    spect_refl_kpt = static_cast<int>(spect_refl_pos/grid_dz) + 1;

    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    mpi_ft_spect_transm = 0;
    mpi_ft_spect_refl = 0;

    if((spect_transm_kpt >= mpi_my_grid_nz1) && (spect_transm_kpt <= mpi_my_grid_nz2))
    {
      mpi_ft_spect_transm = 1;
    }

    if((spect_refl_kpt >= mpi_my_grid_nz1) && (spect_refl_kpt <= mpi_my_grid_nz2))
    {
      mpi_ft_spect_refl = 1;
    }
  } // end if(itransm==1)


  //=========================================================
  // NON-PERIODIC SCATTERING CROSS SECTION CALCULATION
  //=========================================================
  // i. we only set this up if the user has set to calculate the scattering cross section
  // so we can conserve memory

  if(iscat_calc==1)
  {
    scat_k1_ftex = new complex<double>**[spect_npts+1];
    scat_k1_ftey = new complex<double>**[spect_npts+1];
    scat_k1_fthx = new complex<double>**[spect_npts+1];
    scat_k1_fthy = new complex<double>**[spect_npts+1];
    scat_k2_ftex = new complex<double>**[spect_npts+1];
    scat_k2_ftey = new complex<double>**[spect_npts+1];
    scat_k2_fthx = new complex<double>**[spect_npts+1];
    scat_k2_fthy = new complex<double>**[spect_npts+1];

    scat_j1_ftex = new complex<double>**[spect_npts+1];
    scat_j1_ftez = new complex<double>**[spect_npts+1];
    scat_j1_fthx = new complex<double>**[spect_npts+1];
    scat_j1_fthz = new complex<double>**[spect_npts+1];
    scat_j2_ftex = new complex<double>**[spect_npts+1];
    scat_j2_ftez = new complex<double>**[spect_npts+1];
    scat_j2_fthx = new complex<double>**[spect_npts+1];
    scat_j2_fthz = new complex<double>**[spect_npts+1];

    scat_i1_ftez = new complex<double>**[spect_npts+1];
    scat_i1_ftey = new complex<double>**[spect_npts+1];
    scat_i1_fthz = new complex<double>**[spect_npts+1];
    scat_i1_fthy = new complex<double>**[spect_npts+1];
    scat_i2_ftez = new complex<double>**[spect_npts+1];
    scat_i2_ftey = new complex<double>**[spect_npts+1];
    scat_i2_fthz = new complex<double>**[spect_npts+1];
    scat_i2_fthy = new complex<double>**[spect_npts+1];

    scat_aux_ex = new complex<double>*[spect_npts+1];
    scat_aux_ey = new complex<double>*[spect_npts+1];
    scat_aux_hx = new complex<double>*[spect_npts+1];
    scat_aux_hy = new complex<double>*[spect_npts+1];


    for(m = 0; m <= spect_npts; ++m)
    {

      scat_aux_ex[m] = new complex<double> [aux_npts+1];
      scat_aux_ey[m] = new complex<double> [aux_npts+1];
      scat_aux_hx[m] = new complex<double> [aux_npts+1];
      scat_aux_hy[m] = new complex<double> [aux_npts+1];


      for(i = 0; i <= aux_npts; ++i)
      {
        scat_aux_ex[m][i] = COMPLEXZERO;
        scat_aux_ey[m][i] = COMPLEXZERO;
        scat_aux_hx[m][i] = COMPLEXZERO;
        scat_aux_hy[m][i] = COMPLEXZERO;
      }


      // ALLOCATE MEMORY
      scat_k1_ftex[m] = new complex<double>*[my_grid_nx+1];
      scat_k1_ftey[m] = new complex<double>*[my_grid_nx+1];
      scat_k1_fthx[m] = new complex<double>*[my_grid_nx+1];
      scat_k1_fthy[m] = new complex<double>*[my_grid_nx+1];
      scat_k2_ftex[m] = new complex<double>*[my_grid_nx+1];
      scat_k2_ftey[m] = new complex<double>*[my_grid_nx+1];
      scat_k2_fthx[m] = new complex<double>*[my_grid_nx+1];
      scat_k2_fthy[m] = new complex<double>*[my_grid_nx+1];

      scat_j1_ftex[m] = new complex<double>*[my_grid_nx+1];
      scat_j1_ftez[m] = new complex<double>*[my_grid_nx+1];
      scat_j1_fthx[m] = new complex<double>*[my_grid_nx+1];
      scat_j1_fthz[m] = new complex<double>*[my_grid_nx+1];
      scat_j2_ftex[m] = new complex<double>*[my_grid_nx+1];
      scat_j2_ftez[m] = new complex<double>*[my_grid_nx+1];
      scat_j2_fthx[m] = new complex<double>*[my_grid_nx+1];
      scat_j2_fthz[m] = new complex<double>*[my_grid_nx+1];

      // ALLOCATE & ZERO MEMORY
      for(i = 0; i <= my_grid_nx; ++i)
      {
      
        // ALLOCATE MEMORY
        scat_k1_ftex[m][i] = new complex<double> [my_grid_ny+1];
        scat_k1_ftey[m][i] = new complex<double> [my_grid_ny+1];
        scat_k1_fthx[m][i] = new complex<double> [my_grid_ny+1];
        scat_k1_fthy[m][i] = new complex<double> [my_grid_ny+1];
        scat_k2_ftex[m][i] = new complex<double> [my_grid_ny+1];
        scat_k2_ftey[m][i] = new complex<double> [my_grid_ny+1];
        scat_k2_fthx[m][i] = new complex<double> [my_grid_ny+1];
        scat_k2_fthy[m][i] = new complex<double> [my_grid_ny+1];


        // ZERO ARRAYS
        for(j = 0; j <= my_grid_ny; ++j)
        {
          scat_k1_ftex[m][i][j] = COMPLEXZERO;
          scat_k1_ftey[m][i][j] = COMPLEXZERO;
          scat_k1_fthx[m][i][j] = COMPLEXZERO;
          scat_k1_fthy[m][i][j] = COMPLEXZERO;
          scat_k2_ftex[m][i][j] = COMPLEXZERO;
          scat_k2_ftey[m][i][j] = COMPLEXZERO;
          scat_k2_fthx[m][i][j] = COMPLEXZERO;
          scat_k2_fthy[m][i][j] = COMPLEXZERO;
        }


        // ALLOCATE MEMORY
        scat_j1_ftex[m][i] = new complex<double> [my_grid_nz+1];
        scat_j1_ftez[m][i] = new complex<double> [my_grid_nz+1];
        scat_j1_fthx[m][i] = new complex<double> [my_grid_nz+1];
        scat_j1_fthz[m][i] = new complex<double> [my_grid_nz+1];
        scat_j2_ftex[m][i] = new complex<double> [my_grid_nz+1];
        scat_j2_ftez[m][i] = new complex<double> [my_grid_nz+1];
        scat_j2_fthx[m][i] = new complex<double> [my_grid_nz+1];
        scat_j2_fthz[m][i] = new complex<double> [my_grid_nz+1];


        // ZERO ARRAYS
        for(k = 0; k <= my_grid_nz; ++k)
        {
          scat_j1_ftex[m][i][k] = COMPLEXZERO;
          scat_j1_ftez[m][i][k] = COMPLEXZERO;
          scat_j1_fthx[m][i][k] = COMPLEXZERO;
          scat_j1_fthz[m][i][k] = COMPLEXZERO;
          scat_j2_ftex[m][i][k] = COMPLEXZERO;
          scat_j2_ftez[m][i][k] = COMPLEXZERO;
          scat_j2_fthx[m][i][k] = COMPLEXZERO;
          scat_j2_fthz[m][i][k] = COMPLEXZERO;
        }
      } //++i


      scat_i1_ftez[m] = new complex<double>*[my_grid_ny+1];
      scat_i1_ftey[m] = new complex<double>*[my_grid_ny+1];
      scat_i1_fthz[m] = new complex<double>*[my_grid_ny+1];
      scat_i1_fthy[m] = new complex<double>*[my_grid_ny+1];
      scat_i2_ftez[m] = new complex<double>*[my_grid_ny+1];
      scat_i2_ftey[m] = new complex<double>*[my_grid_ny+1];
      scat_i2_fthz[m] = new complex<double>*[my_grid_ny+1];
      scat_i2_fthy[m] = new complex<double>*[my_grid_ny+1];

      // ALLOCATE MEMORY & ZERO ARRAYS
      for(j = 0; j <= my_grid_ny; ++j)
      {
        scat_i1_ftez[m][j] = new complex<double> [my_grid_nz+1];
        scat_i1_ftey[m][j] = new complex<double> [my_grid_nz+1];
        scat_i1_fthz[m][j] = new complex<double> [my_grid_nz+1];
        scat_i1_fthy[m][j] = new complex<double> [my_grid_nz+1];
        scat_i2_ftez[m][j] = new complex<double> [my_grid_nz+1];
        scat_i2_ftey[m][j] = new complex<double> [my_grid_nz+1];
        scat_i2_fthz[m][j] = new complex<double> [my_grid_nz+1];
        scat_i2_fthy[m][j] = new complex<double> [my_grid_nz+1];

        // ZERO ARRAYS
        for(k = 0; k <= my_grid_nz; ++k)
        {
          scat_i1_ftez[m][j][k] = COMPLEXZERO;
          scat_i1_ftey[m][j][k] = COMPLEXZERO;
          scat_i1_fthz[m][j][k] = COMPLEXZERO;
          scat_i1_fthy[m][j][k] = COMPLEXZERO;
          scat_i2_ftez[m][j][k] = COMPLEXZERO;
          scat_i2_ftey[m][j][k] = COMPLEXZERO;
          scat_i2_fthz[m][j][k] = COMPLEXZERO;
          scat_i2_fthy[m][j][k] = COMPLEXZERO;
        }
      } //++j

    } // ++m


    // GET BOUNDING POSITIONS FOR CS BOX
    scat_i1pt = static_cast<int>((scat_xcenter - scat_radius)/grid_dx) + 1;
    scat_i2pt = static_cast<int>((scat_xcenter + scat_radius)/grid_dx) + 1;

    scat_j1pt = static_cast<int>((scat_ycenter - scat_radius)/grid_dy) + 1;
    scat_j2pt = static_cast<int>((scat_ycenter + scat_radius)/grid_dy) + 1;

    scat_k1pt = static_cast<int>((scat_zcenter - scat_radius)/grid_dz) + 1;
    scat_k2pt = static_cast<int>((scat_zcenter + scat_radius)/grid_dz) + 1;


    // INITIALLY SET TO NO RESPONSIBLE TO FT ANY FIELDS
    scat_k1ft = 0;
    scat_k2ft = 0;
    scat_j1ft = 0;
    scat_j2ft = 0;
    scat_i1ft = 0;
    scat_i2ft = 0;


    // K1
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_k1pt >= mpi_my_grid_nz1) && (scat_k1pt <= mpi_my_grid_nz2))
    {

      scat_k1ft = 1;

      if(scat_i1pt < mpi_my_grid_nx1)
      {
        scat_k1x1 = mpi_my_grid_nx1;
      }
      else if((scat_i1pt >= mpi_my_grid_nx1) && (scat_i1pt <= mpi_my_grid_nx2))
      {
        scat_k1x1 = scat_i1pt;
      }
      else
      {
        scat_k1ft = 0;
      }

      if(scat_i2pt > mpi_my_grid_nx2)
      {
        scat_k1x2 = mpi_my_grid_nx2;
      }
      else if((scat_i2pt >= mpi_my_grid_nx1) && (scat_i2pt <= mpi_my_grid_nx2))
      {
        scat_k1x2 = scat_i2pt;
      }
      else
      {
        scat_k1ft = 0;
      }


      if(scat_j1pt < mpi_my_grid_ny1)
      {
        scat_k1y1 = mpi_my_grid_ny1;
      }
      else if((scat_j1pt >= mpi_my_grid_ny1) && (scat_j1pt <= mpi_my_grid_ny2))
      {
        scat_k1y1 = scat_j1pt;
      }
      else
      {
        scat_k1ft = 0;
      }

      if(scat_j2pt > mpi_my_grid_ny2)
      {
        scat_k1y2 = mpi_my_grid_ny2;
      }
      else if((scat_j2pt >= mpi_my_grid_ny1) && (scat_j2pt <= mpi_my_grid_ny2))
      {
        scat_k1y2 = scat_j2pt;
      }
      else
      {
        scat_k1ft = 0;
      }

    }


    // K2
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_k2pt >= mpi_my_grid_nz1) && (scat_k2pt <= mpi_my_grid_nz2))
    {

      scat_k2ft = 1;

      if(scat_i1pt < mpi_my_grid_nx1)
      {
        scat_k2x1 = mpi_my_grid_nx1;
      }
      else if((scat_i1pt >= mpi_my_grid_nx1) && (scat_i1pt <= mpi_my_grid_nx2))
      {
        scat_k2x1 = scat_i1pt;
      }
      else
      {
        scat_k2ft = 0;
      }
 
      if(scat_i2pt > mpi_my_grid_nx2)
      {
        scat_k2x2 = mpi_my_grid_nx2;
      }
      else if((scat_i2pt >= mpi_my_grid_nx1) && (scat_i2pt <= mpi_my_grid_nx2))
      {
        scat_k2x2 = scat_i2pt;
      }
      else
      {
        scat_k2ft = 0;
      }


      if(scat_j1pt < mpi_my_grid_ny1)
      {
        scat_k2y1 = mpi_my_grid_ny1;
      }
      else if((scat_j1pt >= mpi_my_grid_ny1) && (scat_j1pt <= mpi_my_grid_ny2))
      {
        scat_k2y1 = scat_j1pt;
      }
      else
      {
        scat_k2ft = 0;
      }

      if(scat_j2pt > mpi_my_grid_ny2)
      {
        scat_k2y2 = mpi_my_grid_ny2;
      }
      else if((scat_j2pt >= mpi_my_grid_ny1) && (scat_j2pt <= mpi_my_grid_ny2))
      {
        scat_k2y2 = scat_j2pt;
      }
      else
      {
        scat_k2ft = 0;
      }
    }

    // J1
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_j1pt >= mpi_my_grid_ny1) && (scat_j1pt <= mpi_my_grid_ny2))
    {

      scat_j1ft = 1;

      if(scat_i1pt < mpi_my_grid_nx1)
      {
        scat_j1x1 = mpi_my_grid_nx1;
      }
      else if((scat_i1pt >= mpi_my_grid_nx1) && (scat_i1pt <= mpi_my_grid_nx2))
      {
        scat_j1x1 = scat_i1pt;
      }
      else
      {
        scat_j1ft = 0;
      }

      if(scat_i2pt > mpi_my_grid_nx2)
      {
        scat_j1x2 = mpi_my_grid_nx2;
      }
      else if((scat_i2pt >= mpi_my_grid_nx1) && (scat_i2pt <= mpi_my_grid_nx2))
      {
        scat_j1x2 = scat_i2pt;
      }
      else
      {
        scat_j1ft = 0;
      }


      if(scat_k1pt < mpi_my_grid_nz1)
      {
        scat_j1z1 = mpi_my_grid_nz1;
      }
      else if((scat_k1pt >= mpi_my_grid_nz1) && (scat_k1pt <= mpi_my_grid_nz2))
      {
        scat_j1z1 = scat_k1pt;
      }
      else
      {
        scat_j1ft = 0;
      }

      if(scat_k2pt > mpi_my_grid_nz2)
      {
        scat_j1z2 = mpi_my_grid_nz2;
      }
      else if((scat_k2pt >= mpi_my_grid_nz1) && (scat_k2pt <= mpi_my_grid_nz2))
      {
        scat_j1z2 = scat_k2pt;
      }
      else
      {
        scat_j1ft = 0;
      }

    }


    // J2
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_j2pt >= mpi_my_grid_ny1) && (scat_j2pt <= mpi_my_grid_ny2))
    {

      scat_j2ft = 1;

      if(scat_i1pt < mpi_my_grid_nx1)
      {
        scat_j2x1 = mpi_my_grid_nx1;
      }
      else if((scat_i1pt >= mpi_my_grid_nx1) && (scat_i1pt <= mpi_my_grid_nx2))
      {
        scat_j2x1 = scat_i1pt;
      }
      else
      {
        scat_j2ft = 0;
      }

      if(scat_i2pt > mpi_my_grid_nx2)
      {
        scat_j2x2 = mpi_my_grid_nx2;
      }
      else if((scat_i2pt >= mpi_my_grid_nx1) && (scat_i2pt <= mpi_my_grid_nx2))
      {
        scat_j2x2 = scat_i2pt;
      }
      else
      {
        scat_j2ft = 0;
      }


      if(scat_k1pt < mpi_my_grid_nz1)
      {
        scat_j2z1 = mpi_my_grid_nz1;
      }
      else if((scat_k1pt >= mpi_my_grid_nz1) && (scat_k1pt <= mpi_my_grid_nz2))
      {
        scat_j2z1 = scat_k1pt;
      }
      else
      {
        scat_j2ft = 0;
      }

      if(scat_k2pt > mpi_my_grid_nz2)
      {
        scat_j2z2 = mpi_my_grid_nz2;
      }
      else if((scat_k2pt >= mpi_my_grid_nz1) && (scat_k2pt <= mpi_my_grid_nz2))
      {
        scat_j2z2 = scat_k2pt;
      }
      else
      {
        scat_j2ft = 0;
      }

    }


    // I1
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_i1pt >= mpi_my_grid_nx1) && (scat_i1pt <= mpi_my_grid_nx2))
    {

      scat_i1ft = 1;

      if(scat_j1pt < mpi_my_grid_ny1)
      {
        scat_i1y1 = mpi_my_grid_ny1;
      }
      else if((scat_j1pt >= mpi_my_grid_ny1) && (scat_j1pt <= mpi_my_grid_ny2))
      {
        scat_i1y1 = scat_j1pt;
      }
      else
      {
        scat_i1ft = 0;
      }

      if(scat_j2pt > mpi_my_grid_ny2)
      {
        scat_i1y2 = mpi_my_grid_ny2;
      }
      else if((scat_j2pt >= mpi_my_grid_ny1) && (scat_j2pt <= mpi_my_grid_ny2))
      {
        scat_i1y2 = scat_j2pt;
      }
      else
      {
        scat_i1ft = 0;
      }


      if(scat_k1pt < mpi_my_grid_nz1)
      {
        scat_i1z1 = mpi_my_grid_nz1;
      }
      else if((scat_k1pt >= mpi_my_grid_nz1) && (scat_k1pt <= mpi_my_grid_nz2))
      {
        scat_i1z1 = scat_k1pt;
      }
      else
      {
        scat_i1ft = 0;
      }

      if(scat_k2pt > mpi_my_grid_nz2)
      {
        scat_i1z2 = mpi_my_grid_nz2;
      }
      else if((scat_k2pt >= mpi_my_grid_nz1) && (scat_k2pt <= mpi_my_grid_nz2))
      {
        scat_i1z2 = scat_k2pt;
      }
      else
      {
        scat_i1ft = 0;
      }

    }


    // I2
    // NOW FIGURE OUT IF THESE k POINT IS THIS PROCESSORS RESPONSIBILITY
    if((scat_i2pt >= mpi_my_grid_nx1) && (scat_i2pt <= mpi_my_grid_nx2))
    {

      scat_i2ft = 1;

      if(scat_j1pt < mpi_my_grid_ny1)
      {
        scat_i2y1 = mpi_my_grid_ny1;
      }
      else if((scat_j1pt >= mpi_my_grid_ny1) && (scat_j1pt <= mpi_my_grid_ny2))
      {
        scat_i2y1 = scat_j1pt;
      }
      else
      {
        scat_i2ft = 0;
      }

      if(scat_j2pt > mpi_my_grid_ny2)
      {
        scat_i2y2 = mpi_my_grid_ny2;
      }
      else if((scat_j2pt >= mpi_my_grid_ny1) && (scat_j2pt <= mpi_my_grid_ny2))
      {
        scat_i2y2 = scat_j2pt;
      }
      else
      {
        scat_i2ft = 0;
      }


      if(scat_k1pt < mpi_my_grid_nz1)
      {
        scat_i2z1 = mpi_my_grid_nz1;
      }
      else if((scat_k1pt >= mpi_my_grid_nz1) && (scat_k1pt <= mpi_my_grid_nz2))
      {
        scat_i2z1 = scat_k1pt;
      }
      else
      {
        scat_i2ft = 0;
      }

      if(scat_k2pt > mpi_my_grid_nz2)
      {
        scat_i2z2 = mpi_my_grid_nz2;
      }
      else if((scat_k2pt >= mpi_my_grid_nz1) && (scat_k2pt <= mpi_my_grid_nz2))
      {
        scat_i2z2 = scat_k2pt;
      }
      else
      {
        scat_i2ft = 0;
      }

    }

  } // end (iscat_calc==1);


  //=========================================================
  // SIMFILE OUTPUT
  //=========================================================

  if(isimfile==1)
  {
    // WRITE TO SIMULATION FILE
    simfile << "GRID RESPONSIBILITIES: " << endl;
    simfile << " -mpi_my_grid_nx1: " << mpi_my_grid_nx1 << endl;
    simfile << " -mpi_my_grid_nx2: " << mpi_my_grid_nx2 << endl;
    simfile << " -mpi_my_grid_ny1: " << mpi_my_grid_ny1 << endl;
    simfile << " -mpi_my_grid_ny2: " << mpi_my_grid_ny2 << endl;
    simfile << " -mpi_my_grid_nz1: " << mpi_my_grid_nz1 << endl;
    simfile << " -mpi_my_grid_nz2: " << mpi_my_grid_nz2 << endl << endl;
 
    simfile << "TF/SF RESPONSIBILITIES: " << endl;
    simfile << " -mpi_itfsf_k1: " << mpi_itfsf_k1 << endl;
    simfile << " -mpi_tfsf_k1: " << mpi_tfsf_k1 << endl << endl;

    simfile << "GRID RESPONSIBILITIES (TRUE) " << endl;
    simfile << " -xpos_min: " << xpos_min << endl;
    simfile << " -xpos_max: " << xpos_max << endl;
    simfile << " -ypos_min: " << ypos_min << endl;
    simfile << " -ypos_max: " << ypos_max << endl;
    simfile << " -zpos_min: " << zpos_min << endl;
    simfile << " -zpos_max: " << zpos_max << endl << endl;

    simfile << "GRID RESPONSIBILITIES (ADJUSTED FOR TOLERANCES)" << endl;
    simfile << " -xpos_min: " << xpos_min << endl;
    simfile << " -xpos_max: " << xpos_max << endl;
    simfile << " -ypos_min: " << ypos_min << endl;
    simfile << " -ypos_max: " << ypos_max << endl;
    simfile << " -zpos_min: " << zpos_min << endl;
    simfile << " -zpos_max: " << zpos_max << endl << endl;

    // WRITE TO TF/SF FILE
    simfile << "PLANE RESPONSIBILITIES: " << endl;

    for(p = 1; p <= output_nplanes; ++p)
    {
      simfile << " -plane: " << p << endl;
      simfile << "   *mpi_ioutput_plane[" << p << "]: " << mpi_ioutput_plane[p] << endl;
    }

    simfile << endl;

    // WRITE TO TF/SF FILE
    simfile << "SPECTRUM RESPONSIBILITIES: " << endl;
    simfile << " -mpi_ft_spect_transm: " << mpi_ft_spect_transm << endl;
    simfile << " -mpi_ft_spect_refl: " << mpi_ft_spect_refl << endl << endl;
  }


  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void init_pml()
//	DESC:	Set's up the CMPL.
//
//	NOTES: 	i. PML is assigned to all sides of the computational domain, but if
//		there is x or y periodicity then the parameters are chosen so that it 
//		has no effect 
//
// ?? xdepth is what is defined in Taflove pg. 292
//
//========================================================================
//========================================================================
void init_pml()
{

  // local indices
  int i, ii, j, jj, k, kk;

  double xpos, ypos, zpos;

  // PML PARAMETERS
  double cpml_kappa, cpml_sigma, cpml_alpha, cpml_depth, depth, cpml_sigmamax, cpml_sigmamax_top, cpml_sigmamax_bottom;


  //=========================================================
  // ALOCATE MEMORY & INITIALIZE
  //=========================================================
  // i. the coefficients are initialized to free space

  psi_eyx = new double**[my_grid_nx+2];
  psi_ezx = new double**[my_grid_nx+2];
  psi_hyx = new double**[my_grid_nx+2];
  psi_hzx = new double**[my_grid_nx+2];
  psi_exy = new double**[my_grid_nx+2];
  psi_ezy = new double**[my_grid_nx+2];
  psi_hxy = new double**[my_grid_nx+2];
  psi_hzy = new double**[my_grid_nx+2];
  psi_exz = new double**[my_grid_nx+2];
  psi_eyz = new double**[my_grid_nx+2];
  psi_hxz = new double**[my_grid_nx+2];
  psi_hyz = new double**[my_grid_nx+2];

  kedx = new double [my_grid_nx+2];
  khdx = new double [my_grid_nx+2];
  kedy = new double [my_grid_ny+2];
  khdy = new double [my_grid_ny+2];
  kedz = new double [my_grid_nz+2];
  khdz = new double [my_grid_nz+2];

  be_x = new double [my_grid_nx+2];
  ce_x = new double [my_grid_nx+2];
  bh_x = new double [my_grid_nx+2];
  ch_x = new double [my_grid_nx+2];
  be_y = new double [my_grid_ny+2];
  ce_y = new double [my_grid_ny+2];
  bh_y = new double [my_grid_ny+2];
  ch_y = new double [my_grid_ny+2];
  be_z = new double [my_grid_nz+2];
  ce_z = new double [my_grid_nz+2];
  bh_z = new double [my_grid_nz+2];
  ch_z = new double [my_grid_nz+2];

  for(i = 0; i <= (my_grid_nx+1); ++i)
  {
    psi_eyx[i] = new double*[my_grid_ny+2];
    psi_ezx[i] = new double*[my_grid_ny+2];
    psi_hyx[i] = new double*[my_grid_ny+2];
    psi_hzx[i] = new double*[my_grid_ny+2];
    psi_exy[i] = new double*[my_grid_ny+2];
    psi_ezy[i] = new double*[my_grid_ny+2];
    psi_hxy[i] = new double*[my_grid_ny+2];
    psi_hzy[i] = new double*[my_grid_ny+2];
    psi_exz[i] = new double*[my_grid_ny+2];
    psi_eyz[i] = new double*[my_grid_ny+2];
    psi_hxz[i] = new double*[my_grid_ny+2];
    psi_hyz[i] = new double*[my_grid_ny+2];

    for(j = 0; j <= (my_grid_ny+1); ++j)
    {
      psi_eyx[i][j] = new double [my_grid_nz+2];
      psi_ezx[i][j] = new double [my_grid_nz+2];
      psi_hyx[i][j] = new double [my_grid_nz+2];
      psi_hzx[i][j] = new double [my_grid_nz+2];
      psi_exy[i][j] = new double [my_grid_nz+2];
      psi_ezy[i][j] = new double [my_grid_nz+2];
      psi_hxy[i][j] = new double [my_grid_nz+2];
      psi_hzy[i][j] = new double [my_grid_nz+2];
      psi_exz[i][j] = new double [my_grid_nz+2];
      psi_eyz[i][j] = new double [my_grid_nz+2];
      psi_hxz[i][j] = new double [my_grid_nz+2];
      psi_hyz[i][j] = new double [my_grid_nz+2];

      for(k = 0; k <= (my_grid_nz+1); ++k)
      {
        psi_eyx[i][j][k] = 0.0;
        psi_ezx[i][j][k] = 0.0;
        psi_hyx[i][j][k] = 0.0;
        psi_hzx[i][j][k] = 0.0;
        psi_exy[i][j][k] = 0.0;
        psi_ezy[i][j][k] = 0.0;
        psi_hxy[i][j][k] = 0.0;
        psi_hzy[i][j][k] = 0.0;
        psi_exz[i][j][k] = 0.0;
        psi_eyz[i][j][k] = 0.0;
        psi_hxz[i][j][k] = 0.0;
        psi_hyz[i][j][k] = 0.0;      

        be_z[k] = 0.0;
        ce_z[k] = 0.0;
        bh_z[k] = 0.0;
        ch_z[k] = 0.0;

        kedz[k] = 1.0;
        khdz[k] = 1.0;
      }

      be_y[j] = 0.0;
      ce_y[j] = 0.0;
      bh_y[j] = 0.0;
      ch_y[j] = 0.0;

      kedy[j] = 1.0;
      khdy[j] = 1.0;
    } // ++j

    be_x[i] = 0.0;
    ce_x[i] = 0.0;
    bh_x[i] = 0.0;
    ch_x[i] = 0.0;

    kedx[i] = 1.0;
    khdx[i] = 1.0;
  } // ++i 


  //=========================================================
  // CALCULATE COEFFICIENTS
  //=========================================================
  // i. the sigmamax equation is taken from Taflove's book pg. 294

  cpml_sigmamax_bottom = (cpml_sigmamax_coeff*(static_cast<double>(cpml_m + 1)))/(Z0*grid_dx*sqrt(cpml_epsr*cpml_mur)); 
  cpml_sigmamax_top = (cpml_sigmamax_coeff*(static_cast<double>(cpml_m + 1)))/(Z0*grid_dx*sqrt(cpml_epsr*cpml_mur)); 

  //---------------------------------------------------------
  // z-DIRECTION
  //---------------------------------------------------------

  // GET REAL DEPTH OF PML
  cpml_depth = (static_cast<double>(cpml_layers))*grid_dz;

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    // get local k coordinate
    k = kk - mpi_my_grid_nz1 + 1;

    //---------------------------------------------------------
    // E-FIELD
    //---------------------------------------------------------
 
    // get the zpos of this point
    zpos = (static_cast<double>(kk))*grid_dz;

    if(zpos < cpml_depth)
    {
      // get the depth into the pml layer
      depth = cpml_depth - zpos;

      cpml_sigma = cpml_sigmamax_bottom*pow((depth/cpml_depth), cpml_m);
      cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
      cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

      kedz[k] = cpml_kappa*grid_dz;

      be_z[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
      ce_z[k] = cpml_sigma*(be_z[k] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);

    }
    else if(zpos > (grid_zsize - cpml_depth))
    {

      // get the depth into the pml layer
      depth = zpos - (grid_zsize - cpml_depth);

      cpml_sigma = cpml_sigmamax_top*pow((depth/cpml_depth), cpml_m);
      cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
      cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

      kedz[k] = cpml_kappa*grid_dz;

      be_z[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
      ce_z[k] = cpml_sigma*(be_z[k] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
    }
    else
    {
      cpml_kappa = 1.0;
      cpml_sigma = 0.0;
      cpml_alpha = 0.0;

      kedz[k] = cpml_kappa*grid_dz;
      be_z[k] = 0.0;
      ce_z[k] = 0.0;
    }

    //---------------------------------------------------------
    // H-FIELD
    //---------------------------------------------------------
    zpos = grid_dz*(static_cast<double>(kk)) + grid_dz/2.0;

    if(zpos < cpml_depth)
    {
      // get the depth into the pml layer
      depth = cpml_depth - zpos;

      cpml_sigma = cpml_sigmamax_bottom*pow((depth/cpml_depth), cpml_m);
      cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
      cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

      khdz[k] = cpml_kappa*grid_dz;

      bh_z[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
      ch_z[k] = cpml_sigma*(bh_z[k] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
    }
    else if(zpos > (grid_zsize - cpml_depth))
    {
      // get the depth into the pml layer
      depth = zpos - (grid_zsize - cpml_depth);

      cpml_sigma = cpml_sigmamax_top*pow((depth/cpml_depth), cpml_m);
      cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
      cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

      khdz[k] = cpml_kappa*grid_dz;

      bh_z[k] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
      ch_z[k] = cpml_sigma*(bh_z[k] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
    }
    else
    {
      cpml_kappa = 1.0;
      cpml_sigma = 0.0;
      cpml_alpha = 0.0;

      khdz[k] = cpml_kappa*grid_dz;
      bh_z[k] = 0.0;
      ch_z[k] = 0.0;
    }
  }


  //---------------------------------------------------------
  // x-DIRECTION
  //---------------------------------------------------------


  // GET REAL DEPTH OF PML
  cpml_depth = (static_cast<double>(cpml_layers))*grid_dx;

  // SIGMA_MAX
  // i. see Taflove pg. 294
  cpml_sigmamax = (0.8*(static_cast<double>(cpml_m + 1)))/(Z0*grid_dx*sqrt(cpml_epsr*cpml_mur));  


  for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
  {
    // get local i coordinate
    i = ii - mpi_my_grid_nx1 + 1;

    // ---------------------- x-non-periodic -----------------------
    if(ixperiodic == 0)
    {

      //---------------------------------------------------------
      // E - FIELD
      //---------------------------------------------------------
      // get the zpos of this point
      xpos = (static_cast<double>(ii))*grid_dx;

      if(xpos < cpml_depth)
      {
        // get the depth into the pml layer
        depth = cpml_depth - xpos;

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        kedx[i] = cpml_kappa*grid_dx;

        be_x[i] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ce_x[i] = cpml_sigma*(be_x[i] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);

      }
      else if(xpos > (grid_xsize - cpml_depth))
      {
        // get the depth into the pml layer
        depth = xpos - (grid_xsize - cpml_depth);

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        kedx[i] = cpml_kappa*grid_dx;

        be_x[i] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ce_x[i] = cpml_sigma*(be_x[i] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else
      {
        cpml_kappa = 1.0;
        cpml_sigma = 0.0;
        cpml_alpha = 0.0;

        kedx[i] = cpml_kappa*grid_dx;
        be_x[i] = 0.0;
        ce_x[i] = 0.0;
      }

      //---------------------------------------------------------
      // H - FIELD
      //---------------------------------------------------------
      xpos = grid_dx*(static_cast<double>(ii)) + grid_dx/2.0;

      if(xpos < cpml_depth)
      {
        // get the depth into the pml layer
        depth = cpml_depth - xpos;

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        khdx[i] = cpml_kappa*grid_dx;

        bh_x[i] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ch_x[i] = cpml_sigma*(bh_x[i] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else if(xpos > (grid_xsize - cpml_depth))
      {
        // get the depth into the pml layer
        depth = xpos - (grid_xsize - cpml_depth);

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        khdx[i] = cpml_kappa*grid_dx;

        bh_x[i] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ch_x[i] = cpml_sigma*(bh_x[i] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else
      {
        cpml_kappa = 1.0;
        cpml_sigma = 0.0;
        cpml_alpha = 0.0;

        khdx[i] = cpml_kappa*grid_dx;
        bh_x[i] = 0.0;
        ch_x[i] = 0.0;
      }
    }
    // ---------------------- x-periodic -----------------------
    else
    {

      cpml_kappa = 1.0;
      cpml_sigma = 0.0;
      cpml_alpha = 0.0;


      kedx[i] = cpml_kappa*grid_dx;
      be_x[i] = 0.0;
      ce_x[i] = 0.0;

      khdx[i] = cpml_kappa*grid_dx;
      bh_x[i] = 0.0;
      ch_x[i] = 0.0;
    }

  } // ++ii



  //---------------------------------------------------------
  // y-DIRECTION
  //---------------------------------------------------------

  // GET REAL DEPTH OF PML
  cpml_depth = (static_cast<double>(cpml_layers))*grid_dy;

  // SIGMA_MAX
  // i. see Taflove pg. 294
  cpml_sigmamax = (0.8*(static_cast<double>(cpml_m + 1)))/(Z0*grid_dy*sqrt(cpml_epsr*cpml_mur));  


  for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
  {
    // get local i coordinate
    j = jj - mpi_my_grid_ny1 + 1;

    // ---------------------- y-non-periodic -----------------------
    if(iyperiodic == 0)
    {
      //---------------------------------------------------------
      // E - FIELD
      //---------------------------------------------------------
      // get the ypos of this point
      ypos = (static_cast<double>(jj))*grid_dy;

      if(ypos < cpml_depth)
      {
        // get the depth into the pml layer
        depth = cpml_depth - ypos;

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        kedy[j] = cpml_kappa*grid_dy;

        be_y[j] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ce_y[j] = cpml_sigma*(be_y[j] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);

      }
      else if(ypos > (grid_ysize - cpml_depth))
      {
        // get the depth into the pml layer
        depth = ypos - (grid_ysize - cpml_depth);

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        kedy[j] = cpml_kappa*grid_dy;

        be_y[j] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ce_y[j] = cpml_sigma*(be_y[j] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else
      {
        cpml_kappa = 1.0;
        cpml_sigma = 0.0;
        cpml_alpha = 0.0;

        kedy[j] = cpml_kappa*grid_dy;
        be_y[j] = 0.0;
        ce_y[j] = 0.0;
      }

      //---------------------------------------------------------
      // H - FIELD
      //---------------------------------------------------------
      // get the ypos of this point
      ypos = (static_cast<double>(jj))*grid_dy + grid_dy/2.0;

      if(ypos < cpml_depth)
      {
        // get the depth into the pml layer
        depth = cpml_depth - ypos;

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        khdy[j] = cpml_kappa*grid_dy;

        bh_y[j] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ch_y[j] = cpml_sigma*(bh_y[j] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else if(ypos > (grid_ysize - cpml_depth))
      {
        // get the depth into the pml layer
        depth = ypos - (grid_ysize - cpml_depth);

        cpml_sigma = cpml_sigmamax*pow((depth/cpml_depth), cpml_m);
        cpml_kappa = 1.0 + (cpml_kappamax - 1.0)*pow((depth/cpml_depth), cpml_m);
        cpml_alpha = cpml_alphamax*pow(((cpml_depth-depth)/cpml_depth), cpml_ma);

        khdy[j] = cpml_kappa*grid_dy;

        bh_y[j] = exp(0.0 - (cpml_sigma/cpml_kappa + cpml_alpha)*time_dt/EPS0);
        ch_y[j] = cpml_sigma*(bh_y[j] - 1.0)/(cpml_sigma*cpml_kappa + cpml_kappa*cpml_kappa*cpml_alpha);
      }
      else
      {
        cpml_kappa = 1.0;
        cpml_sigma = 0.0;
        cpml_alpha = 0.0;

        khdy[j] = cpml_kappa*grid_dy;
        bh_y[j] = 0.0;
        ch_y[j] = 0.0;
      }

    }
    // ---------------------- y-periodic -----------------------
    else
    {
      cpml_kappa = 1.0;
      cpml_sigma = 0.0;
      cpml_alpha = 0.0;

      kedy[j] = cpml_kappa*grid_dy;
      be_y[j] = 0.0;
      ce_y[j] = 0.0;

      khdy[j] = cpml_kappa*grid_dy;
      bh_y[j] = 0.0;
      ch_y[j] = 0.0;
    }

  } // ++ii


  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void structure_initialize()
//	DESC:	set up structure and PML
//
//	NOTES: 	i. this is called by all processes
//		ii. note the staggering again:
//
//			Ex (i+1/2, j, k)
//			Ey (i, j+1/2, k)
//			Ez (i, j, k+1/2)
//
//========================================================================
//========================================================================
void structure_initialize()
{

  // "local" indices
  int i, ii, j, jj, k, kk;

  // generic positions and distances
  double xpos, ypos, zpos;


  //=========================================================
  // INITIALIZATION ALL MATERIALS TO FREE SPACE
  //=========================================================

  for(k = 0; k <= my_grid_nz+1; ++k)
  {
    for(j = 0; j <= my_grid_ny+1; ++j)
    {
      for(i = 0; i <= my_grid_nx+1; ++i)
      {
        grid_material_ex[i][j][k] = grid_material_ey[i][j][k] = grid_material_ez[i][j][k] = 99;

        permittivity_ex[i][j][k] = permittivity_ey[i][j][k] = permittivity_ez[i][j][k] = EPS0;
        permeability_hx[i][j][k] = permeability_hy[i][j][k] = permeability_hz[i][j][k] = MU0;
      }
    }
  }


  //=========================================================
  // FILL COMPUTATIONAL GRID WITH MATERIALS AND SET
  //	FIELD COEFFICIENTS
  //=========================================================

  double tempd;

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        // Ex (i+1/2, j, k)
        xpos = (static_cast<double>(ii))*grid_dx + grid_dx/2.0;
        ypos = (static_cast<double>(jj))*grid_dy;
        zpos = (static_cast<double>(kk))*grid_dz;

        get_structure(xpos, ypos, zpos, grid_material_ex[i][j][k], permittivity_ex[i][j][k], tempd);

        // Ey (i, j+1/2, k)
        xpos = (static_cast<double>(ii))*grid_dx;
        ypos = (static_cast<double>(jj))*grid_dy + grid_dy/2.0;
        zpos = (static_cast<double>(kk))*grid_dz;

        get_structure(xpos, ypos, zpos, grid_material_ey[i][j][k], permittivity_ey[i][j][k], tempd);

        // Ez (i, j, k+1/2)
        xpos = (static_cast<double>(ii))*grid_dx;
        ypos = (static_cast<double>(jj))*grid_dy;
        zpos = (static_cast<double>(kk))*grid_dz + grid_dz/2.0;

        get_structure(xpos, ypos, zpos, grid_material_ez[i][j][k], permittivity_ez[i][j][k], tempd);

        // Hx (i, j+1/2, k+1/2)
        //xpos = (static_cast<double>(ii))*grid_dx;
        //ypos = (static_cast<double>(jj))*grid_dy + grid_dy/2.0;
        //zpos = (static_cast<double>(kk))*grid_dz + grid_dz/2.0;

        // Hy (i+1/2, j, k+1/2)
        //xpos = (static_cast<double>(ii))*grid_dx + grid_dx/2.0;
        //ypos = (static_cast<double>(jj))*grid_dy;
        //zpos = (static_cast<double>(kk))*grid_dz + grid_dz/2.0;

        // Hz (i+1/2, j+1/2, k)
        //xpos = (static_cast<double>(ii))*grid_dx + grid_dx/2.0;
        //ypos = (static_cast<double>(jj))*grid_dy + grid_dy/2.0;
        //zpos = (static_cast<double>(kk))*grid_dz;

      } // ++ii
    } // ++jj
  } // ++kk


  //=========================================================
  // GIVE SIMFILE OUTPUT
  //=========================================================
  if(isimfile==1)
  {
    simfile << " * structure initialized " << endl << endl;
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void perform_fdtd()
//
//	DESC:	performs the actual time stepping
//
//	NOTES: 	i. this is called by all processes
//		ii.
//
//
//========================================================================
//========================================================================
void perform_fdtd()
{

  // "local" indices
  int j;

  //---------------------------------------------------------
  // SET UP SOME OUTPUT STUFF
  //---------------------------------------------------------
  int iunit = 0;
  int ilook = 0; 
  int noutput_steps;

  if(noutput == 0)
  {
    noutput_steps = nsteps + 1;
  }
  else
  { 
    noutput_steps = nsteps/noutput;
  }


  //=========================================================
  // PERFORM TIME-STEPPING
  //=========================================================
  for(istep = 1; istep <= nsteps; ++istep)
  {

    //---------------------------------------------------------
    // UPDATE FIELDS
    //---------------------------------------------------------
    // i.) the order of these updates is important:
    //	i.a) all grids are being updated with time_current = (istep-1)*time_dt
    //	i.b) the e field updates must come before the h field updates because 
    //  of the time staggering (e.g. e goes from t-1/2 to t+1/2 and then h is 
    //	updated with the e t+1/2 information)
    //	i.c) it does not matter if the auxillary or "real" grid updates come 
    //	first because the e relies on h, etc 


    // !!! we start at time_current = 0
    // e should be updated from n-1/2 -> n+1/2 (n-1/2 does not exist)
    // h should be then be updated from n -> n+1
    // the auxillary grids should be updated last because they are used as previous time steps
    // to update h and e


    // update the "true" grid
    field_e_update();

    // update auxillary grid
    field_e_aux_update();

    // update the "true" grid
    field_h_update();

    // update the auxillary grid
    field_h_aux_update();

    // FT FIELDS FOR TRANSMITTANCE & REFLECTION IF WE ARE SUPPOSED TO
    if(itransm == 1)
    {
      if((mpi_ft_spect_transm == 1) || (mpi_ft_spect_refl == 1))
      {
        ft_fields_spect();
      }
    }

    // SCATTERING FT
    if(iscat_calc == 1)
    {
      ft_fields_scat();
    }

    // NSOM FT
    if(nsom_ift == 1)
    {
      ft_field_nsom();
    }

    //---------------------------------------------------------
    // OUTPUT
    //---------------------------------------------------------

    // FT FIELDS FOR OUTPUT IF WE ARE SUPPOSED TO
    if(mpi_ift == 1)
    {
      ft_fields();
    }

    // INCREMENT NUMBER OF STEPS SINCE OUTPUT
    ilook++;

    if(ilook == noutput_steps)
    {
      iunit++;

      for(j = 1; j <= output_nplanes; ++j)
      {
        if(mpi_ioutput_plane[j] == 1)
        {
          output(1, j, iunit);
        }     
      }

      // RESET STEPS SINCE OUTPUT
      ilook = 0;
    }


    //---------------------------------------------------------
    // UPDATE THE TIME 
    //---------------------------------------------------------
    // i. we do this at the end, and the first pass through the time is 0.0.  In 
    // reality the time is advanced immediately after updating the fields.  But the 
    // subroutines (Fourier transforming, etc) take this into account.
    time_current+=time_dt;  

  }


  return;

}





//========================================================================
//========================================================================
//
//	NAME:	void field_e_update()
//	DESC:	Updates the E-field.
//
//	NOTES: 	i.
//
//
//========================================================================
//========================================================================
void field_e_update()
{

  // local indices
  int i, ii, j, jj, k, kk, l;

  int material;

  double metal_update_term, e_n, *j_n;
  j_n = new double [2];


  //=========================================================
  // UPDATE ENTIRE GRID
  //=========================================================

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        //=========================================================
        // Ex
        //=========================================================

        // SET THE MATERIAL
        material = grid_material_ex[i][j][k];
 
        // ---------------------------- LOSSY MATERIAL ---------------------------- 
        if(material != 99)
        {

          //---------------------------------------------------------
          // STORE E & J AT (t-1/2)
          //---------------------------------------------------------
          e_n = field_ex[i][j][k];

          for(l=0;l<2;++l)
          {
            j_n[l] = lorentz_currentp_x[l][i][j][k];
          }   

          //---------------------------------------------------------
          // CALCULATE THE CURRENTS (metal_update_term) AT (t-1/2) AND (t-1-1/2)
          //---------------------------------------------------------
          // i. these were derived outside of Taflove, but see eqs (9.62) and (9.52)
          // ii. this terms have the appearance that the currents are placed on
          // same time grid as E(t-1/2).  However, in the derivation of these the 
          // currents actually exist at (t)

          // set the metal_update_term to the Drude part
          metal_update_term = (1.0 + drude_kp[material])*drude_current_x[i][j][k];

          // now add in Lorentz terms
          for(l=0;l<2;++l)
          {
            metal_update_term += (1.0 + lorentz_alphap[material][l])*lorentz_currentp_x[l][i][j][k] + lorentz_zetap[material][l]*lorentz_currentp_x_nm1[l][i][j][k];
          }  

          // now we need to multiply by -1/2
          // i. see the derivation if unclear
          metal_update_term *= -0.5;

          //---------------------------------------------------------
          // UPDATE THE E FIELD TO (t+1/2)
          //---------------------------------------------------------
          field_ex[i][j][k] = (1.0/metal_c1[material])*((field_hz[i][j][k] - field_hz[i][j-1][k])/kedy[j] - (field_hy[i][j][k] - field_hy[i][j][k-1])/kedz[k] + psi_exy[i][j][k] - psi_exz[i][j][k] + metal_update_term + metal_c2[material]*e_n + metal_c3[material]*field_ex_nm1[i][j][k]);

          //---------------------------------------------------------
          // UPDATE THE CURRENTS TO (t+1/2)
          //---------------------------------------------------------
          // i. these were derived outside of Taflove, but see eqs (9.58) and (9.48)
          // ii. these are updated to the same time grid as E
   
          // Drude term
          drude_current_x[i][j][k] = drude_kp[material]*drude_current_x[i][j][k] + drude_betap[material]*(field_ex[i][j][k] + e_n);

          // Lorentz terms
          for(l=0;l<2;++l)
          {
            lorentz_currentp_x[l][i][j][k] = lorentz_alphap[material][l]*j_n[l] + lorentz_zetap[material][l]*lorentz_currentp_x_nm1[l][i][j][k] + lorentz_gammap[material][l]*(field_ex[i][j][k] - field_ex_nm1[i][j][k])/(2.0*time_dt);
          }  

          //---------------------------------------------------------
          // UPDATE THE STORED (t-1/2 next step) CURRENTS
          //---------------------------------------------------------
          field_ex_nm1[i][j][k] = e_n;

          for(l=0;l<2;++l)
          {
            lorentz_currentp_x_nm1[l][i][j][k] = j_n[l];
          } 

        }
        // ---------------------------- REAL MATERIAL ---------------------------- 
        else
        {
          //---------------------------------------------------------
          // UPDATE THE E FIELD TO (t+1/2)
          //---------------------------------------------------------
          field_ex[i][j][k] += time_dt*((field_hz[i][j][k] - field_hz[i][j-1][k])/kedy[j] - (field_hy[i][j][k] - field_hy[i][j][k-1])/kedz[k] + psi_exy[i][j][k] - psi_exz[i][j][k])/permittivity_ex[i][j][k];
        }


        //=========================================================
        // Ey
        //=========================================================
 
        // SET THE MATERIAL
        material = grid_material_ey[i][j][k];

        // ---------------------------- LOSSY MATERIAL ---------------------------- 
        if(material != 99)
        {
          //---------------------------------------------------------
          // STORE E & J AT (t-1/2)
          //---------------------------------------------------------
          e_n = field_ey[i][j][k];

          for(l=0;l<2;++l)
          {
            j_n[l] = lorentz_currentp_y[l][i][j][k];
          }   

          //---------------------------------------------------------
          // CALCULATE THE CURRENTS (metal_update_term) AT (t-1/2) AND (t-1-1/2)
          //---------------------------------------------------------
 
          // set the metal_update_term to the Drude part
          metal_update_term = (1.0 + drude_kp[material])*drude_current_y[i][j][k];

          // now add in Lorentz terms
          for(l=0;l<2;++l)
          {
            metal_update_term += (1.0 + lorentz_alphap[material][l])*lorentz_currentp_y[l][i][j][k] + lorentz_zetap[material][l]*lorentz_currentp_y_nm1[l][i][j][k];
          }  

          metal_update_term *= -0.5;

          //---------------------------------------------------------
          // UPDATE THE E FIELD TO (t+1/2)
          //---------------------------------------------------------
          field_ey[i][j][k] = (1.0/metal_c1[material])*((field_hx[i][j][k] - field_hx[i][j][k-1])/kedz[k] - (field_hz[i][j][k] - field_hz[i-1][j][k])/kedx[i] + psi_eyz[i][j][k] - psi_eyx[i][j][k] + metal_update_term + metal_c2[material]*e_n + metal_c3[material]*field_ey_nm1[i][j][k]);

          //---------------------------------------------------------
          // UPDATE THE CURRENTS TO (t+1/2)
          //---------------------------------------------------------

          // Drude term
          drude_current_y[i][j][k] = drude_kp[material]*drude_current_y[i][j][k] + drude_betap[material]*(field_ey[i][j][k] + e_n);

          // Lorentz terms
          for(l=0;l<2;++l)
          {
            lorentz_currentp_y[l][i][j][k] = lorentz_alphap[material][l]*j_n[l] + lorentz_zetap[material][l]*lorentz_currentp_y_nm1[l][i][j][k] + lorentz_gammap[material][l]*(field_ey[i][j][k] - field_ey_nm1[i][j][k])/(2.0*time_dt);
          }  

          //---------------------------------------------------------
          // UPDATE THE STORED (t-1/2 next step) CURRENTS
          //---------------------------------------------------------
          field_ey_nm1[i][j][k] = e_n;

          for(l=0;l<2;++l)
          {
            lorentz_currentp_y_nm1[l][i][j][k] = j_n[l];
          } 

        }
        // ---------------------------- REAL MATERIAL ---------------------------- 
        else
        {
          //---------------------------------------------------------
          // UPDATE THE E FIELD TO (t+1/2)
          //---------------------------------------------------------
          field_ey[i][j][k] += time_dt*((field_hx[i][j][k] - field_hx[i][j][k-1])/kedz[k] - (field_hz[i][j][k] - field_hz[i-1][j][k])/kedx[i] + psi_eyz[i][j][k] - psi_eyx[i][j][k])/permittivity_ey[i][j][k];
        }

        //=========================================================
        // Ez
        //=========================================================
 
        // SET THE MATERIAL
        material = grid_material_ez[i][j][k];

        // ---------------------------- LOSSY MATERIAL ---------------------------- 
        if(grid_material_ez[i][j][k] != 99)
        {
          //---------------------------------------------------------
          // STORE E & J AT (t-1/2)
          //---------------------------------------------------------
          e_n = field_ez[i][j][k];

          for(l=0;l<2;++l)
          {
            j_n[l] = lorentz_currentp_z[l][i][j][k];
          }   

          //---------------------------------------------------------
          // CALCULATE THE CURRENTS (metal_update_term) AT (t-1/2) AND (t-1-1/2)
          //---------------------------------------------------------
 
          // set the metal_update_term to the Drude part
          metal_update_term = (1.0 + drude_kp[material])*drude_current_z[i][j][k];

          // now add in Lorentz terms
          for(l=0;l<2;++l)
          {
            metal_update_term += (1.0 + lorentz_alphap[material][l])*lorentz_currentp_z[l][i][j][k] + lorentz_zetap[material][l]*lorentz_currentp_z_nm1[l][i][j][k];
          }  

          metal_update_term *= -0.5;

          //---------------------------------------------------------
          // UPDATE THE E FIELD TO (t+1/2)
          //---------------------------------------------------------
          field_ez[i][j][k] = (1.0/metal_c1[material])*((field_hy[i][j][k] - field_hy[i-1][j][k])/kedx[i] - (field_hx[i][j][k] - field_hx[i][j-1][k])/kedy[j] + psi_ezx[i][j][k] - psi_ezy[i][j][k] + metal_update_term + metal_c2[material]*e_n + metal_c3[material]*field_ez_nm1[i][j][k]);

          //---------------------------------------------------------
          // UPDATE THE CURRENTS TO (t+1/2)
          //---------------------------------------------------------

          // Drude term
          drude_current_z[i][j][k] = drude_kp[material]*drude_current_z[i][j][k] + drude_betap[material]*(field_ez[i][j][k] + e_n);

          // Lorentz terms
          for(l=0;l<2;++l)
          {
            lorentz_currentp_z[l][i][j][k] = lorentz_alphap[material][l]*j_n[l] + lorentz_zetap[material][l]*lorentz_currentp_z_nm1[l][i][j][k] + lorentz_gammap[material][l]*(field_ez[i][j][k] - field_ez_nm1[i][j][k])/(2.0*time_dt);
          }  

          //---------------------------------------------------------
          // UPDATE THE STORED (t-1/2 next step) CURRENTS
          //---------------------------------------------------------
          field_ez_nm1[i][j][k] = e_n;

          for(l=0;l<2;++l)
          {
            lorentz_currentp_z_nm1[l][i][j][k] = j_n[l];
          } 

        }
        // ---------------------------- REAL MATERIAL ---------------------------- 
        else
        {
          field_ez[i][j][k] += time_dt*((field_hy[i][j][k] - field_hy[i-1][j][k])/kedx[i] - (field_hx[i][j][k] - field_hx[i][j-1][k])/kedy[j] + psi_ezx[i][j][k] - psi_ezy[i][j][k])/permittivity_ez[i][j][k];
        }


      } // i++
    } // j++
  } // k++


  //=========================================================
  // UPDATE PSI TERMS
  //=========================================================
  // i. update the psi terms used for the H updates because the these need to be 
  // on the same time step as E

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        psi_hxy[i][j][k] = bh_y[j]*psi_hxy[i][j][k] + ch_y[j]*(field_ez[i][j+1][k] - field_ez[i][j][k])/grid_dy;
        psi_hxz[i][j][k] = bh_z[k]*psi_hxz[i][j][k] + ch_z[k]*(field_ey[i][j][k+1] - field_ey[i][j][k])/grid_dz;
        psi_hyx[i][j][k] = bh_x[i]*psi_hyx[i][j][k] + ch_x[i]*(field_ez[i+1][j][k] - field_ez[i][j][k])/grid_dx;
        psi_hyz[i][j][k] = bh_z[k]*psi_hyz[i][j][k] + ch_z[k]*(field_ex[i][j][k+1] - field_ex[i][j][k])/grid_dz;
        psi_hzx[i][j][k] = bh_x[i]*psi_hzx[i][j][k] + ch_x[i]*(field_ey[i+1][j][k] - field_ey[i][j][k])/grid_dx;
        psi_hzy[i][j][k] = bh_y[j]*psi_hzy[i][j][k] + ch_y[j]*(field_ex[i][j+1][k] - field_ex[i][j][k])/grid_dy;
      }
    }
  }


  //=========================================================
  // TF/SF UPDATE
  //=========================================================
  // i. for normal incidence we only need to update Ex & Ey:
  //     field_ex[i][j][k](T) += field_ex_coeff[i][j][k]*(0.0 - (field_hy[i][j][k](T) - field_hy[i][j][k-1](S))/grid_dz);
  //     field_ey[i][j][k](T) += field_ey_coeff[i][j][k]*((field_hx[i][j][k](T) - field_hx[i][j][k-1](S))/grid_dz);	
  //     field_ex[i][j][k](T) += field_ex_coeff[i][j][k]*Inc_hy(k-1)/grid_dz;
  //     field_ey[i][j][k](T) += field_ey_coeff[i][j][k]*(0.0 - Inc_hx(k-1))/grid_dz;	
  // !fix: should this be placed before the psi updates?

  if(mpi_itfsf_k1 == 1)
  {
    k = mpi_tfsf_k1 - mpi_my_grid_nz1 + 1;

    double hx_inc, hy_inc;
    hx_inc = hx_aux[aux_tfsf_k1-1];
    hy_inc = hy_aux[aux_tfsf_k1-1];

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        field_ex[i][j][k] += (time_dt/permittivity_ex[i][j][k])*hy_inc/kedz[k];
        field_ey[i][j][k] += (  (time_dt/permittivity_ey[i][j][k])*(0.0 - hx_inc)/kedz[k] );	
      }
    }

  } // end if


  //=========================================================
  // SEND/RECEIVE BOUNDARY PLANES
  //=========================================================
  // i. we need to send the lowest plane values to the processor directly below,
  // which will be received as the highest plane values
  // ii. what is needed are:
  //
  //      SURFACE 1 (k SURFACE): {ex(i,j, k+1) & ey(i,j, k+1)} surface vals  
  //      SURFACE 2 (j SURFACE): {ex(i,j+1, k) & ez(i,j+1, k)} surface vals  
  //      SURFACE 3 (i SURFACE): {ey(i+1,j, k) & ez(i+1,j, k)} surface vals  

  // counter for filling the surface values
  int counter;

  // COMMUNICATION PARAMETERS 
  MPI_Request mpi_send_request_iey, mpi_recv_request_iey, mpi_send_request_iez, mpi_recv_request_iez, mpi_send_request_jex, mpi_recv_request_jex, mpi_send_request_jez, mpi_recv_request_jez, mpi_send_request_kex, mpi_recv_request_kex, mpi_send_request_key, mpi_recv_request_key;
  MPI_Status mpi_status;
  int mpi_tag;

  // SET THE PLANE TO SEND TO BE THE LOWEST
  int plane=1;

  //---------------------------------------------------------  
  // FILL SURFACES
  //---------------------------------------------------------  

/*
  // i-SURFACE
  if( (ixperiodic == 0) && (mpi_my_grid_coords[0] == 0) )
  {
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;
      
        isurface_y_send[counter] = 0.0;
        isurface_z_send[counter] = 0.0;
        counter++;
      }
    }

  }
  else
  {
*/
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;
      
        isurface_y_send[counter] = field_ey[plane][j][k];
        isurface_z_send[counter] = field_ez[plane][j][k];
        counter++;
      }
    }
//  }

/*
  // j-SURFACE
  if( (iyperiodic == 0) && (mpi_my_grid_coords[1] == 0) )
  {
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        jsurface_x_send[counter] = 0.0;
        jsurface_z_send[counter] = 0.0;
        counter++;
      }
    }
  }
  else
  {
*/
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        jsurface_x_send[counter] = field_ex[i][plane][k];
        jsurface_z_send[counter] = field_ez[i][plane][k];
        counter++;
      }
    }
//  }


/*
  // k-SURFACE
  if( (izperiodic == 0) && (mpi_my_grid_coords[2] == 0) )
  {
    counter = 0;
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        ksurface_x_send[counter] = 0.0;
        ksurface_y_send[counter] = 0.0;
        counter++;
      }
    }
  }
  else
  {
*/
    counter = 0;
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        ksurface_x_send[counter] = field_ex[i][j][plane];
        ksurface_y_send[counter] = field_ey[i][j][plane];
        counter++;
      }
    }
//  }


  //---------------------------------------------------------  
  // SEND & RECEIVE SURFACES
  //---------------------------------------------------------  

  if(isendtype==1)
  {
    // i-SURFACE
    // Ey
    mpi_tag = 11;
    MPI_Isend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm, &mpi_send_request_iey);
    MPI_Irecv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_iey);
    // Ez
    mpi_tag = 12;
    MPI_Isend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm, &mpi_send_request_iez);
    MPI_Irecv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_iez);

    // j-SURFACE
    // Ex
    mpi_tag = 13;
    MPI_Isend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm, &mpi_send_request_jex);
    MPI_Irecv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_jex);
    // Ez
    mpi_tag = 14;
    MPI_Isend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm, &mpi_send_request_jez);
    MPI_Irecv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_jez);

    // k-SURFACE
    // Ex
    mpi_tag = 15;
    MPI_Isend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm, &mpi_send_request_kex);
    MPI_Irecv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_kex);
    // Ey
    mpi_tag = 16;
    MPI_Isend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm, &mpi_send_request_key);
    MPI_Irecv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_key);

    // MAKE SURE ALL SEND/RECEIVES ARE DONE
    // !fix: not sure optimal place to put these or if I should use MPI_Test instead
    // !fix: instead of doing this it may be better to use arrays and issue MPI_Waitall command
    // !fix: also not sure if it would be better to first wait for sends
    // then do receieves, or do send/recv/send/recv in the order they were issued
    // i SURFACE
    MPI_Wait(&mpi_send_request_iey, &mpi_status);
    MPI_Wait(&mpi_send_request_iez, &mpi_status);
    MPI_Wait(&mpi_recv_request_iey, &mpi_status);
    MPI_Wait(&mpi_recv_request_iez, &mpi_status);
    // j SURFACE
    MPI_Wait(&mpi_send_request_jex, &mpi_status);
    MPI_Wait(&mpi_send_request_jez, &mpi_status);
    MPI_Wait(&mpi_recv_request_jex, &mpi_status);
    MPI_Wait(&mpi_recv_request_jez, &mpi_status);
    // k SURFACE
    MPI_Wait(&mpi_send_request_kex, &mpi_status);
    MPI_Wait(&mpi_send_request_key, &mpi_status);
    MPI_Wait(&mpi_recv_request_kex, &mpi_status);
    MPI_Wait(&mpi_recv_request_key, &mpi_status);
  }
  else if(isendtype==2)
  {
    // i-SURFACE
    // Ey
    mpi_tag = 11;
    MPI_Sendrecv(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag,mpi_grid_comm, &mpi_status);
    // Ez
    mpi_tag = 12;
    MPI_Sendrecv(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag,mpi_grid_comm, &mpi_status);

    // j-SURFACE
    // Ex
    mpi_tag = 13;
    MPI_Sendrecv(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag,mpi_grid_comm, &mpi_status);
    // Ez
    mpi_tag = 14;
    MPI_Sendrecv(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag,mpi_grid_comm, &mpi_status);

    // k-SURFACE
    // Ex
    mpi_tag = 15;
    MPI_Sendrecv(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag,mpi_grid_comm, &mpi_status);
    // Ey
    mpi_tag = 16;
    MPI_Sendrecv(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev

    // i-SURFACE
    if(mpi_my_grid_coords[0] % 2 == 0)
    {
      // Ey
      mpi_tag = 11;
      MPI_Ssend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm);
      MPI_Recv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
      // Ez
      mpi_tag = 12;
      MPI_Ssend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm);
      MPI_Recv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      // Ey
      mpi_tag = 11;
      MPI_Recv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm);
      // Ez
      mpi_tag = 12;
      MPI_Recv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_efield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_efield_isend, mpi_tag, mpi_grid_comm);
    }

    // j-SURFACE
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[1] % 2 == 0)
    {
      // Ex
      mpi_tag = 13;
      MPI_Ssend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm);
      MPI_Recv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
      // Ez
      mpi_tag = 14;
      MPI_Ssend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm);
      MPI_Recv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      // Ex
      mpi_tag = 13;
      MPI_Recv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm);
      // Ez
      mpi_tag = 14;
      MPI_Recv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_efield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_efield_jsend, mpi_tag, mpi_grid_comm);
    }

    // k-SURFACE
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[2] % 2 == 0)
    {
      // Ex
      mpi_tag = 15;
      MPI_Ssend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm);
      MPI_Recv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
      // Ey
      mpi_tag = 16;
      MPI_Ssend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm);
      MPI_Recv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      // Ex
      mpi_tag = 15;
      MPI_Recv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm);
      // Ey
      mpi_tag = 16;
      MPI_Recv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_efield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_efield_ksend, mpi_tag, mpi_grid_comm);
    }

  }


  //=========================================================
  // ADD IN RECEIVED VALUES
  //=========================================================
  // i. the values need to be added in at 1 above the highest range

  // i-SURFACE
  counter = 0;
  j = mpi_my_grid_ny2 - mpi_my_grid_ny1 + 2;
  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
    {
      i = ii - mpi_my_grid_nx1 + 1;

      field_ex[i][j][k] = jsurface_x_recv[counter];
      field_ez[i][j][k] = jsurface_z_recv[counter];
      counter++;
    }
  }

  // j-SURFACE
  counter = 0;
  i = mpi_my_grid_nx2 - mpi_my_grid_nx1 + 2;
  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      field_ey[i][j][k] = isurface_y_recv[counter];
      field_ez[i][j][k] = isurface_z_recv[counter];
      counter++;
    }
  }

  // k-SURFACE
  counter = 0;
  k = mpi_my_grid_nz2 - mpi_my_grid_nz1 + 2;
  for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
  {
    j = jj - mpi_my_grid_ny1 + 1;

    for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
    {
      i = ii - mpi_my_grid_nx1 + 1;

      field_ex[i][j][k] = ksurface_x_recv[counter];
      field_ey[i][j][k] = ksurface_y_recv[counter];
      counter++;
    }
  }


  //=========================================================
  // CLEANUP & RETURN
  //=========================================================
  delete [] j_n;


  return;
}




//========================================================================
//========================================================================
//
//	NAME:	void field_e_aux_update()
//	DESC:	Updates the e field auxillary grid
//
//	NOTES: 	i. this is called by all processes
//
//		ii. recall that the (time) staggering is:
//
//			E(i){t-1/2 --> t+1/2}
//			H(i+1/2){t --> t+1}
//
//		iii. recall that the (space) staggering is:
//
//			Ex (i+1/2, j, k)
//			Ey (i, j+1/2, k)
//
//			Hx (i, j+1/2, k+1/2)
//			Hy (i+1/2, j, k+1/2)
//
//		iii. after some testing the following points have become clear:
//
//			ii.a. the following two sets of components should be used for the auxillary grid
//			(two sets are used so we can use circular polarization):
//
//				dHx/dt = (1/mu)*dEy/dz
//				dEy/dt = (1/eps)*dHx/dz
//
//				dHy/dt = -(1/mu)*dEx/dz
//				dEx/dt = -(1/eps)*dHy/dz				
//
//			iii.b. a hard source should be used with + sign should be used for both grids
//			(this will not really matter since components will propagate in both directions):
//
//				Ey(tfsf-2) = Ey source
//				Ex(tfsf-2) = Ex source
//
//				iii.b.1. note that for circular polarization there should be a pi/2 amplitude
//				difference in Ey and Ex
//
//				iii.b.1. using this hard source field components will end up
//				propagating in both +x and -x with the right amplitude
//					iii.b.1.a. using both Ey and Hz hard sources I still get waves
//					propagating in both directions, with the right amplitude
//					iii.b.1.b. using Ey and Hz soft sources I get propagation in the
//					right direction, but with the wrong amplitude
//
//			iii.c. !!! due to the sign of the equations in ii.a. the add in signs for the TF/SF
//			should be (+Ex, +Hy) and (+Ey, -Hx) 
//
//			iv. ! this subroutine is set up ONLY for normal incidence (i.e. no z-component fields)
//
//========================================================================
//========================================================================
void field_e_aux_update()
{

  // local indices
  int k, m;

  double ex_hard_src, ey_hard_src;


  //---------------------------------------------------------
  // UPDATE ENTIRE 1D FDTD GRID
  //---------------------------------------------------------
  // i. the first and last points are never updated, but the fields
  // will never propagate that far anyway 
  // ii. we also update with grid_dz
  for(k = 2; k <= (aux_npts-1); ++k)
  {
    ex_aux[k] -= e_aux_coeff*(hy_aux[k] - hy_aux[k-1])/grid_dz;
    ey_aux[k] += e_aux_coeff*(hx_aux[k] - hx_aux[k-1])/grid_dz;
  }


  //---------------------------------------------------------
  // ADD IN THE HARD SOURCE
  //---------------------------------------------------------
  // i. we are adding this in 2 points behind what we are taking
  // to be the tfsf k1 surface in line with what Taflove describes

  // get the effective time
  // !!! I am pretty sure this should be t+dt/2
  double t_effective = time_current + time_dt/2.0;

  // now get the gaussian param
  double gauss_param = (t_effective - source_gauss_center)*(t_effective - source_gauss_center)/(source_gauss_width*source_gauss_width);

 
  if(src_ieypol==1)
  {
    ex_hard_src = 0.0;
    ey_hard_src = source_intensity*exp(0.0-gauss_param)*sin(source_wavenumber*t_effective);
  }
  else if(src_iexpol==1)
  {
    ey_hard_src = 0.0;
    ex_hard_src = source_intensity*exp(0.0-gauss_param)*sin(source_wavenumber*t_effective);
  } 
  else if(src_i45pol==1)
  {
    ex_hard_src=ey_hard_src=source_intensity*exp(0.0-gauss_param)*sin(source_wavenumber*t_effective)/sqrt(2.0);
  } 
  // i. if we have circular polarization we have both components that are off-centered by pi/2
  else if(src_icircpol==1)
  {
    ex_hard_src = source_intensity*exp(0.0-gauss_param)*sin(source_wavenumber*t_effective + PI/2.0);
    ey_hard_src = source_intensity*exp(0.0-gauss_param)*sin(source_wavenumber*t_effective);
  } 
  else
  {
    cout << "no source used in computation!" << endl;
    ex_hard_src=ey_hard_src=0.0;
  }
 

  // now add in the hard source to the e grid
  ex_aux[aux_tfsf_k1 - 2] = ex_hard_src;
  ey_aux[aux_tfsf_k1 - 2] = ey_hard_src;


  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR THE TRANSMISSION SPECTRUM
  //---------------------------------------------------------
  // i. this is at the aux_spect_transm_kpt
  // ii. with normal incidence we will never have ez, but with general
  // incidence we would 
  if( (itransm==1) || (iscat_calc == 1) )
  {
    for(m = 1; m <= spect_npts; ++m)
    {
      aux_ftex_transm[m] += dft_multiplier_ht[m][istep]*ex_aux[aux_spect_transm_kpt];
      aux_ftey_transm[m] += dft_multiplier_ht[m][istep]*ey_aux[aux_spect_transm_kpt];
    }
  }

  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR SCATTERING CALC
  //---------------------------------------------------------

  if(iscat_calc == 1)
  {

    int k1 = aux_tfsf_k1 + (scat_k1pt - mpi_tfsf_k1) - 5;
    int k2 = aux_tfsf_k1 + (scat_k2pt - mpi_tfsf_k1) + 5;

    if(k2>aux_npts)
    {
      cout << "error code 652976." << endl;
    }

    for(m = 1; m <= spect_npts; ++m)
    {

      // only ft enough points to cover the domain of our simulation
      for(k = k1; k <= k2; ++k)
      {
        scat_aux_ex[m][k] += dft_multiplier_ht[m][istep]*ex_aux[k];
        scat_aux_ey[m][k] += dft_multiplier_ht[m][istep]*ey_aux[k];
      }

    }
  }


  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR THE FT SPECTRA
  //---------------------------------------------------------

  for(m=0;m<nft_planes;++m)
  {
    aux_ftex[m] += exp(0.0 - COMPLEXJ*ft_wavenumbers[m]*(time_current + time_dt/2.0))*ex_aux[ft_aux_coord[m]];
    aux_ftey[m] += exp(0.0 - COMPLEXJ*ft_wavenumbers[m]*(time_current + time_dt/2.0))*ey_aux[ft_aux_coord[m]];
  }


  //---------------------------------------------------------
  // FT FOR NSOM
  //---------------------------------------------------------
  // i. !!! this is not exactly correct because I am doing FT at the spectra transmission point
  //aux_ftex_nsom += 0.0;
//  aux_ftey_nsom += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current + time_dt/2.0))*ey_aux[aux_spect_transm_kpt];
  //aux_ftez_nsom += 0.0;


  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void field_h_update()
//	DESC:	updates the h-field
//
//	NOTES: 	i.) this is called by all processes, but each process only
//		operates on its own grid section
//
//
//========================================================================
//========================================================================
void field_h_update()
{
 
  // "local" indices for loops
  int i, ii, j, jj, k, kk;

  //=========================================================================
  // UPDATE ENTIRE GRID
  //=========================================================================
  // i.) !!! note in this update we do not include the points on the boundaries that
  // rely on other processors
  // ii.) note we update things differently depending on if we are a PML processor
  // or an interior processor


  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        field_hx[i][j][k] += time_dt*((field_ey[i][j][k+1] - field_ey[i][j][k])/khdz[k] - (field_ez[i][j+1][k] - field_ez[i][j][k])/khdy[j] + psi_hxz[i][j][k] - psi_hxy[i][j][k] )/permeability_hx[i][j][k];
        field_hy[i][j][k] += time_dt*((field_ez[i+1][j][k] - field_ez[i][j][k])/khdx[i] - (field_ex[i][j][k+1] - field_ex[i][j][k])/khdz[k] + psi_hyx[i][j][k] - psi_hyz[i][j][k] )/permeability_hy[i][j][k];
        field_hz[i][j][k] += time_dt*((field_ex[i][j+1][k] - field_ex[i][j][k])/khdy[j] - (field_ey[i+1][j][k] - field_ey[i][j][k])/khdx[i] + psi_hzy[i][j][k] - psi_hzx[i][j][k] )/permeability_hz[i][j][k];

      } // ii++
    } // jj++
  } // kk++



  //=========================================================================
  // UPDATE PSI TERMS
  //=========================================================================
  // i. now that we have updated the H fields we need to update the psi terms to use in the 
  // E updates because we need these psi terms on the same footing as the H

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        psi_exy[i][j][k] = be_y[j]*psi_exy[i][j][k] + ce_y[j]*(field_hz[i][j][k] - field_hz[i][j-1][k])/grid_dy;
        psi_exz[i][j][k] = be_z[k]*psi_exz[i][j][k] + ce_z[k]*(field_hy[i][j][k] - field_hy[i][j][k-1])/grid_dz;

        psi_eyx[i][j][k] = be_x[i]*psi_eyx[i][j][k] + ce_x[i]*(field_hz[i][j][k] - field_hz[i-1][j][k])/grid_dx;
        psi_eyz[i][j][k] = be_z[k]*psi_eyz[i][j][k] + ce_z[k]*(field_hx[i][j][k] - field_hx[i][j][k-1])/grid_dz;

        psi_ezx[i][j][k] = be_x[i]*psi_ezx[i][j][k] + ce_x[i]*(field_hy[i][j][k] - field_hy[i-1][j][k])/grid_dx;
        psi_ezy[i][j][k] = be_y[j]*psi_ezy[i][j][k] + ce_y[j]*(field_hx[i][j][k] - field_hx[i][j-1][k])/grid_dy;
      }
    }
  }


  //=========================================================
  // UPDATE THE BOTTOM TF/SF SURFACE
  //=========================================================

  // IF THIS PROCESSOR IS RESPONSIBLE FOR BOTTOM TF/SF 
  // ??? i.) for this region the [k+1] value needs to be updates (it is a TF and 
  // needs to be updated to a SF)
  // ??? ii.) TF = SF + IF

//			Ex (i+1/2, j, k)
//			Ey (i, j+1/2, k)
//
//			Hx (i, j+1/2, k+1/2)
//			Hy (i+1/2, j, k+1/2)

// ** at tfsf_z1 all are TF
// ** at tfsf_z1+any all are TF
// **!!! at tfsf_z1-1 H are SF and E[k+1] are TF

//field_hx[i][j][k] += hcoefficient*((field_ey[i][j][k+1] - field_ey[i][j][k])/grid_dz);
//field_hy[i][j][k] += hcoefficient*(0.0 - (field_ex[i][j][k+1] - field_ex[i][j][k])/grid_dz);

//field_hx[i][j][k](S) += hcoefficient*(((field_ey[i][j][k+1](T) - Inc_ey) - field_ey[i][j][k](S))/grid_dz);
//field_hy[i][j][k](S) += hcoefficient*(0.0 - ((field_ex[i][j][k+1](T)-Inc_ex) - field_ex[i][j][k](S))/grid_dz);

//field_hx[i][j][k](S) += hcoefficient*(0.0 - Inc_ey(k = tfsf_z1-1+1))/grid_dz;
//field_hy[i][j][k](S) += hcoefficient*Inc_ex(k = tfsf_z1-1+1)/grid_dz;


  if(mpi_itfsf_k1 == 1)
  {

    // upml !!! this needs to be changed to be more generic (with free space around tf/sf it will be ok though)
    double hcoefficient = time_dt/MU0;

    // i.) we need to update at at tfsf_z1-1
    // ii.) the commented out expression is more explanatory
    //k = mpi_tfsf_k1 - 1 - mpi_my_grid_nz1 + 1;
    k = mpi_tfsf_k1 - mpi_my_grid_nz1;

    // GET THE INCIDENT FIELD
    // i.) we get the incident field from the auxillary grid
    // iii.) remember we need the field [at the k1 surface], but 
    // update [at the k1-1 surface]
    double ex_inc, ey_inc;

    ex_inc = ex_aux[aux_tfsf_k1];    
    ey_inc = ey_aux[aux_tfsf_k1];


    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;
	
        // i.) note if we only want one type of polarization we don't need
        // to update both of these 
        field_hx[i][j][k] += ( hcoefficient*(0.0 - ey_inc)/khdz[k] );
        // !!! the polarization below is taken by looking at the MW equations
        // and the hx above
        field_hy[i][j][k] += hcoefficient*ex_inc/khdz[k];
      }
    }

  } // end if




  //=========================================================
  // SEND/RECEIVE SURFACES
  //=========================================================
  // to do all of this we note the following points

  // i.) there are 6 surfaces to for each process(or)
  // ii.) we only need to send out 3 of the 6 surfaces 
  // iii.) we thus only need to receive 3 of the 6 surfaces
  // iv.) for these h updates we need the lowest value from the received e-fields


  // BELOW IS WHAT WE NEED FOR THE h FIELDS
  //
  // WHAT WE NOW NEED TO DO IS SEND OUR HIGHEST VALUES TO THE PROCESS(OR) THAT 
  // IS RANKED RIGHT ABOVE US AND ASSIGN IT TO THE LOWEST VALUE
  //
  //-----------------------------------------------------------
  // SURFACE 1
  // we need the {hx(i,j, k-1) & hy(i,j, k-1)} surface vals  

  // SURFACE 2
  // we need the {hx(i,j-1, k) & hz(i,j-1, k)} surface vals  

  // SURFACE 3
  // we need the {hy(i-1,j, k) & hz(i-1,j, k)} surface vals  
  //-----------------------------------------------------------

//field_ex[i][j][k] += field_ex_coeff[i][j][k]*((field_hz[i][j][k] - field_hz[i][j-1][k])/grid_dy - (field_hy[i][j][k] - field_hy[i][j][k-1])/grid_dz);
//field_ey[i][j][k] += field_ey_coeff[i][j][k]*((field_hx[i][j][k] - field_hx[i][j][k-1])/grid_dz - (field_hz[i][j][k] - field_hz[i-1][j][k])/grid_dx);	
//field_ez[i][j][k] += field_ez_coeff[i][j][k]*((field_hy[i][j][k] - field_hy[i-1][j][k])/grid_dx - (field_hx[i][j][k] - field_hx[i][j-1][k])/grid_dy);


  // MISC 

  // counter for filling the surface values
  int counter;

  //---------------------------------------------------------  
  // MPI COMMUNICATION PARAMETERS
  //--------------------------------------------------------- 

  // SENDING PARAMETERS 

  // MISC
  
  MPI_Status mpi_status;

  int mpi_tag;

  // !!! instead of doing this it may be better to use
  // arrays and issue MPI_Waitall command

  MPI_Request mpi_send_request_ihy, mpi_recv_request_ihy;
  MPI_Request mpi_send_request_ihz, mpi_recv_request_ihz;

  MPI_Request mpi_send_request_jhx, mpi_recv_request_jhx;
  MPI_Request mpi_send_request_jhz, mpi_recv_request_jhz;

  MPI_Request mpi_send_request_khx, mpi_recv_request_khx;
  MPI_Request mpi_send_request_khy, mpi_recv_request_khy;

  //=========================================================
  // SEND AND RECEIVE SURFACES
  //=========================================================


  //---------------------------------------------------------  
  // i SURFACE
  //---------------------------------------------------------  

  // GET INFO READY FOR SEND
/*
  // reset the counter
  if( (ixperiodic == 0) && (mpi_my_grid_coords[0] == (mpi_nxprocs-1)) )
  {
    i = mpi_my_grid_nx2 - mpi_my_grid_nx1 + 1;

    // fill up the array to send
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        isurface_y_send[counter] = 0.0;
        isurface_z_send[counter] = 0.0;

        counter++;
      }
    }

  }
  else
  {
*/
    i = mpi_my_grid_nx2 - mpi_my_grid_nx1 + 1;

    // fill up the array to send
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        isurface_y_send[counter] = field_hy[i][j][k];
        isurface_z_send[counter] = field_hz[i][j][k];

        counter++;
      }
    }
//  }

  // PERFROM NON-BLOCKING SEND AND RECEIVES
  // sendtest

  // hy
  mpi_tag = 21;

  if(isendtype==1)
  {

    MPI_Isend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm, &mpi_send_request_ihy);
    MPI_Irecv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_ihy);
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[0] % 2 == 0)
    {
      MPI_Ssend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm);
      MPI_Recv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(isurface_y_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(isurface_y_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm);
    }
  }


  // hz
  mpi_tag = 22;

  if(isendtype==1)
  {
    MPI_Isend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm, &mpi_send_request_ihz);
    MPI_Irecv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_ihz);;
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[0] % 2 == 0)
    {
      MPI_Ssend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm);
      MPI_Recv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(isurface_z_recv, isurface_nvals, MPI_DOUBLE, mpi_hfield_irecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(isurface_z_send, isurface_nvals, MPI_DOUBLE, mpi_hfield_isend, mpi_tag, mpi_grid_comm);
    }
  }


  //---------------------------------------------------------
  // DO THE j SURFACE
  //---------------------------------------------------------

/*
  // GET INFO READY FOR SEND
  if( (iyperiodic == 0) && (mpi_my_grid_coords[1] == (mpi_nyprocs-1)) )
  {
    j = mpi_my_grid_ny2 - mpi_my_grid_ny1 + 1;

    // fill up the array to send
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        jsurface_x_send[counter] = 0.0;
        jsurface_z_send[counter] = 0.0;

        counter++;
      }
    }
  }
  else
  {
*/
    j = mpi_my_grid_ny2 - mpi_my_grid_ny1 + 1;

    // fill up the array to send
    counter = 0;
    for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
    {
      k = kk - mpi_my_grid_nz1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        jsurface_x_send[counter] = field_hx[i][j][k];
        jsurface_z_send[counter] = field_hz[i][j][k];

        counter++;
      }
    }
//  }

  // PERFROM NON-BLOCKING SEND AND RECEIVES
  // sendtest

  // hx
  mpi_tag = 23;

  if(isendtype==1)
  {
    MPI_Isend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm, &mpi_send_request_jhx);
    MPI_Irecv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_jhx);
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[1] % 2 == 0)
    {
      MPI_Ssend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm);
      MPI_Recv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(jsurface_x_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(jsurface_x_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm);
    }
  }

  // hz
  mpi_tag = 24;

  if(isendtype==1)
  {
    MPI_Isend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm, &mpi_send_request_jhz);
    MPI_Irecv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_jhz);
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[1] % 2 == 0)
    {
      MPI_Ssend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm);
      MPI_Recv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(jsurface_z_recv, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jrecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(jsurface_z_send, jsurface_nvals, MPI_DOUBLE, mpi_hfield_jsend, mpi_tag, mpi_grid_comm);
    }
  }


  //---------------------------------------------------------
  // DO THE k SURFACE
  //---------------------------------------------------------

/*
  if( (izperiodic == 0) && (mpi_my_grid_coords[2] == (mpi_nzprocs-1)) )
  {
    k = mpi_my_grid_nz2 - mpi_my_grid_nz1 + 1;

    // fill up the array to send
    counter = 0;
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        ksurface_x_send[counter] = 0.0;
        ksurface_y_send[counter] = 0.0;

        counter++;
      }
    }
  }
  else
  {
*/
    k = mpi_my_grid_nz2 - mpi_my_grid_nz1 + 1;

    // fill up the array to send
    counter = 0;
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        ksurface_x_send[counter] = field_hx[i][j][k];
        ksurface_y_send[counter] = field_hy[i][j][k];

        counter++;
      }
    }
//  }

  // PERFROM NON-BLOCKING SEND AND RECEIVES
  // sendtest

  // hx
  mpi_tag = 25;

  if(isendtype==1)
  {
    MPI_Isend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm, &mpi_send_request_khx);
    MPI_Irecv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_khx);
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[2] % 2 == 0)
    {
      MPI_Ssend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm);
      MPI_Recv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(ksurface_x_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(ksurface_x_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm);
    }
  }


  // hy
  mpi_tag = 26;

  if(isendtype==1)
  {
    MPI_Isend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm, &mpi_send_request_khy);
    MPI_Irecv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_recv_request_khy);
  }
  else if(isendtype==2)
  {
    MPI_Sendrecv(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag,mpi_grid_comm, &mpi_status);
  }
  else if(isendtype==3)
  {
    // with syncronous send we must be careful of sending order because the syncronous send will wait for beginning to receiev
    if(mpi_my_grid_coords[2] % 2 == 0)
    {
      MPI_Ssend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm);
      MPI_Recv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
    }
    else
    {
      MPI_Recv(ksurface_y_recv, ksurface_nvals, MPI_DOUBLE, mpi_hfield_krecv, mpi_tag, mpi_grid_comm, &mpi_status);
      MPI_Ssend(ksurface_y_send, ksurface_nvals, MPI_DOUBLE, mpi_hfield_ksend, mpi_tag, mpi_grid_comm);
    }
  }



  //---------------------------------------------------------
  // WAIT FOR ALL COMMUNICATIONS TO FINISH
  //---------------------------------------------------------
  // sendtest

  // now as a final thing we should finalize the send and receive calls with MPI_Wait
  // !!! not sure optimal place to put these or if I should use MPI_Test instead

  // !!! instead of doing this it may be better to use
  // arrays and issue MPI_Waitall command


  // FIRST MAKE SURE ALL SENDS ARE DONE

  if(isendtype==1)
  {
    // i SURFACE
    MPI_Wait(&mpi_send_request_ihy, &mpi_status);
    MPI_Wait(&mpi_send_request_ihz, &mpi_status);
    // j SURFACE
    MPI_Wait(&mpi_send_request_jhx, &mpi_status);
    MPI_Wait(&mpi_send_request_jhz, &mpi_status);
    // k SURFACE
    MPI_Wait(&mpi_send_request_khx, &mpi_status);
    MPI_Wait(&mpi_send_request_khy, &mpi_status);


    // NOW MAKE SURE ALL RECEIEVES ARE DONE

    // i SURFACE
    MPI_Wait(&mpi_recv_request_ihy, &mpi_status);
    MPI_Wait(&mpi_recv_request_ihz, &mpi_status);
    // j SURFACE
    MPI_Wait(&mpi_recv_request_jhx, &mpi_status);
    MPI_Wait(&mpi_recv_request_jhz, &mpi_status);
    // k SURFACE
    MPI_Wait(&mpi_recv_request_khx, &mpi_status);
    MPI_Wait(&mpi_recv_request_khy, &mpi_status);
  }

  //=========================================================
  // ADD IN RECEIVED VALUES
  //=========================================================
  // i.) we are going to be adding in the lowest position to our
  // local arrays, which will be 0


  // NOW ADD IN THE RECEIVED VALUES
  counter = 0;

  i = 0;

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      field_hy[i][j][k] = isurface_y_recv[counter];
      field_hz[i][j][k] = isurface_z_recv[counter];

      // increase counter
      counter += 1;
    }
  }


  // NOW ADD IN THE RECEIVED VALUES
  counter = 0;

  j = 0;

  for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
  {
    k = kk - mpi_my_grid_nz1 + 1;

    for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
    {
      i = ii - mpi_my_grid_nx1 + 1;

      field_hx[i][j][k] = jsurface_x_recv[counter];
      field_hz[i][j][k] = jsurface_z_recv[counter];

      // increase counter
      counter += 1;
    }
  }


  // NOW ADD IN THE RECEIVED VALUES
  // !!! although if we are the lowest processor we may need to figure
  // out something else to do because of no periodic BC in z
  counter = 0;

  k = 0;

  for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
  {
    j = jj - mpi_my_grid_ny1 + 1;

    for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
    {
      i = ii - mpi_my_grid_nx1 + 1;

      field_hx[i][j][k] = ksurface_x_recv[counter];
      field_hy[i][j][k] = ksurface_y_recv[counter];

      // increase counter
      counter += 1;
    }
  }


  //=========================================================
  // CLEANUP
  //=========================================================



  return;

}



//========================================================================
//========================================================================
//
//	NAME:	void field_h_aux_update()
//	DESC:	Update the h_field auxillary grid
//
//	NOTES: 	i. this is called by all processes
//
//		ii. recall that the (time) staggering is:
//
//			E(i){t-1/2 --> t+1/2}
//			H(i+1/2){t --> t+1}
//
//		iii. recall that the (space) staggering is:
//
//			Ex (i+1/2, j, k)
//			Ey (i, j+1/2, k)
//
//			Hx (i, j+1/2, k+1/2)
//			Hy (i+1/2, j, k+1/2)
//
//		ii. see the notes in e_aux_update() for descriptions on why what is done is done 
//
//========================================================================
//========================================================================
void field_h_aux_update()
{

  // indices
  int k, m;

  //---------------------------------------------------------
  // UPDATE ENTIRE 1D FDTD GRID
  //---------------------------------------------------------
  // i. the first and last points are never updated, but the fields
  // will never propagate that far anyway 
  // ii. we also update with grid_dz
  for(k = 2; k <= (aux_npts-1); ++k)
  {
    hx_aux[k] += h_aux_coeff*(ey_aux[k+1] - ey_aux[k])/grid_dz;
    hy_aux[k] -= h_aux_coeff*(ex_aux[k+1] - ex_aux[k])/grid_dz;
  }


  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR THE TRANSMISSION SPECTRUM
  //---------------------------------------------------------
  // i. this is at the aux_spect_transm_kpt
  // ii. with normal incidence we will never have hz, but with general
  // incidence we would
  if( (itransm==1) || (iscat_calc == 1) )
  { 
    for(m = 1; m <= spect_npts; ++m)
    {
      aux_fthx_transm[m] += dft_multiplier_t[m][istep]*hx_aux[aux_spect_transm_kpt];
      aux_fthy_transm[m] += dft_multiplier_t[m][istep]*hy_aux[aux_spect_transm_kpt];
    }
  }

  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR SCATTERING CALC
  //---------------------------------------------------------

  int k1 = aux_tfsf_k1 + (scat_k1pt - mpi_tfsf_k1) - 5;
  int k2 = aux_tfsf_k1 + (scat_k2pt - mpi_tfsf_k1) + 5;

  if(iscat_calc == 1)
  {

    for(m = 1; m <= spect_npts; ++m)
    {

      // only ft enough points to cover the domain of our simulation
      for(k = k1; k <= k2; ++k)
      {
        scat_aux_hx[m][k] += dft_multiplier_t[m][istep]*hx_aux[k];
        scat_aux_hy[m][k] += dft_multiplier_t[m][istep]*hy_aux[k];
      }

    }
  }

  //---------------------------------------------------------
  // NOW GET THE INCIDENT FIELD FOR THE FT SPECTRA
  //---------------------------------------------------------

  for(m=0;m<nft_planes;++m)
  {
    aux_fthx[m] += exp(0.0 - COMPLEXJ*ft_wavenumbers[m]*(time_current + time_dt))*hx_aux[ft_aux_coord[m]];
    aux_fthy[m] += exp(0.0 - COMPLEXJ*ft_wavenumbers[m]*(time_current + time_dt))*hy_aux[ft_aux_coord[m]];
  }

  //---------------------------------------------------------
  // FT FOR NSOM
  //---------------------------------------------------------
  // i. !!! this is not exactly correct because I am doing FT at the spectra transmission point
  aux_fthx_nsom += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current + time_dt))*hx_aux[aux_spect_transm_kpt];
  //aux_fthy_nsom += 0.0;


  return;
}



//========================================================================
//========================================================================
//
//	NAME:	void ft_fields()
//	DESC:	FT the fields on a surface
//
//	NOTES: 	i.) this is called by processes that contain the FT surface
//
//        // A X B = [a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1]
        // POYNTING VECTOR = 0.5[E X H*]
//
//
//	S_z = a1*b2 - a2*b1 = Ex*conj(Hy) - Ey*conj(Hx)
//		ii. because the dielectric function evaluates to a NEGATIVE imaginary part
//		the Fourier transforms have to be such that they are done with exp(-i*omega*t),
//		and the time component has exp(i*omega*t)
//
//
//========================================================================
//========================================================================
void ft_fields_spect()
{

  // GET "LOCAL" INDICES
  int i, ii, j, jj, k, m;

  //=========================================================
  // TRANSMISSION
  //=========================================================
  if(mpi_ft_spect_transm == 1)
  {

    k = spect_transm_kpt - mpi_my_grid_nz1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          i = ii - mpi_my_grid_nx1 + 1;

          // !!! we also note that since we are projecting onto the z direction only we don't need
          // to store all of the components but can just do the sum below
          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          transm_ftex[m][i][j] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          transm_ftey[m][i][j] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];

          transm_fthx[m][i][j] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          transm_fthy[m][i][j] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m

  } // end if(mpi_ft_spect_transm == 1)


  //=========================================================
  // REFLECTION
  //=========================================================
  if(mpi_ft_spect_refl == 1)
  {

    k = spect_refl_kpt - mpi_my_grid_nz1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          i = ii - mpi_my_grid_nx1 + 1;

          // !!! we also note that since we are projecting onto the z direction only we don't need
          // to store all of the components but can just do the sum below
          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          refl_ftex[m][i][j] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          refl_ftey[m][i][j] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];

          refl_fthx[m][i][j] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          refl_fthy[m][i][j] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m

  } // end if(mpi_ft_spect_refl == 1)


  return;
}






//========================================================================
//========================================================================
//
//	NAME:	void output(int ioutput_type, int iplane, int iunit)
//	DESC:	outputs a file
//
//	NOTES: 	i.) this is called by a process(or) responsible for 
//		outputting
//
//	output_field[iplane] is the plane we wish to output (xy, xz, or yz)
//
// 		ii. types:
// 
//			1 == e_total
//			2 == ex
//			3 == ey
//			4 == ez
//			5 == h_total
//			6 == hx
//			7 == hy
//			8 == hz
//
//			9 == |ft_etotal|^2/|ft_einc|^2 - !!! only xy,xz
//			10 == |ft_ex|^2/|ft_einc|^2 - !!! only xy,xz
//			11 == |ft_ey|^2/|ft_einc|^2 - !!! only xy,xz
//			12 == |ft_ez|^2/|ft_einc|^2 - !!! only xy,xz
//
//			21 == |ft_btotal|^2/|ft_binc|^2 - !!! only yz,xz
//
//			13 == ft_px - !!! not done
//
//	output_plane[iplane] is the type of the field we want to output
//	(ex, ey, ez or ee)
//
//
//		iii. types of planes:
//			1 == xy
//			2 == xz
//			3 == yz
//
//		iv. I have removed field averaging because in order to average end points 
//		of processor responsibilities we would need info from neighboring processors
//
//		v. !!! the setting up file names, etc I think is very sloppy
//		
//		!!! we can probably just use the ft e_inc from the aux grid (at the appropriate wavenumber)
//		for the FT planes
//		vi. H fields are used for B outputs by multiplication by the permeability
//
//
//========================================================================
//========================================================================
void output(int ioutput_type, int iplane, int iunit)
{

  // output fields
  if(ioutput_type == 1)
  {
    int i, ii, j, jj, k, kk;

    double e_total;
    complex<double> ft_e_total, ft_e_inc;


    //=========================================================
    // SET-UP FILE (FILENAME, ETC)
    //=========================================================
    // i. the filename will be set differently depending on the type of output format:
    //      1 == gnuplot: e_total.plane#.iunit.grid_coord1.grid_coord2
    //      2 == OpenDX: e_total.plane#.iunit.grid_coordx.grid_coordy.grid_coordz

    char filename_1[] = "./output/e_field.";
    char filename_2[5];
    char filename_3[5];
    char filename_4[5];
    char filename_5[5];
    char filename_6[5];

    char filename[31];	

    // SET UP iplane AND iunit PART OF FILENAME
    sprintf(filename_2,"%i",iplane); 
    sprintf(filename_3,"%i",iunit);
   
    // SET UP COORDINATE PORTIONS
    // Gnuplot
    if(output_format[iplane] == 1)
    {
      // xy
      if(output_plane[iplane] == 1) 
      {
        sprintf(filename_4,"%i",mpi_my_grid_coords[0]);  
        sprintf(filename_5,"%i",mpi_my_grid_coords[1]);
      }
      // xz
      else if(output_plane[iplane] == 2)
      {
        sprintf(filename_4,"%i",mpi_my_grid_coords[0]);  
        sprintf(filename_5,"%i",mpi_my_grid_coords[2]);   
      }
      // yz
      else if(output_plane[iplane] == 3)
      {
        sprintf(filename_4,"%i",mpi_my_grid_coords[1]);  
        sprintf(filename_5,"%i",mpi_my_grid_coords[2]);   
      }
    }
    // OpenDX
    else if(output_format[iplane] == 2)
    {
      sprintf(filename_4,"%i",mpi_my_grid_coords[0]);  
      sprintf(filename_5,"%i",mpi_my_grid_coords[1]); 
      sprintf(filename_6,"%i",mpi_my_grid_coords[2]); 
    }

    // NOW COPY ALL FILENAMES TOGETHER
    strcpy(filename, filename_1);
    strcat(filename, filename_2);
    strcat(filename, ".");
    strcat(filename, filename_3);
    strcat(filename, ".");
    strcat(filename, filename_4);
    strcat(filename, ".");
    strcat(filename, filename_5);

    // OpenDX specific
    if(output_format[iplane] == 2)
    {
      strcat(filename, ".");
      strcat(filename, filename_6);
    }

    // OPEN FILE
    ofstream file(filename,ios::app);


    //=========================================================
    // WRITE FILE
    //=========================================================


    // Gnuplot
    if(output_format[iplane] == 1)
    {

      // -------------------------------------- xy -------------------------------------
      if(output_plane[iplane] == 1) 
      {
        // get the local k coordinate
        k = output_plane_local_coord[iplane];

        // i.) note we start at the lower value + 1 because we are taking
        // an average and need the value right before
        // i.) these loops are backwards from all other loops because
        // of how we want the first to be const and the second to vary
        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
          {
            // get the "local" k and i grid points
            j = jj - mpi_my_grid_ny1 + 1;
            i = ii - mpi_my_grid_nx1 + 1;


            // GET THE FIELD WE WANT FOR OUTPUT

            // total e-field
            if(output_field[iplane] == 1)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k] + field_ey[i][j][k]*field_ey[i][j][k]+ field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // ex
            else if(output_field[iplane] == 2)
            {
            e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k])/source_intensity;
            }
            // ey
            else if(output_field[iplane] == 3)
            {
              e_total = sqrt(field_ey[i][j][k]*field_ey[i][j][k])/source_intensity;
            }
            // ez
            else if(output_field[iplane] == 4)
            {
              e_total = sqrt(field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // total h-field
            else if(output_field[iplane] == 5)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k] + field_hy[i][j][k]*field_hy[i][j][k]+ field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // hx
            else if(output_field[iplane] == 6)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k])/source_intensity;
            }
            // hy
            else if(output_field[iplane] == 7)
            {
              e_total = sqrt(field_hy[i][j][k]*field_hy[i][j][k])/source_intensity;
            }
            // hz
            else if(output_field[iplane] == 8)
            {
              e_total = sqrt(field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // |E_total|^2/|E_inc|^2
            else if(output_field[iplane] == 9)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_ex[output_to_ft[iplane]][i][j]) + field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_ey[output_to_ft[iplane]][i][j]) + field_ft_ez[output_to_ft[iplane]][i][j]*conj(field_ft_ez[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ex|^2/|E_inc|^2
            else if(output_field[iplane] == 10)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_ex[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ey|^2/|E_inc|^2
            else if(output_field[iplane] == 11)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_ey[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ez|^2/|E_inc|^2
            else if(output_field[iplane] == 12)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ez[output_to_ft[iplane]][i][j]*conj(field_ft_ez[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |H_total|^2/|H_inc|^2
            else if(output_field[iplane] == 13)
            {
              ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]) + field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) + field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hx|^2/|H_inc|^2
            else if(output_field[iplane] == 14)
            {
              ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hy|^2/|H_inc|^2
            else if(output_field[iplane] == 15)
            {
              ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hzl|^2/|H_inc|^2
            else if(output_field[iplane] == 16)
            {
              ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |B_total|^2/|B_inc|^2
            // i. remember we multiply H by the permeability to get B
            else if(output_field[iplane] == 21)
            {
              ft_e_inc = (MU0*aux_mur)*(MU0*aux_mur)*(aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = permeability_hx[i][j][k]*permeability_hx[i][j][k]*field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]) + permeability_hy[i][j][k]*permeability_hy[i][j][k]*field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) + permeability_hz[i][j][k]*permeability_hz[i][j][k]*field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]) ;
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
// j12
            // Pz
            else if(output_field[iplane] == 51)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) - aux_ftey[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) ;
              ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) - field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }

            // write the file
            file << (static_cast<double>(ii))*grid_dx << "     " << (static_cast<double>(jj))*grid_dy << "     " << e_total << endl; 
          }

          // leave an extra line for gnuplot
          file << endl;
        }

      }
      // -------------------------------------- xz ---------------------------------------
      else if(output_plane[iplane] == 2)
      {
        // GET THE LOCAL y-POINT
        j = output_plane_local_coord[iplane];
 
        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
          {
            // get the "local" k and i grid points
            k = kk - mpi_my_grid_nz1 + 1;
            i = ii - mpi_my_grid_nx1 + 1;

            // GET THE FIELD WE WANT FOR OUTPUT

            // total e-field
            if(output_field[iplane] == 1)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k] + field_ey[i][j][k]*field_ey[i][j][k]+ field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // ex
            else if(output_field[iplane] == 2)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k])/source_intensity;
            }
            // ey
            else if(output_field[iplane] == 3)
            {
              e_total = sqrt(field_ey[i][j][k]*field_ey[i][j][k])/source_intensity;
            }
            // ez
            else if(output_field[iplane] == 4)
            {
              e_total = sqrt(field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // total h-field
            else if(output_field[iplane] == 5)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k] + field_hy[i][j][k]*field_hy[i][j][k]+ field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // hx
            else if(output_field[iplane] == 6)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k])/source_intensity;
            }
            // hy
            else if(output_field[iplane] == 7)
            {
              e_total = sqrt(field_hy[i][j][k]*field_hy[i][j][k])/source_intensity;
            }
            // hz
            else if(output_field[iplane] == 8)
            {
              e_total = sqrt(field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // |ft_e|^2
            else if(output_field[iplane] == 9)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ex[output_to_ft[iplane]][i][k]*conj(field_ft_ex[output_to_ft[iplane]][i][k]) + field_ft_ey[output_to_ft[iplane]][i][k]*conj(field_ft_ey[output_to_ft[iplane]][i][k]) + field_ft_ez[output_to_ft[iplane]][i][k]*conj(field_ft_ez[output_to_ft[iplane]][i][k]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |ft_ex|^2
            else if(output_field[iplane] == 10)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ex[output_to_ft[iplane]][i][k]*conj(field_ft_ex[output_to_ft[iplane]][i][k]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |ft_ey|^2
            else if(output_field[iplane] == 11)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ey[output_to_ft[iplane]][i][k]*conj(field_ft_ey[output_to_ft[iplane]][i][k]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |ft_ez|^2
            else if(output_field[iplane] == 12)
            {
              ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              ft_e_total = field_ft_ez[output_to_ft[iplane]][i][k]*conj(field_ft_ez[output_to_ft[iplane]][i][k]);
              e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |B_total|^2/|B_inc|^2
            // i. remember we multiply H by the permeability to get B
            else if(output_field[iplane] == 21)
            {
              ft_e_inc = (MU0*aux_mur)*(MU0*aux_mur)*(aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = permeability_hx[i][j][k]*permeability_hx[i][j][k]*field_ft_hx[output_to_ft[iplane]][i][k]*conj(field_ft_hx[output_to_ft[iplane]][i][k]) + permeability_hy[i][j][k]*permeability_hy[i][j][k]*field_ft_hy[output_to_ft[iplane]][i][k]*conj(field_ft_hy[output_to_ft[iplane]][i][k]) + permeability_hz[i][j][k]*permeability_hz[i][j][k]*field_ft_hz[output_to_ft[iplane]][i][k]*conj(field_ft_hz[output_to_ft[iplane]][i][k]) ;
              e_total = real(ft_e_total)/real(ft_e_inc);
            }

            // write the file
            file << (static_cast<double>(ii))*grid_dx << "     " << (static_cast<double>(kk))*grid_dz << "     " << e_total << endl; 
          }

          // leave an extra line for gnuplot
          file << endl;
        }

      }
      // ------------------------------- yz ----------------------------
      else if(output_plane[iplane] == 3)
      {
        // GET THE LOCAL x-POINT
        i = output_plane_local_coord[iplane];

        for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
        {
          for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
          {
            // get the "local" k and j grid points
            k = kk - mpi_my_grid_nz1 + 1;
            j = jj - mpi_my_grid_ny1 + 1;


            // GET THE FIELD WE WANT FOR OUTPUT

            // total e-field
            if(output_field[iplane] == 1)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k] + field_ey[i][j][k]*field_ey[i][j][k]+ field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // ex
            else if(output_field[iplane] == 2)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k])/source_intensity;
            }
            // ey
            else if(output_field[iplane] == 3)
            {
              e_total = sqrt(field_ey[i][j][k]*field_ey[i][j][k])/source_intensity;
            }
            // ez
            else if(output_field[iplane] == 4)
            {
              e_total = sqrt(field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // total h-field
            else if(output_field[iplane] == 5)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k] + field_hy[i][j][k]*field_hy[i][j][k]+ field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // hx
            else if(output_field[iplane] == 6)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k])/source_intensity;
            }
            // hy
            else if(output_field[iplane] == 7)
            {
              e_total = sqrt(field_hy[i][j][k]*field_hy[i][j][k])/source_intensity;
            }
            // hz
            else if(output_field[iplane] == 8)
            {
              e_total = sqrt(field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // ft total e-field
            else if(output_field[iplane] == 9)
            {

            }
            else if(output_field[iplane] == 10)
            {

            }
            else if(output_field[iplane] == 15)
            {

            }
            // |B_total|^2/|B_inc|^2
            // i. remember we multiply H by the permeability to get B
            else if(output_field[iplane] == 21)
            {
              ft_e_inc = (MU0*aux_mur)*(MU0*aux_mur)*(aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              ft_e_total = permeability_hx[i][j][k]*permeability_hx[i][j][k]*field_ft_hx[output_to_ft[iplane]][j][k]*conj(field_ft_hx[output_to_ft[iplane]][j][k]) + permeability_hy[i][j][k]*permeability_hy[i][j][k]*field_ft_hy[output_to_ft[iplane]][j][k]*conj(field_ft_hy[output_to_ft[iplane]][j][k]) + permeability_hz[i][j][k]*permeability_hz[i][j][k]*field_ft_hz[output_to_ft[iplane]][j][k]*conj(field_ft_hz[output_to_ft[iplane]][j][k]) ;
              e_total = real(ft_e_total)/real(ft_e_inc);
            }

            // write the file
            file << (static_cast<double>(jj))*grid_dy << "     " << (static_cast<double>(kk))*grid_dz << "     " << e_total << endl; 
          }

          // leave an extra line for gnuplot
          file << endl;
        }
      }

    }
    // OpenDX
    else if(output_format[iplane] == 2)
    {


      // WRITE THE DIMENSIONS FOR OUTPUT
      file << mpi_my_grid_nx1 << "   " << mpi_my_grid_nx2 << endl;       
      file << mpi_my_grid_ny1 << "   " << mpi_my_grid_ny2 << endl;  
      file << mpi_my_grid_nz1 << "   " << mpi_my_grid_nz2 << endl;  

      // i.) note we start at the lower value + 1 because we are taking
      // an average and need the value right before
      // i.) these loops are backwards from all other loops because
      // of how we want the first to be const and the second to vary
      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
        {
          for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
          {
            // get the "local" k and i grid points
            j = jj - mpi_my_grid_ny1 + 1;
            i = ii - mpi_my_grid_nx1 + 1;
            k = kk - mpi_my_grid_nz1 + 1;

            // GET THE FIELD WE WANT FOR OUTPUT

            // total e-field
            if(output_field[iplane] == 1)
            {
              e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k] + field_ey[i][j][k]*field_ey[i][j][k]+ field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // ex
            else if(output_field[iplane] == 2)
            {
            e_total = sqrt(field_ex[i][j][k]*field_ex[i][j][k])/source_intensity;
            }
            // ey
            else if(output_field[iplane] == 3)
            {
              e_total = sqrt(field_ey[i][j][k]*field_ey[i][j][k])/source_intensity;
            }
            // ez
            else if(output_field[iplane] == 4)
            {
              e_total = sqrt(field_ez[i][j][k]*field_ez[i][j][k])/source_intensity;
            }
            // total h-field
            else if(output_field[iplane] == 5)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k] + field_hy[i][j][k]*field_hy[i][j][k]+ field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // hx
            else if(output_field[iplane] == 6)
            {
              e_total = sqrt(field_hx[i][j][k]*field_hx[i][j][k])/source_intensity;
            }
            // hy
            else if(output_field[iplane] == 7)
            {
              e_total = sqrt(field_hy[i][j][k]*field_hy[i][j][k])/source_intensity;
            }
            // hz
            else if(output_field[iplane] == 8)
            {
              e_total = sqrt(field_hz[i][j][k]*field_hz[i][j][k])/source_intensity;
            }
            // |E_total|^2/|E_inc|^2
            else if(output_field[iplane] == 9)
            {
              //ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              //ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_ex[output_to_ft[iplane]][i][j]) + field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_ey[output_to_ft[iplane]][i][j]) + field_ft_ez[output_to_ft[iplane]][i][j]*conj(field_ft_ez[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ex|^2/|E_inc|^2
            else if(output_field[iplane] == 10)
            {
              //ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              //ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_ex[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ey|^2/|E_inc|^2
            else if(output_field[iplane] == 11)
            {
              //ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              //ft_e_total = field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_ey[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Ez|^2/|E_inc|^2
            else if(output_field[iplane] == 12)
            {
              //ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_ftex[output_to_ft[iplane]]) + aux_ftey[output_to_ft[iplane]]*conj(aux_ftey[output_to_ft[iplane]]) + aux_ftez[output_to_ft[iplane]]*conj(aux_ftez[output_to_ft[iplane]]);
              //ft_e_total = field_ft_ez[output_to_ft[iplane]][i][j]*conj(field_ft_ez[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |H_total|^2/|H_inc|^2
            else if(output_field[iplane] == 13)
            {
              //ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              //ft_e_total = field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]) + field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) + field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hx|^2/|H_inc|^2
            else if(output_field[iplane] == 14)
            {
              //ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              //ft_e_total = field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hy|^2/|H_inc|^2
            else if(output_field[iplane] == 15)
            {
              //ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              //ft_e_total = field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |Hzl|^2/|H_inc|^2
            else if(output_field[iplane] == 16)
            {
              //ft_e_inc = (aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              //ft_e_total = field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
            // |B_total|^2/|B_inc|^2
            // i. remember we multiply H by the permeability to get B
            else if(output_field[iplane] == 21)
            {
              //ft_e_inc = (MU0*aux_mur)*(MU0*aux_mur)*(aux_fthx[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) + aux_fthy[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) + aux_fthz[output_to_ft[iplane]]*conj(aux_fthz[output_to_ft[iplane]]) );
              //ft_e_total = permeability_hx[i][j][k]*permeability_hx[i][j][k]*field_ft_hx[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]) + permeability_hy[i][j][k]*permeability_hy[i][j][k]*field_ft_hy[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) + permeability_hz[i][j][k]*permeability_hz[i][j][k]*field_ft_hz[output_to_ft[iplane]][i][j]*conj(field_ft_hz[output_to_ft[iplane]][i][j]) ;
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }
// j12
            // Pz
            else if(output_field[iplane] == 51)
            {
              //ft_e_inc = aux_ftex[output_to_ft[iplane]]*conj(aux_fthy[output_to_ft[iplane]]) - aux_ftey[output_to_ft[iplane]]*conj(aux_fthx[output_to_ft[iplane]]) ;
              //ft_e_total = field_ft_ex[output_to_ft[iplane]][i][j]*conj(field_ft_hy[output_to_ft[iplane]][i][j]) - field_ft_ey[output_to_ft[iplane]][i][j]*conj(field_ft_hx[output_to_ft[iplane]][i][j]);
              //e_total = real(ft_e_total)/real(ft_e_inc);
            }

            // WRITE THE FIELD
            file << e_total << endl; 
          } // ++kk
        } // ++jj
      } // ++ii


    }


    // CLOSE FILE
    file.close();
  }
  // field stitching data
  else if(ioutput_type == 2)
  {

    //=========================================================
    // STITCH FILE
    //=========================================================

    //---------------------------------------------------------
    // OPEN FILE
    //---------------------------------------------------------
    ofstream stitch_file("./output/field_stitch.dat", ios::app);


    //---------------------------------------------------------
    // WRITE FILE
    //---------------------------------------------------------
  
    // leave an initial spacer
    stitch_file << endl;

    // we first need to output a little bit of system info

    // output number of processors
    stitch_file << mpi_nxprocs << endl;    
    stitch_file << mpi_nyprocs << endl; 
    stitch_file << mpi_nzprocs << endl << endl; 

    // output number of processors
    stitch_file << grid_nx << endl;    
    stitch_file << grid_ny << endl; 
    stitch_file << grid_nz << endl << endl; 

    stitch_file << grid_dx << endl;    
    stitch_file << grid_dy << endl; 
    stitch_file << grid_dz << endl << endl; 

    // write number of output files per plane
    stitch_file << noutput << endl;

    // leave a spacer
    stitch_file << endl;

    // write number of planes
    stitch_file << output_nplanes << endl;

    // now loop over the planes
    for(int p = 1; p <= output_nplanes; ++p)
    {
      // OUTPUT THE FORMAT OF THE PLANE
      stitch_file << output_format[p] << endl;

      // OUTPUT THE TYPE OF PLANE
      stitch_file << output_plane[p] << endl;
    }


    // leave a spacer
    stitch_file << endl;


    //---------------------------------------------------------
    // CLOSE FILE
    //---------------------------------------------------------
    stitch_file.close();
 

  }
  // simulation data
  else if(ioutput_type == 3)
  {


    //=========================================================
    // SIMULATION FILE
    //=========================================================

    //---------------------------------------------------------
    // OPEN FILE
    //---------------------------------------------------------
    ofstream info_file("./output/sim_info.dat", ios::app);


    //---------------------------------------------------------
    // WRITE FILE
    //---------------------------------------------------------
  
    // leave an initial spacer
    info_file << endl;

    // SYSTEM INFO
    info_file << " ** SIMULATION INFO **" << endl;
    info_file << "*******************************" << endl << endl;  
    info_file << "total sim time: " << systime_total << "s" << endl << endl;  

    // output number of processors
    info_file << "NUMBER OF PROCESSORS" << endl;   
    info_file << "nx processors: " << mpi_nxprocs << endl;    
    info_file << "ny processors: " << mpi_nyprocs << endl;   
    info_file << "nz processors: " << mpi_nzprocs << endl << endl;  


    // SIMULATION INFO
    info_file << " ** SIMULATION INFO **" << endl; 
    info_file << endl; 

    // grid info
    info_file << "GRID" << endl; 
    info_file << "grid xsize: " << grid_xsize << "m" << endl;   
    info_file << "grid ysize: " << grid_ysize << "m" << endl;   
    info_file << "grid zsize: " << grid_zsize << "m" << endl << endl; 

    info_file << "grid dx: " << grid_dx << "m" << endl;   
    info_file << "grid dy: " << grid_dy << "m" << endl;    
    info_file << "grid dz: " << grid_dz << "m" << endl << endl; 

    // time info
    info_file << "TIME" << endl; 
    info_file << "total time: " << time_total << "s" << endl; 
    info_file << "dt: " << time_dt << "s" << endl << endl; 

    info_file << "nsteps: " << nsteps << endl << endl; 

    // source info
    info_file << "SOURCE" << endl; 
//    info_file << "incident angle: " << source_inc_angle << "deg" << endl; 
    info_file << "wavelength: " << source_wavelength << "m" << endl; 

    // PLANE INFO
    info_file << " ** PLANE INFO **" << endl << endl; 

    // write number of planes
    info_file << "number of output planes: " << output_nplanes << endl << endl;


    // now loop over the planes
    for(int p = 1; p <= output_nplanes; ++p)
    {
      // output the plane number
      info_file << "PLANE: " << p << endl; 
      info_file << endl;

      // output the type of plane we want to output
      info_file << "plane type is: ";

      if(output_plane[p] == 1)
      {
        info_file << "xy (z-plane)" << endl; 
        info_file << endl;
      }
      else if(output_plane[p] == 2)
      {
        info_file << "xz (y-plane)" << endl; 
        info_file << endl;
      }
      else if(output_plane[p] == 3)
      {
        info_file << "yz (x-plane)" << endl; 
        info_file << endl;
      }

      // output the mid position of this plane
      info_file << "plane is at: " << output_plane_midpos[p] << endl;
      info_file << endl;



      // output the type of plane
      info_file << "field type is: "; 

      if(output_field[p] == 1)
      {
        info_file << "e_total" << endl; 
      }
      else if(output_field[p] == 2)
      {
        info_file << "ex" << endl; 
      }
      else if(output_field[p] == 3)
      {
        info_file << "ey" << endl; 
      }
      else if(output_field[p] == 4)
      {
        info_file << "ez" << endl; 
      }
      else if(output_field[output_to_ft[p]] == 9)
      {
        info_file << "ft field e total" << endl; 
        info_file << "ft wavelength: " << ft_wavelengths[p] << endl;  
      } 
      else if(output_field[output_to_ft[p]] == 10)
      {
        info_file << "ft snapshot e total" << endl; 
        info_file << "ft wavelength: " << ft_wavelengths[p] << endl;  
      } 
      else if(output_field[output_to_ft[p]] == 11)
      {
        info_file << "ft snapshot Ez" << endl; 
        info_file << "ft wavelength: " << ft_wavelengths[p] << endl;  
      } 
  
      // leave a spacer
      info_file << endl;
    }


    //---------------------------------------------------------
    // CLOSE FILE
    //---------------------------------------------------------
    info_file.close();
 
  }



  return;

}

//========================================================================
//========================================================================
//
//	NAME:	void ft_fields()
//	DESC:	FT the fields on a surface
//
//	NOTES: 	i.) this is called by processes that contain the FT surface
//
//        // A X B = [a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1]
        // POYNTING VECTOR = 0.5[E X H*]
//
//
//	S_z = a1*b2 - a2*b1 = Ex*conj(Hy) - Ey*conj(Hx)
//
          // i.) we should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it
//
//========================================================================
//========================================================================
void ft_fields()
{

  // local indices
  int i, ii, j, jj, k, kk, p;

  complex<double> exp_ht, exp_t;


  for(p=0;p<nft_planes;++p)
  {

    // GET THE EXPONENTIAL PREFACTORS
    // i. !!! it would probably be a lot faster to declare this in an array like the
    // scattering exponential prefactors
    exp_ht = exp(0.0 - COMPLEXJ*ft_wavenumbers[p]*(time_current+time_dt/2.0));
    exp_t = exp(0.0 - COMPLEXJ*ft_wavenumbers[p]*(time_current+time_dt));



    // ----------------------- xy ------------------------------
    if(output_plane[ft_to_output[p]] == 1)
    {
      // GET THE LOCAL z-POINT
      k = output_plane_local_coord[ft_to_output[p]];

      for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          i = ii - mpi_my_grid_nx1 + 1;

          field_ft_ex[p][i][j] += exp_ht*field_ex[i][j][k];
          field_ft_ey[p][i][j] += exp_ht*field_ey[i][j][k];
          field_ft_ez[p][i][j] += exp_ht*field_ez[i][j][k];

          field_ft_hx[p][i][j] += exp_t*field_hx[i][j][k];
          field_ft_hy[p][i][j] += exp_t*field_hy[i][j][k];
          field_ft_hz[p][i][j] += exp_t*field_hz[i][j][k];
        }
      }


    }
    // -------------------------- xz -------------------------
    else if(output_plane[ft_to_output[p]] == 2)
    {
      // GET THE LOCAL y-POINT
      j = output_plane_local_coord[ft_to_output[p]];


      for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
      {
        k = kk - mpi_my_grid_nz1 + 1;

        for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
        {
          i = ii - mpi_my_grid_nx1 + 1;

          field_ft_ex[p][i][k] += exp_ht*field_ex[i][j][k];
          field_ft_ey[p][i][k] += exp_ht*field_ey[i][j][k];
          field_ft_ez[p][i][k] += exp_ht*field_ez[i][j][k];

          field_ft_hx[p][i][k] += exp_t*field_hx[i][j][k];
          field_ft_hy[p][i][k] += exp_t*field_hy[i][j][k];
          field_ft_hz[p][i][k] += exp_t*field_hz[i][j][k];
        }
      }

    }
    // -------------------------- yz -------------------------
    else if(output_plane[ft_to_output[p]] == 3)
    {
      // GET THE LOCAL x-POINT
      i = output_plane_local_coord[ft_to_output[p]];

      for(kk = mpi_my_grid_nz1; kk <= mpi_my_grid_nz2; ++kk)
      {
        k = kk - mpi_my_grid_nz1 + 1;

        for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
        {
          j = jj - mpi_my_grid_ny1 + 1;

          field_ft_ex[p][j][k] += exp_ht*field_ex[i][j][k];
          field_ft_ey[p][j][k] += exp_ht*field_ey[i][j][k];
          field_ft_ez[p][j][k] += exp_ht*field_ez[i][j][k];

          field_ft_hx[p][j][k] += exp_t*field_hx[i][j][k];
          field_ft_hy[p][j][k] += exp_t*field_hy[i][j][k];
          field_ft_hz[p][j][k] += exp_t*field_hz[i][j][k];
        }
      }

    }

  } // ++p



  return;
}




//========================================================================
//========================================================================
//
//	NAME:	void ft_field_nsom()
//	DESC:	FT the fields on a surface
//
//	NOTES: 	i.) this is called by processes that contain the FT surface
//
//		i. should multiply the integrals by time_dt, but since this occurs in the 
//		incident fields too we won't worry about it
//========================================================================
//========================================================================
void ft_field_nsom()
{
  // "local" indices
  int n;


  //=========================================================
  // TRANSMISSION
  //=========================================================

  // loop over all of the wavelength we need to FT
  for(n = 1; n <= transm_nsom_npts; ++n)
  {
    transm_nsom_ex[n] += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current+time_dt/2.0))*field_ex[transm_nsom_ptoi[n]][transm_nsom_ptoj[n]][transm_nsom_ptok];
    transm_nsom_ey[n] += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current+time_dt/2.0))*field_ey[transm_nsom_ptoi[n]][transm_nsom_ptoj[n]][transm_nsom_ptok];

    transm_nsom_hx[n] += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current+time_dt))*field_hx[transm_nsom_ptoi[n]][transm_nsom_ptoj[n]][transm_nsom_ptok];
    transm_nsom_hy[n] += exp(0.0 - COMPLEXJ*nsom_wavenumber*(time_current+time_dt))*field_hy[transm_nsom_ptoi[n]][transm_nsom_ptoj[n]][transm_nsom_ptok];
  }



  return;

}



//========================================================================
//========================================================================
//
//	NAME:	void ft_fields_scat()
//	DESC:	Fourier transform the fields around the cross section box
//
//	NOTES: 	i.
//
//
//========================================================================
//========================================================================
void ft_fields_scat()
{

  // local indices
  int m, i, ii, j, jj, k, kk;

  //=========================================================
  // i1
  //=========================================================
  if(scat_i1ft == 1)
  {
  
    i = scat_i1pt - mpi_my_grid_nx1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {


      for(jj = scat_i1y1; jj <= scat_i1y2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(kk = scat_i1z1; kk <= scat_i1z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_i1_ftey[m][j][k] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];
          scat_i1_ftez[m][j][k] += dft_multiplier_ht[m][istep]*field_ez[i][j][k];

          scat_i1_fthy[m][j][k] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
          scat_i1_fthz[m][j][k] += dft_multiplier_t[m][istep]*field_hz[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m

  }

  
  //=========================================================
  // i2
  //=========================================================
  if(scat_i2ft == 1)
  {


    i = scat_i2pt - mpi_my_grid_nx1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(jj = scat_i2y1; jj <= scat_i2y2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(kk = scat_i2z1; kk <= scat_i2z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_i2_ftey[m][j][k] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];
          scat_i2_ftez[m][j][k] += dft_multiplier_ht[m][istep]*field_ez[i][j][k];

          scat_i2_fthy[m][j][k] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
          scat_i2_fthz[m][j][k] += dft_multiplier_t[m][istep]*field_hz[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m


  }


  //=========================================================
  // j1
  //=========================================================
  if(scat_j1ft == 1)
  {


    j = scat_j1pt - mpi_my_grid_ny1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(ii = scat_j1x1; ii <= scat_j1x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(kk = scat_j1z1; kk <= scat_j1z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_j1_ftex[m][i][k] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          scat_j1_ftez[m][i][k] += dft_multiplier_ht[m][istep]*field_ez[i][j][k];

          scat_j1_fthx[m][i][k] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          scat_j1_fthz[m][i][k] += dft_multiplier_t[m][istep]*field_hz[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m



  }


  //=========================================================
  // j2
  //=========================================================
  if(scat_j2ft == 1)
  {


    j = scat_j2pt - mpi_my_grid_ny1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(ii = scat_j2x1; ii <= scat_j2x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(kk = scat_j2z1; kk <= scat_j2z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_j2_ftex[m][i][k] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          scat_j2_ftez[m][i][k] += dft_multiplier_ht[m][istep]*field_ez[i][j][k];

          scat_j2_fthx[m][i][k] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          scat_j2_fthz[m][i][k] += dft_multiplier_t[m][istep]*field_hz[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m



  }


  //=========================================================
  // k1
  //=========================================================
  if(scat_k1ft == 1)
  {

    k = scat_k1pt - mpi_my_grid_nz1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(ii = scat_k1x1; ii <= scat_k1x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(jj = scat_k1y1; jj <= scat_k1y2; ++jj)
        {
          j = jj - mpi_my_grid_ny1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_k1_ftex[m][i][j] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          scat_k1_ftey[m][i][j] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];

          scat_k1_fthx[m][i][j] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          scat_k1_fthy[m][i][j] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m


  }

  //=========================================================
  // k2
  //=========================================================
  if(scat_k2ft == 1)
  {

    k = scat_k2pt - mpi_my_grid_nz1 + 1;


    // loop over all of the wavelength we need to FT
    for(m = 1; m <= spect_npts; ++m)
    {

      for(ii = scat_k2x1; ii <= scat_k2x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(jj = scat_k2y1; jj <= scat_k2y2; ++jj)
        {
          j = jj - mpi_my_grid_ny1 + 1;

          // i.) should multiply by time_dt, but since this occurs in the incident fields too we won't
          // worry about it

          scat_k2_ftex[m][i][j] += dft_multiplier_ht[m][istep]*field_ex[i][j][k];
          scat_k2_ftey[m][i][j] += dft_multiplier_ht[m][istep]*field_ey[i][j][k];

          scat_k2_fthx[m][i][j] += dft_multiplier_t[m][istep]*field_hx[i][j][k];
          scat_k2_fthy[m][i][j] += dft_multiplier_t[m][istep]*field_hy[i][j][k];
        } // ++ii
      } // ++jj

    } // ++m

  }


  return;
}




//========================================================================
//========================================================================
//
//	NAME:	void calc_cs()
//	DESC:	Calculate the cross sections
//
//	NOTES: 	i.
//
//		ii. the calculations here assume normal incidence so not ez or
//		hz components
//
//
//========================================================================
//========================================================================
void calc_cs()
{


  // local indices
  int i, ii, j, jj, k, kk, m, p;
  int kaux;

  // poynting vector & h scattered component
  complex<double> poynting_integral;
  complex<double> etemp, etemp2, htemp, htemp2;

  double *pintegral;
  pintegral = new double [spect_npts+1];


  //=========================================================
  // CALCULATE SURFACE INTEGRAL OF S FOR EACH WAVELENGTH
  //=========================================================

  // for each m ...
  for(m = 1; m <= spect_npts; ++m)
  {

    // zero the poynting integral
    poynting_integral = COMPLEXZERO;


    //---------------------------------------------------------
    // i1 (normal S = -x)
    //---------------------------------------------------------
    // i. s(-x) = -(ey*hz-ez*hy)

    if(scat_i1ft == 1)
    {
  
      i = scat_i1pt - mpi_my_grid_nx1 + 1;

      for(jj = scat_i1y1; jj <= scat_i1y2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(kk = scat_i1z1; kk <= scat_i1z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          kaux = aux_tfsf_k1 + (kk - mpi_tfsf_k1);
          etemp = scat_i1_ftey[m][j][k]-scat_aux_ey[m][kaux];
          htemp = scat_i1_fthy[m][j][k]-scat_aux_hy[m][kaux];
          poynting_integral += grid_dy*grid_dz*(  0.0-(etemp*conj(scat_i1_fthz[m][j][k]) - scat_i1_ftez[m][j][k]*conj(htemp))  );

        } // ++kk
      } // ++jj

    }


  
    //---------------------------------------------------------
    // i2 (normal S = x)
    //---------------------------------------------------------
    // i. s(x) = ey*hz-ez*hy

    if(scat_i2ft == 1)
    {

      i = scat_i2pt - mpi_my_grid_nx1 + 1;

      for(jj = scat_i2y1; jj <= scat_i2y2; ++jj)
      {
        j = jj - mpi_my_grid_ny1 + 1;

        for(kk = scat_i2z1; kk <= scat_i2z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          kaux = aux_tfsf_k1 + (kk - mpi_tfsf_k1);
          etemp = scat_i2_ftey[m][j][k]-scat_aux_ey[m][kaux];
          htemp = scat_i2_fthy[m][j][k]-scat_aux_hy[m][kaux];
          poynting_integral += grid_dy*grid_dz*(etemp*conj(scat_i2_fthz[m][j][k]) - scat_i2_ftez[m][j][k]*conj(htemp));
        } // ++ii
      } // ++jj

    }



    //---------------------------------------------------------
    // j1 (normal S = -y)
    //---------------------------------------------------------
    // i. s(-y) = -(ez*hx-ex*hz)

    if(scat_j1ft == 1)
    {

      j = scat_j1pt - mpi_my_grid_ny1 + 1;

      for(ii = scat_j1x1; ii <= scat_j1x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(kk = scat_j1z1; kk <= scat_j1z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          kaux = aux_tfsf_k1 + (kk - mpi_tfsf_k1);
          etemp = scat_j1_ftex[m][i][k]-scat_aux_ex[m][kaux];
          htemp = scat_j1_fthx[m][i][k]-scat_aux_hx[m][kaux];

          poynting_integral += grid_dx*grid_dz*( 0.0-(scat_j1_ftez[m][i][k]*conj(htemp) - etemp*conj(scat_j1_fthz[m][i][k])) );
        } // ++ii
      } // ++jj


    }



    //---------------------------------------------------------
    // j2 (normal S = y)
    //---------------------------------------------------------
    // i. s(y) = (ez*hx-ex*hz)

    if(scat_j2ft == 1)
    {

      j = scat_j2pt - mpi_my_grid_ny1 + 1;

      for(ii = scat_j2x1; ii <= scat_j2x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(kk = scat_j2z1; kk <= scat_j2z2; ++kk)
        {
          k = kk - mpi_my_grid_nz1 + 1;

          kaux = aux_tfsf_k1 + (kk - mpi_tfsf_k1);
          etemp = scat_j2_ftex[m][i][k]-scat_aux_ex[m][kaux];
          htemp = scat_j2_fthx[m][i][k]-scat_aux_hx[m][kaux];
          poynting_integral += grid_dx*grid_dz*(scat_j2_ftez[m][i][k]*conj(htemp) - etemp*conj(scat_j2_fthz[m][i][k]));
        } // ++ii
      } // ++jj

    }



    //---------------------------------------------------------
    // k1 (normal S = -z)
    //---------------------------------------------------------
    // i. s(-z) = -(ex*hy - ey*hx)

    if(scat_k1ft == 1)
    {

      k = scat_k1pt - mpi_my_grid_nz1 + 1;

      kaux = aux_tfsf_k1 + (scat_k1pt - mpi_tfsf_k1);

      for(ii = scat_k1x1; ii <= scat_k1x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(jj = scat_k1y1; jj <= scat_k1y2; ++jj)
        {
          j = jj - mpi_my_grid_ny1 + 1;

          etemp = scat_k1_ftex[m][i][j]-scat_aux_ex[m][kaux];
          etemp2 = scat_k1_ftey[m][i][j]-scat_aux_ey[m][kaux];
          htemp = scat_k1_fthx[m][i][j]-scat_aux_hx[m][kaux];
          htemp2 = scat_k1_fthy[m][i][j]-scat_aux_hy[m][kaux];
          poynting_integral += grid_dx*grid_dy*(  0.0-(etemp*conj(htemp2) - etemp2*conj(htemp)) );
        } // ++ii
      } // ++jj

    }

    //---------------------------------------------------------
    // k2 (normal S = z)
    //---------------------------------------------------------
    // i. s(z) = (ex*hy - ey*hx)
    if(scat_k2ft == 1)
    {

      k = scat_k2pt - mpi_my_grid_nz1 + 1;

      kaux = aux_tfsf_k1 + (scat_k2pt - mpi_tfsf_k1);

      for(ii = scat_k2x1; ii <= scat_k2x2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        for(jj = scat_k2y1; jj <= scat_k2y2; ++jj)
        {
          j = jj - mpi_my_grid_ny1 + 1;

          etemp = scat_k2_ftex[m][i][j]-scat_aux_ex[m][kaux];
          etemp2 = scat_k2_ftey[m][i][j]-scat_aux_ey[m][kaux];
          htemp = scat_k2_fthx[m][i][j]-scat_aux_hx[m][kaux];
          htemp2 = scat_k2_fthy[m][i][j]-scat_aux_hy[m][kaux];
          poynting_integral += grid_dx*grid_dy*(etemp*conj(htemp2) - etemp2*conj(htemp));
        } // ++ii
      } // ++jj

    }

    // take 1/2 of the real part of the integral
    pintegral[m] = 0.5*real(poynting_integral);

    // also normalize to nm
    pintegral[m] /= ((1.0e-9)*(1.0e-9));

  } // ++m


  //=========================================================
  // SEND ALL INFO TO MASTER PROCESS
  //=========================================================
  // i.) note we are using blocking send/recv here
  // ii.) we are also accessing the MPI_COMM_WORLD communicator IDs not the grid one

  int mpi_source, mpi_destination;
  MPI_Status mpi_status;
  int mpi_tag;

  // SEND VALUES TO MASTER PROCESSOR
  if(mpi_my_id != mpi_master_id)
  {

    // set the destination to the master process
    mpi_destination = mpi_master_id;

    // SEND THE POYNTING VECTOR
    mpi_tag = 41;
    MPI_Send(pintegral, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
  }



  //=========================================================
  // RECEIVE INFO AT MASTER PROCESS AND CALCULATE PROPERTIES
  //=========================================================

  if(mpi_my_id == mpi_master_id)
  {


    //---------------------------------------------------------
    // RECEIVE POYNTING VECTORS
    //---------------------------------------------------------
    double *pintegral_recv;
    pintegral_recv = new double [spect_npts+1];


    // NOW RECEIVE AND SUM VALUES
    for(p = 1; p <= (mpi_nprocessors - 1); ++p) 
    {
      // SET THE SOURCE TO THE CURRENT PROCESS(OR)
      mpi_source = p;
 
      // RECEIVE POYNTING VECTOR
      mpi_tag = 41;
      MPI_Recv(pintegral_recv, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);


      // ADD IN THE RECEIVED VALUES TO OUR SUMMATION
      for(m = 1; m <= spect_npts; ++m)
      {
        pintegral[m] += pintegral_recv[m];
      } // ++m


    }  // ++p


    // CLEANUP
    delete [] pintegral_recv;


    //---------------------------------------------------------
    // NORMALIZE
    //---------------------------------------------------------
    // i. the incident flux is purely in the +z direction
    // ii. s(z) = (ex*hy - ey*hx)

    // incident poynting vector
    complex<double> pzinc;

    // loop over the transmission points
    for(m = 1; m <= spect_npts; ++m)
    {
      // get the incident poynting vector
      pzinc = 0.5*(  aux_ftex_transm[m]*conj(aux_fthy_transm[m]) - aux_ftey_transm[m]*conj(aux_fthx_transm[m])  );

      // normalize the poynting integral we have calculated
      pintegral[m] /= real(pzinc);
    } // ++m


    //---------------------------------------------------------
    // WRITE FILES
    //---------------------------------------------------------

    // OPEN FILES
    ofstream scatt_file("./output/scatt.dat",ios::app);

    // WRITE FILES
    double ft_wave;
    

    for(m = 1; m <= spect_npts; ++m)
    {
      // get the wavelength at this wavevector
      ft_wave = 2.0*PI*CSPEED/spect_wavenumbers[m];

      // write the normalized values
      scatt_file << ft_wave << "     " << pintegral[m] << endl;
    }


    // CLOSE FILES
    scatt_file.close();


  } // end if master id



  //=========================================================
  // CLEANUP
  //=========================================================
  delete [] pintegral;


  return;

}



//========================================================================
//========================================================================
//
//	NAME:	void calc_transm()
//	DESC:	calculates the transmittance spectrum
//
//	NOTES: 	i. this is called by all processes
//
//		ii. the cross product is defined as:
//
//			a x b = (ay*bz - az*by)i + 
//				(az*bx - ax*bz)j + 
//				(ax*by - ay*bx)k 
//
//		iii. we are interested in only zero-order (z) scattering 
//		so we will be interested in only the k component:
//
//			(a x b)z = (ax*by - ay*bx)k			
//
//
//========================================================================
//========================================================================
void calc_transm()
{

  // loop indices
  int i, ii, j, jj, m, p;


  //=========================================================
  // SUM FIELD COMPONENTS ON CALCULATION PLANES
  //=========================================================
  // i. we allocate so much memory because we have to send real and imaginary 
  // parts separately with MPI

  complex<double> ftexsum_transm, fteysum_transm, fthxsum_transm, fthysum_transm;
  complex<double> ftexsum_refl, fteysum_refl, fthxsum_refl, fthysum_refl;

  double *ftexsum_transm_real, *fteysum_transm_real, *fthxsum_transm_real, *fthysum_transm_real;
  double *ftexsum_transm_imag, *fteysum_transm_imag, *fthxsum_transm_imag, *fthysum_transm_imag;
  double *ftexsum_refl_real, *fteysum_refl_real, *fthxsum_refl_real, *fthysum_refl_real;
  double *ftexsum_refl_imag, *fteysum_refl_imag, *fthxsum_refl_imag, *fthysum_refl_imag;

  ftexsum_transm_real = new double [spect_npts+1];
  ftexsum_transm_imag = new double [spect_npts+1];

  fteysum_transm_real = new double [spect_npts+1];
  fteysum_transm_imag = new double [spect_npts+1];

  fthxsum_transm_real = new double [spect_npts+1];
  fthxsum_transm_imag = new double [spect_npts+1];

  fthysum_transm_real = new double [spect_npts+1];
  fthysum_transm_imag = new double [spect_npts+1];

  ftexsum_refl_real = new double [spect_npts+1];
  ftexsum_refl_imag = new double [spect_npts+1];

  fteysum_refl_real = new double [spect_npts+1];
  fteysum_refl_imag = new double [spect_npts+1];

  fthxsum_refl_real = new double [spect_npts+1];
  fthxsum_refl_imag = new double [spect_npts+1];

  fthysum_refl_real = new double [spect_npts+1];
  fthysum_refl_imag = new double [spect_npts+1];



  // do this for each transmission point
  for(m = 1; m <= spect_npts; ++m)
  {
    // zero the sums
    ftexsum_transm = COMPLEXZERO;
    fteysum_transm = COMPLEXZERO;
    fthxsum_transm = COMPLEXZERO;
    fthysum_transm = COMPLEXZERO;

    ftexsum_refl = COMPLEXZERO;
    fteysum_refl = COMPLEXZERO;
    fthxsum_refl = COMPLEXZERO;
    fthysum_refl = COMPLEXZERO;

    // INTEGRATE THE FIELDS OVER THE CALCULATION SURFACES
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;


        if( (ixperiodic == 0) && (iyperiodic == 0) )
        {

          if( (ii > cpml_layers) && (ii < grid_nx-cpml_layers) && (jj > cpml_layers) && (jj < grid_ny-cpml_layers) )
          {

            // transmitted fields
            ftexsum_transm += transm_ftex[m][i][j];
            fteysum_transm += transm_ftey[m][i][j];
            fthxsum_transm += transm_fthx[m][i][j];
            fthysum_transm += transm_fthy[m][i][j];

            // relected fields
            ftexsum_refl += refl_ftex[m][i][j];
            fteysum_refl += refl_ftey[m][i][j];
            fthxsum_refl += refl_fthx[m][i][j];
            fthysum_refl += refl_fthy[m][i][j];

          }


        }
        else 
        {
          // transmitted fields
          ftexsum_transm += transm_ftex[m][i][j];
          fteysum_transm += transm_ftey[m][i][j];
          fthxsum_transm += transm_fthx[m][i][j];
          fthysum_transm += transm_fthy[m][i][j];

          // relected fields
          ftexsum_refl += refl_ftex[m][i][j];
          fteysum_refl += refl_ftey[m][i][j];
          fthxsum_refl += refl_fthx[m][i][j];
          fthysum_refl += refl_fthy[m][i][j];
        }

      } // ++ii
    } // ++jj



    // ASSIGN THE FIELDS TO AN ARRAY

    // transmitted fields
    ftexsum_transm_real[m] = real(ftexsum_transm);
    ftexsum_transm_imag[m] = imag(ftexsum_transm);

    fteysum_transm_real[m] = real(fteysum_transm);
    fteysum_transm_imag[m] = imag(fteysum_transm);

    fthxsum_transm_real[m] = real(fthxsum_transm);
    fthxsum_transm_imag[m] = imag(fthxsum_transm);

    fthysum_transm_real[m] = real(fthysum_transm);
    fthysum_transm_imag[m] = imag(fthysum_transm);

    // reflected fields
    ftexsum_refl_real[m] = real(ftexsum_refl);
    ftexsum_refl_imag[m] = imag(ftexsum_refl);

    fteysum_refl_real[m] = real(fteysum_refl);
    fteysum_refl_imag[m] = imag(fteysum_refl);

    fthxsum_refl_real[m] = real(fthxsum_refl);
    fthxsum_refl_imag[m] = imag(fthxsum_refl);

    fthysum_refl_real[m] = real(fthysum_refl);
    fthysum_refl_imag[m] = imag(fthysum_refl);

  } // ++m


  //=========================================================
  // SEND ALL INFO TO MASTER PROCESS
  //=========================================================
  // i.) note we are using blocking send/recv here
  // ii.) we are also accessing the MPI_COMM_WORLD communicator IDs not the grid one

  int mpi_source, mpi_destination;
  MPI_Status mpi_status;
  int mpi_tag;

  // SEND VALUES TO MASTER PROCESSOR
  if(mpi_my_id != mpi_master_id)
  {

    // set the destination to the master process
    mpi_destination = mpi_master_id;

    // SEND THE TRANSMITTED FIELD

    // now send the ex (real)
    mpi_tag = 41;
    MPI_Send(ftexsum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the ex (imag)
    mpi_tag = 42;
    MPI_Send(ftexsum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the ey inc (real)
    mpi_tag = 43;
    MPI_Send(fteysum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the ey inc (imag)
    mpi_tag = 44;
    MPI_Send(fteysum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the hx inc (real)
    mpi_tag = 45;
    MPI_Send(fthxsum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the hx inc (imag)
    mpi_tag = 46;
    MPI_Send(fthxsum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the hy inc (real)
    mpi_tag = 47;
    MPI_Send(fthysum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the hy inc (imag)
    mpi_tag = 48;
    MPI_Send(fthysum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // SEND THE REFLECTED FIELD

    // now send the ex (real)
    mpi_tag = 51;
    MPI_Send(ftexsum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the ex (imag)
    mpi_tag = 52;
    MPI_Send(ftexsum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the ey inc (real)
    mpi_tag = 53;
    MPI_Send(fteysum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the ey inc (imag)
    mpi_tag = 54;
    MPI_Send(fteysum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the hx inc (real)
    mpi_tag = 55;
    MPI_Send(fthxsum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the hx inc (imag)
    mpi_tag = 56;
    MPI_Send(fthxsum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);

    // now send the hy inc (real)
    mpi_tag = 57;
    MPI_Send(fthysum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
    // now send the hy inc (imag)
    mpi_tag = 58;
    MPI_Send(fthysum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
  }



  //=========================================================
  // RECEIVE INFO AT MASTER PROCESS AND CALCULATE PROPERTIES
  //=========================================================
  // TW l.1160 shows pzsums defined


  if(mpi_my_id == mpi_master_id)
  {

    //---------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------

    // fields
    complex<double> *mast_ftexsum_transm, *mast_fteysum_transm, *mast_fthxsum_transm, *mast_fthysum_transm;
    complex<double> *mast_ftexsum_refl, *mast_fteysum_refl, *mast_fthxsum_refl, *mast_fthysum_refl;

    mast_ftexsum_transm = new complex<double> [spect_npts+1];
    mast_fteysum_transm = new complex<double> [spect_npts+1];
    mast_fthxsum_transm = new complex<double> [spect_npts+1];
    mast_fthysum_transm = new complex<double> [spect_npts+1];

    mast_ftexsum_refl = new complex<double> [spect_npts+1];
    mast_fteysum_refl = new complex<double> [spect_npts+1];
    mast_fthxsum_refl = new complex<double> [spect_npts+1];
    mast_fthysum_refl = new complex<double> [spect_npts+1];

    // Poynting vectors
    // i. we will calculate the Poynting vectors in 2 ways
    // ii. how I use pzsum_abs is not exactly correct as it
    // signifies just the absorption, not a Poynting vector
    double *pzsum_transm, *pzsum_refl, *pzsum_abs;

    pzsum_transm = new double [spect_npts+1];
    pzsum_refl = new double [spect_npts+1];
    pzsum_abs = new double [spect_npts+1];


    // ZERO ARRAYS
    for(m = 1; m <= spect_npts; ++m)
    {
      // transmitted fields
      mast_ftexsum_transm[m] = COMPLEXZERO;
      mast_fteysum_transm[m] = COMPLEXZERO;
      mast_fthxsum_transm[m] = COMPLEXZERO;
      mast_fthysum_transm[m] = COMPLEXZERO;

      // reflected fields
      mast_ftexsum_refl[m] = COMPLEXZERO;
      mast_fteysum_refl[m] = COMPLEXZERO;
      mast_fthxsum_refl[m] = COMPLEXZERO;
      mast_fthysum_refl[m] = COMPLEXZERO;

      // Poynting vectors
      pzsum_transm[m] = 0.0;
      pzsum_refl[m] = 0.0;
    }

    //---------------------------------------------------------
    // SUM FIELDS
    //---------------------------------------------------------

    // FIRST ADD IN MASTER PROCESS VALUES
    for(m = 1; m <= spect_npts; ++m)
    {
      // transmitted fields
      mast_ftexsum_transm[m] = complex<double> (ftexsum_transm_real[m], ftexsum_transm_imag[m]);
      mast_fteysum_transm[m] = complex<double> (fteysum_transm_real[m], fteysum_transm_imag[m]);
      mast_fthxsum_transm[m] = complex<double> (fthxsum_transm_real[m], fthxsum_transm_imag[m]);
      mast_fthysum_transm[m] = complex<double> (fthysum_transm_real[m], fthysum_transm_imag[m]);

      // reflected fields
      mast_ftexsum_refl[m] = complex<double> (ftexsum_refl_real[m], ftexsum_refl_imag[m]);
      mast_fteysum_refl[m] = complex<double> (fteysum_refl_real[m], fteysum_refl_imag[m]);
      mast_fthxsum_refl[m] = complex<double> (fthxsum_refl_real[m], fthxsum_refl_imag[m]);
      mast_fthysum_refl[m] = complex<double> (fthysum_refl_real[m], fthysum_refl_imag[m]);
    }

    // NOW RECEIVE AND SUM VALUES
    for(p = 1; p <= (mpi_nprocessors - 1); ++p) 
    {
      // SET THE SOURCE TO THE CURRENT PROCESS(OR)
      mpi_source = p;
 
      // RECEIVE TRANSMITTED FIELDS

      // now receive the ex inc (real)
      mpi_tag = 41;
      MPI_Recv(ftexsum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the ex inc (imag)
      mpi_tag = 42;
      MPI_Recv(ftexsum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the ey inc (real)
      mpi_tag = 43;
      MPI_Recv(fteysum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the ey inc (imag)
      mpi_tag = 44;
      MPI_Recv(fteysum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the hx inc (real)
      mpi_tag = 45;
      MPI_Recv(fthxsum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the hx inc (imag)
      mpi_tag = 46;
      MPI_Recv(fthxsum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the hy inc (real)
      mpi_tag = 47;
      MPI_Recv(fthysum_transm_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the hy inc (imag)
      mpi_tag = 48;
      MPI_Recv(fthysum_transm_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);


      // RECEIVE REFLECTED FIELDS

      // now receive the ex inc (real)
      mpi_tag = 51;
      MPI_Recv(ftexsum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the ex inc (imag)
      mpi_tag = 52;
      MPI_Recv(ftexsum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the ey inc (real)
      mpi_tag = 53;
      MPI_Recv(fteysum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the ey inc (imag)
      mpi_tag = 54;
      MPI_Recv(fteysum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the hx inc (real)
      mpi_tag = 55;
      MPI_Recv(fthxsum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the hx inc (imag)
      mpi_tag = 56;
      MPI_Recv(fthxsum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // now receive the hy inc (real)
      mpi_tag = 57;
      MPI_Recv(fthysum_refl_real, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);
      // now receive the hy inc (imag)
      mpi_tag = 58;
      MPI_Recv(fthysum_refl_imag, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);



      // ADD IN THE RECEIVED VALUES TO OUR SUMMATION
      for(m = 1; m <= spect_npts; ++m)
      {
        // transmitted fields
        mast_ftexsum_transm[m] += complex<double> (ftexsum_transm_real[m], ftexsum_transm_imag[m]);
        mast_fteysum_transm[m] += complex<double> (fteysum_transm_real[m], fteysum_transm_imag[m]);
        mast_fthxsum_transm[m] += complex<double> (fthxsum_transm_real[m], fthxsum_transm_imag[m]);
        mast_fthysum_transm[m] += complex<double> (fthysum_transm_real[m], fthysum_transm_imag[m]);

        // reflected fields
        mast_ftexsum_refl[m] += complex<double> (ftexsum_refl_real[m], ftexsum_refl_imag[m]);
        mast_fteysum_refl[m] += complex<double> (fteysum_refl_real[m], fteysum_refl_imag[m]);
        mast_fthxsum_refl[m] += complex<double> (fthxsum_refl_real[m], fthxsum_refl_imag[m]);
        mast_fthysum_refl[m] += complex<double> (fthysum_refl_real[m], fthysum_refl_imag[m]);
      } // ++m

    }  // ++p

 
    //---------------------------------------------------------
    // CALCULATE & NORMALIZE POYNTING VECTORS
    //---------------------------------------------------------
    // i. remember: (a x b)z = (ax*by - ay*bx)k	
    // TW. L. 1141

    // general Poynting vector
    complex<double> pzval;

    double tmp0, tmp1;
    complex<double> cabs_einc;
    double abs_einc;

    if( (ixperiodic == 0) && (iyperiodic == 0) )
    {

      // i. the easiest thing to do here would be to put all of the grid_xx inside
      // a single double conversion, but this could possibly lead to an integer overflow
      // (this happened in Tae-Woo's code, but everything was left as integer)
      tmp0 = (static_cast<double>(grid_nx-cpml_layers))*(static_cast<double>(grid_ny-cpml_layers))*(static_cast<double>(grid_nx-cpml_layers))*(static_cast<double>(grid_ny-cpml_layers)); 

    }
    else
    {

      // i. the easiest thing to do here would be to put all of the grid_xx inside
      // a single double conversion, but this could possibly lead to an integer overflow
      // (this happened in Tae-Woo's code, but everything was left as integer)
      tmp0 = (static_cast<double>(grid_nx))*(static_cast<double>(grid_ny))*(static_cast<double>(grid_nx))*(static_cast<double>(grid_ny));

    }

// j12
    // loop over the transmission points
    for(m = 1; m <= spect_npts; ++m)
    {

      // TRANSMISSION
      // i. this is give by: Pz = (1/2)*Re[E_tot x H_tot*]
      // ii. since light is incident from bottom -> top we want the +z component
      // of the Poynting vector for transmission
      pzval = mast_ftexsum_transm[m]*conj(mast_fthysum_transm[m]) - mast_fteysum_transm[m]*conj(mast_fthxsum_transm[m]);
      pzsum_transm[m] = 0.5*real(pzval);

      // REFLECTION
      // i. this is give by: Pz = -(1/2)*Re[E_tot x H_tot*]
      // ii. since light is incident from bottom -> top we want the -z component
      // of the Poynting vector for reflection
      pzval = 0.0 - (mast_ftexsum_refl[m]*conj(mast_fthysum_refl[m]) - mast_fteysum_refl[m]*conj(mast_fthxsum_refl[m]));
      pzsum_refl[m] = 0.5*real(pzval);


      // NORMALIZE BY INCIDENT POWER

      cabs_einc = aux_ftex_transm[m]*conj(aux_ftex_transm[m]) + aux_ftey_transm[m]*conj(aux_ftey_transm[m]) + aux_ftez_transm[m]*conj(aux_ftez_transm[m]);
      // !!! should this sqrt actually be here
      abs_einc = sqrt(real(cabs_einc));

      // now get tmp1
      // from TW code: tmp1=cabs(eincdft(m))*cabs(eincdft(m))/imped0/2.0*tmp0
      tmp1 = abs_einc*abs_einc*tmp0/(Z0*2.0);
      // i. !!! the sqrt below should be the epsilon of where the source is located
      tmp1 *= sqrt(aux_epsr);

      // normalize zero-order Poynting vectors
      // i. both the transmission and reflection are normalized in the same way
      pzsum_transm[m] /= tmp1;
      pzsum_refl[m] /= tmp1;

      // we can now use the relation [1 = T + R + A] to get the absorption
      // ii. even though it says pzsum this is absorption, not a Poynting vector
      pzsum_abs[m] = 1.0 - pzsum_transm[m] - pzsum_refl[m];

    } // ++m


    //---------------------------------------------------------
    // WRITE FILES
    //---------------------------------------------------------
    // i.) OUTPUT TRANSMITTANCE AND REFLECTANCE

    // OPEN FILES
    ofstream spect_transm_file("./output/spect_transm.dat",ios::app);
    ofstream spect_refl_file("./output/spect_refl.dat",ios::app);
    ofstream spect_abs_file("./output/spect_abs.dat",ios::app);


    // WRITE FILES
    double ft_wave;
    

    for(m = 1; m <= spect_npts; ++m)
    {
      // get the wavelength at this wavevector
      ft_wave = 2.0*PI*CSPEED/spect_wavenumbers[m];

      // write the normalized values
      spect_transm_file << ft_wave << "     " << pzsum_transm[m] << endl;
      spect_refl_file << ft_wave << "     " << pzsum_refl[m] << endl;
      spect_abs_file << ft_wave << "     " << pzsum_abs[m] << endl;
    }


    // CLOSE FILES
    spect_transm_file.close();
    spect_refl_file.close();
    spect_abs_file.close();


    //---------------------------------------------------------
    // CLEANUP
    //---------------------------------------------------------

    // delete sums
    delete [] mast_ftexsum_transm;
    delete [] mast_fthysum_transm;
    delete [] mast_fteysum_transm;
    delete [] mast_fthxsum_transm;

    delete [] mast_ftexsum_refl;
    delete [] mast_fthysum_refl;
    delete [] mast_fteysum_refl;
    delete [] mast_fthxsum_refl;


    // delete Poynting vectors
    delete [] pzsum_transm;
    delete [] pzsum_refl;
    delete [] pzsum_abs;

  } // end if master id




  //=========================================================
  // CLEANUP
  //=========================================================

  // delete the zero-order FT (components)
  delete [] ftexsum_transm_real;
  delete [] ftexsum_transm_imag;

  delete [] fteysum_transm_real;
  delete [] fteysum_transm_imag;

  delete [] fthxsum_transm_real;
  delete [] fthxsum_transm_imag;

  delete [] fthysum_transm_real;
  delete [] fthysum_transm_imag;

  delete [] ftexsum_refl_real;
  delete [] ftexsum_refl_imag;

  delete [] fteysum_refl_real;
  delete [] fteysum_refl_imag;

  delete [] fthxsum_refl_real;
  delete [] fthxsum_refl_imag;

  delete [] fthysum_refl_real;
  delete [] fthysum_refl_imag;


  return;

}



//========================================================================
//========================================================================
//
//	NAME:	void calc_transm_nf()
//	DESC:	calculates the transmittance spectrum in the near field
//
//	NOTES: 	i. this is called by all processes
//
//		ii. the cross product is defined as:
//
//			a x b = (ay*bz - az*by)i + 
//				(az*bx - ax*bz)j + 
//				(ax*by - ay*bx)k 
//
//		iii. we are interested in only zero-order (z) scattering 
//		so we will be interested in only the k component:
//
//			(a x b)z = (ax*by - ay*bx)k
//
//		iv. ! remember that the aperture size is the diameter of the tip
//		so we divide by 2 because we compare radius' below			
//
//
//========================================================================
//========================================================================
void calc_transm_nf()
{

  // loop indices
  int i, ii, j, jj, m, p;

  //=========================================================
  // INTEGRATE POYNTING VECTOR IF IT IS IN THE APERTURE
  //=========================================================

  double xpos, ypos, xposp, yposp, rad;
  int inside_ap;

  complex<double> pzsum;

  double *pzsum_real;

  pzsum_real = new double [spect_npts+1];


  // do this for each transmission point
  for(m = 1; m <= spect_npts; ++m)
  {
    // zero the sums
    pzsum = COMPLEXZERO;

    // INTEGRATE THE FIELDS OVER THE CALCULATION SURFACES
    for(jj = mpi_my_grid_ny1; jj <= mpi_my_grid_ny2; ++jj)
    {
      j = jj - mpi_my_grid_ny1 + 1;

      // GET THE y-POSITION (CENTERED)
      ypos = (static_cast<double>(jj))*grid_dy;


      for(ii = mpi_my_grid_nx1; ii <= mpi_my_grid_nx2; ++ii)
      {
        i = ii - mpi_my_grid_nx1 + 1;

        // GET THE x-POSITION (CENTERED)
        xpos = (static_cast<double>(ii))*grid_dx;

        
        // FIND OUT IF THE POSITION IS INSIDE THE NEAR FIELD APERTURE
        // i. due to periodicity we need to check all permutations of the periodicity (there are
        // 9 of them) to find out if the aperture exists there

        // set a flag to not inside aperture (this will be switched if one of our permutations lies
        // inside the aperture)
        inside_ap = 0;

        // 1
        yposp = ypos;
        xposp = xpos;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 2
        yposp = ypos + grid_ysize;
        xposp = xpos;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 3
        yposp = ypos - grid_ysize;
        xposp = xpos;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 4
        yposp = ypos;
        xposp = xpos + grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 5
        yposp = ypos;
        xposp = xpos - grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 6
        yposp = ypos - grid_ysize;
        xposp = xpos - grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 7
        yposp = ypos + grid_ysize;
        xposp = xpos - grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 8
        yposp = ypos - grid_ysize;
        xposp = xpos + grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }

        // 9
        yposp = ypos + grid_ysize;
        xposp = xpos + grid_xsize;
        rad = sqrt( (xposp-nf_aperture_xcenter)*(xposp-nf_aperture_xcenter) +  (yposp-nf_aperture_ycenter)*(yposp-nf_aperture_ycenter));

        if(rad <= nf_aperture_size/2.0)
        {
          inside_ap = 1;
        }


        // if never ended up in the aperture from any permutation then continue
        if(inside_ap == 0)
        {
          continue;
        }


        // INTEGRATE THE POYNTING VECTOR
        pzsum += transm_ftex[m][i][j]*conj(transm_fthy[m][i][j]) - transm_ftey[m][i][j]*conj(transm_fthx[m][i][j]);

      } // ++ii
    } // ++jj


    // ASSIGN THE FIELDS TO AN ARRAY
    pzsum_real[m] = 0.5*real(pzsum);

  } // ++m



  //=========================================================
  // SEND ALL INFO TO MASTER PROCESS
  //=========================================================
  // i.) note we are using blocking send/recv here
  // ii.) we are also accessing the MPI_COMM_WORLD communicator IDs not the grid one

  int mpi_source, mpi_destination;
  MPI_Status mpi_status;
  int mpi_tag;

  // SEND VALUES TO MASTER PROCESSOR
  if(mpi_my_id != mpi_master_id)
  {

    // set the destination to the master process
    mpi_destination = mpi_master_id;

    // SEND THE TRANSMITTED FIELD

    // now send the ex (real)
    mpi_tag = 321;
    MPI_Send(pzsum_real, (spect_npts+1), MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
  }



  //=========================================================
  // RECEIVE INFO AT MASTER PROCESS AND CALCULATE PROPERTIES
  //=========================================================
  // TW l.1160 shows pzsums defined


  if(mpi_my_id == mpi_master_id)
  {

    //---------------------------------------------------------
    // INITIALIZE
    //---------------------------------------------------------
    double *pz_recv;

    pz_recv = new double [spect_npts+1];


    // NOW RECEIVE AND SUM VALUES
    for(p = 1; p <= (mpi_nprocessors - 1); ++p) 
    {
      // SET THE SOURCE TO THE CURRENT PROCESS(OR)
      mpi_source = p;
 
      // RECEIVE TRANSMITTED FIELDS

      // now receive the ex inc (real)
      mpi_tag = 321;
      MPI_Recv(pz_recv, (spect_npts+1), MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);


      // ADD IN THE RECEIVED VALUES TO OUR SUMMATION
      for(m = 1; m <= spect_npts; ++m)
      {
        // transmitted fields
        pzsum_real[m] += pz_recv[m]; 
      } // ++m

    }  // ++p




    //---------------------------------------------------------
    // NORMALIZE BY THE INCIDENT POYNTING VECTOR
    //---------------------------------------------------------
    // i. remember: (a x b)z = (ax*by - ay*bx)k	 

    // GET THE INCIDENT POYNTING VECTOR
    complex<double> pzinc;
    double pzinc_real;

    // loop over the transmission points
    for(m = 1; m <= spect_npts; ++m)
    {

      pzinc = aux_ftex_transm[m]*conj(aux_fthy_transm[m]) - aux_ftey_transm[m]*conj(aux_fthx_transm[m]);
      pzinc_real = 0.5*real(pzinc);

      // NORMALIZE THE NORMAL VALUE
      pzsum_real[m] /= pzinc_real;

    } // ++m


    //---------------------------------------------------------
    // WRITE FILES
    //---------------------------------------------------------
    // i.) OUTPUT TRANSMITTANCE AND REFLECTANCE

    // OPEN FILES
    ofstream spect_transm_file("./output/transm_nf.dat",ios::app);

    // WRITE FILES
    double ft_wave;
    
    for(m = 1; m <= spect_npts; ++m)
    {
      // get the wavelength at this wavevector
      ft_wave = 2.0*PI*CSPEED/spect_wavenumbers[m];

      // write the normalized values
      spect_transm_file << ft_wave << "     " << pzsum_real[m] << endl;
    }


    // CLOSE FILES
    spect_transm_file.close();


    //---------------------------------------------------------
    // CLEANUP
    //---------------------------------------------------------

    // delete Poynting vectors
    delete [] pz_recv;

  } // end if master id


  //=========================================================
  // CLEANUP
  //=========================================================

  delete [] pzsum_real;


  return;

}




//========================================================================
//========================================================================
//
//	NAME:	void calc_transm()
//	DESC:	calculates the transmittance spectrum
//
//	NOTES: 	i. this is called by all processes
//
//
//========================================================================
//========================================================================
void calc_transm_nsom()
{

  // local indices
  int n, p;


  //=========================================================
  // GET THE NSOM POYNTING VECTOR CALCULATED BY THIS PROCESS
  //=========================================================

  complex<double> pzsum = COMPLEXZERO;
  double pzsum_real;

  // INTEGRATE THE POYNTING VECTOR OVER THE ENTIRE NSOM PROBE
  for(n = 1; n <= transm_nsom_npts; ++n)
  {
    pzsum += transm_nsom_ex[n]*conj(transm_nsom_hy[n]) - transm_nsom_ey[n]*conj(transm_nsom_hx[n]);
  }

  // ASSIGN THE REAL POYNTING VECTOR TO THE COMPLEX ONE AND MULTIPLY BY 0.5
  pzsum_real = 0.5*real(pzsum);


  //=========================================================
  // SEND SUMS TO MASTER PROCESS
  //=========================================================
  // i.) note we are using blocking send/recv here
  // ii.) we are also accessing the MPI_COMM_WORLD communicator IDs not the grid one

  int mpi_source, mpi_destination;
  MPI_Status mpi_status;
  int mpi_tag;

  // SEND VALUES TO MASTER PROCESSOR
  if(mpi_my_id != mpi_master_id)
  {

    // set the destination to the master process
    mpi_destination = mpi_master_id;

    // send the real part of the poynting vector
    mpi_tag = 151;
    MPI_Send((void *)&pzsum_real, 1, MPI_DOUBLE, mpi_destination, mpi_tag, MPI_COMM_WORLD);
  }



  //=========================================================
  // RECEIVE POYNTING VECTOR AT MASTER PROCESS AND WRITE FILE
  //=========================================================

  if(mpi_my_id == mpi_master_id)
  {

    double pz_recv;

    // NOW RECEIVE AND SUM VALUES
    for(p = 1; p <= (mpi_nprocessors - 1); ++p) 
    {
      // SET THE SOURCE TO THE CURRENT PROCESS(OR)
      mpi_source = p;
 
      // receive the real part of the poynting vector
      mpi_tag = 151;
      MPI_Recv((void *)&pz_recv, 1, MPI_DOUBLE, mpi_source, mpi_tag, MPI_COMM_WORLD, &mpi_status);

      // ADD IN THE RECEIVED VALUE TO OUR SUMMATION
      pzsum_real += pz_recv;
    }  // ++p


    //---------------------------------------------------------
    // NORMALIZE BY THE INCIDENT POYNTING VECTOR
    //---------------------------------------------------------
    // i. remember: (a x b)z = (ax*by - ay*bx)k	 

    // GET THE INCIDENT POYNTING VECTOR
    complex<double> pzinc = aux_ftex_nsom*conj(aux_fthy_nsom) - aux_ftey_nsom*conj(aux_fthx_nsom);
    double pzinc_real = 0.5*real(pzinc);

    // NORMALIZE THE NORMAL VALUE
    pzsum_real /= pzinc_real;



    //---------------------------------------------------------
    // WRITE FILE
    //---------------------------------------------------------

    // OPEN FILES
    ofstream nsom_file("./output/nsom.dat",ios::app);

    // WRITE
    nsom_file << "nsom_xcenter: " << nsom_xcenter << endl;
    nsom_file << "nsom_ycenter: " << nsom_ycenter << endl;
    nsom_file << "nsom_zcenter: " << nsom_zcenter << endl;
    nsom_file << "nsom wavelength: " << nsom_wavelength << endl;
    nsom_file << "Pz sum: " << pzsum_real << endl;

    // CLOSE FILE
    nsom_file.close();

  } // end if master id


  return;

}



//========================================================================
//========================================================================
//
//	NAME:	void free_memory()
//	DESC:	Releases all allocated (global) dynamic memory.
//
//	NOTES: 	i. 
//
//
//========================================================================
//========================================================================
void free_memory()
{

  // "local" indices
  int i, j, l, m, p;



  // OUTPUT
  delete [] output_format;
  delete [] output_field;
  delete [] output_plane;
  delete [] output_plane_midpos;
  delete [] ft_wavelengths;


  delete [] ft_wavenumbers;
  delete [] ft_aux_coord;
  delete [] output_to_ft;
  delete [] ft_to_output;
  delete [] mpi_ioutput_plane;
  delete [] output_plane_local_coord;

  if(nft_planes!=0)
  {
    delete [] aux_ftex;
    delete [] aux_ftey;
    delete [] aux_ftez;

    delete [] aux_fthx;
    delete [] aux_fthy;
    delete [] aux_fthz;
  }


  for(p=0;p<nft_planes;++p)
  {
    // ----------------------- xy ------------------------------
    if(output_plane[ft_to_output[p]] == 1)
    {
      for(i=0;i<=my_grid_nx+1; ++i)
      {
        delete [] field_ft_ex[p][i];
        delete [] field_ft_ey[p][i];
        delete [] field_ft_ez[p][i];
        delete [] field_ft_hx[p][i];
        delete [] field_ft_hy[p][i];
        delete [] field_ft_hz[p][i];
      }
      delete [] field_ft_ex[p];
      delete [] field_ft_ey[p];
      delete [] field_ft_ez[p];
      delete [] field_ft_hx[p];
      delete [] field_ft_hy[p];
      delete [] field_ft_hz[p];
    }
    // ----------------------- xz ------------------------------
    else if(output_plane[ft_to_output[p]] == 2)
    {
      for(i=0;i<=my_grid_nx+1; ++i)
      {
        delete [] field_ft_ex[p][i];
        delete [] field_ft_ey[p][i];
        delete [] field_ft_ez[p][i];
        delete [] field_ft_hx[p][i];
        delete [] field_ft_hy[p][i];
        delete [] field_ft_hz[p][i];
      }
      delete [] field_ft_ex[p];
      delete [] field_ft_ey[p];
      delete [] field_ft_ez[p];
      delete [] field_ft_hx[p];
      delete [] field_ft_hy[p];
      delete [] field_ft_hz[p];
    }
    // ----------------------- yz ------------------------------
    else if(output_plane[ft_to_output[p]] == 2)
    {
      for(j=0;j<=my_grid_ny+1; ++j)
      {
        delete [] field_ft_ex[p][j];
        delete [] field_ft_ey[p][j];
        delete [] field_ft_ez[p][j];
        delete [] field_ft_hx[p][j];
        delete [] field_ft_hy[p][j];
        delete [] field_ft_hz[p][j];
      }
      delete [] field_ft_ex[p];
      delete [] field_ft_ey[p];
      delete [] field_ft_ez[p];
      delete [] field_ft_hx[p];
      delete [] field_ft_hy[p];
      delete [] field_ft_hz[p];
    }

  } //++p

  if(nft_planes!=0)
  {
    delete [] field_ft_ex;
    delete [] field_ft_ey;
    delete [] field_ft_ez;
    delete [] field_ft_hx;
    delete [] field_ft_hy;
    delete [] field_ft_hz;
  }

  // OK LINE ***********************************************************


  //=========================================================
  // USER (OR RELATED) PARAMETERS
  //=========================================================


  for(i = 0; i < d2l_libsize; ++i)
  {
    delete [] lorentz_alphap[i];
    delete [] lorentz_zetap[i];
    delete [] lorentz_gammap[i];
  }

  delete [] lorentz_alphap;
  delete [] lorentz_zetap;
  delete [] lorentz_gammap;

  //=========================================================
  // MPI
  //=========================================================

  //---------------------------------------------------------
  // COMMUNICATION
  //---------------------------------------------------------

  delete [] ksurface_x_send;
  delete [] ksurface_y_send;

  delete [] jsurface_x_send;
  delete [] jsurface_z_send;

  delete [] isurface_y_send;
  delete [] isurface_z_send;

  delete [] isurface_y_recv;
  delete [] isurface_z_recv;

  delete [] jsurface_x_recv;
  delete [] jsurface_z_recv;

  delete [] ksurface_x_recv;
  delete [] ksurface_y_recv;


  //=========================================================
  // MAIN FDTD
  //=========================================================

  //---------------------------------------------------------
  // FIELDS & CURRENTS
  //---------------------------------------------------------
  for(i=0;i<=my_grid_nx+1;++i)
  {
    for(j=0;j<=my_grid_ny+1;++j)
    {
      delete [] field_ex[i][j];
      delete [] field_ey[i][j];
      delete [] field_ez[i][j];
      delete [] field_ex_nm1[i][j];
      delete [] field_ey_nm1[i][j];
      delete [] field_ez_nm1[i][j];
      delete [] field_hx[i][j];
      delete [] field_hy[i][j];
      delete [] field_hz[i][j];

      delete [] grid_material_ex[i][j];
      delete [] grid_material_ey[i][j];
      delete [] grid_material_ez[i][j];

      delete [] permittivity_ex[i][j];
      delete [] permittivity_ey[i][j];
      delete [] permittivity_ez[i][j];
      delete [] permeability_hx[i][j];
      delete [] permeability_hy[i][j];
      delete [] permeability_hz[i][j];

      delete [] drude_current_x[i][j];
      delete [] drude_current_y[i][j];
      delete [] drude_current_z[i][j];
    }
    delete [] field_ex[i];
    delete [] field_ey[i];
    delete [] field_ez[i];
    delete [] field_ex_nm1[i];
    delete [] field_ey_nm1[i];
    delete [] field_ez_nm1[i];
    delete [] field_hx[i];
    delete [] field_hy[i];
    delete [] field_hz[i];

    delete [] grid_material_ex[i];
    delete [] grid_material_ey[i];
    delete [] grid_material_ez[i];

    delete [] permittivity_ex[i];
    delete [] permittivity_ey[i];
    delete [] permittivity_ez[i];
    delete [] permeability_hx[i];
    delete [] permeability_hy[i];
    delete [] permeability_hz[i];

    delete [] drude_current_x[i];
    delete [] drude_current_y[i];
    delete [] drude_current_z[i];
  }

  delete [] field_ex;
  delete [] field_ey;
  delete [] field_ez;
  delete [] field_ex_nm1;
  delete [] field_ey_nm1;
  delete [] field_ez_nm1;
  delete [] field_hx;
  delete [] field_hy;
  delete [] field_hz;

  delete [] grid_material_ex;
  delete [] grid_material_ey;
  delete [] grid_material_ez;

  delete [] permittivity_ex;
  delete [] permittivity_ey;
  delete [] permittivity_ez;
  delete [] permeability_hx;
  delete [] permeability_hy;
  delete [] permeability_hz;
  
  delete [] drude_current_x;
  delete [] drude_current_y;
  delete [] drude_current_z;

  for(l=0;l<2;++l)
  {
    for(i=0;i<=my_grid_nx+1;++i)
    {
      for(j=0;j<=my_grid_ny+1;++j)
      {
        delete [] lorentz_currentp_x[l][i][j];
        delete [] lorentz_currentp_y[l][i][j];
        delete [] lorentz_currentp_z[l][i][j];
        delete [] lorentz_currentp_x_nm1[l][i][j];
        delete [] lorentz_currentp_y_nm1[l][i][j];
        delete [] lorentz_currentp_z_nm1[l][i][j];
      }

      delete [] lorentz_currentp_x[l][i];
      delete [] lorentz_currentp_y[l][i];
      delete [] lorentz_currentp_z[l][i];
      delete [] lorentz_currentp_x_nm1[l][i];
      delete [] lorentz_currentp_y_nm1[l][i];
      delete [] lorentz_currentp_z_nm1[l][i];
    }

    delete [] lorentz_currentp_x[l];
    delete [] lorentz_currentp_y[l];
    delete [] lorentz_currentp_z[l];
    delete [] lorentz_currentp_x_nm1[l];
    delete [] lorentz_currentp_y_nm1[l];
    delete [] lorentz_currentp_z_nm1[l];
  }

  delete [] lorentz_currentp_x;
  delete [] lorentz_currentp_y;
  delete [] lorentz_currentp_z;
  delete [] lorentz_currentp_x_nm1;
  delete [] lorentz_currentp_y_nm1;
  delete [] lorentz_currentp_z_nm1;

  //---------------------------------------------------------
  // LOSSY MATERIAL TERMS
  //---------------------------------------------------------
  delete [] drude_kp;
  delete [] drude_betap;

  delete [] metal_c1;
  delete [] metal_c2;
  delete [] metal_c3;

  //---------------------------------------------------------
  // PML
  //---------------------------------------------------------
  delete [] kedx;
  delete [] kedy;
  delete [] kedz;
  delete [] khdx;
  delete [] khdy;
  delete [] khdz;

  delete [] be_x;
  delete [] ce_x;
  delete [] bh_x;
  delete [] ch_x;

  delete [] be_y;
  delete [] ce_y;
  delete [] bh_y;
  delete [] ch_y;

  delete [] be_z;
  delete [] ce_z;
  delete [] bh_z;
  delete [] ch_z;

  for(i = 0; i <= (my_grid_nx+1); ++i)
  {
    for(j = 0; j <= (my_grid_ny+1); ++j)
    {
      delete [] psi_eyx[i][j];
      delete [] psi_ezx[i][j];
      delete [] psi_hyx[i][j];
      delete [] psi_hzx[i][j];

      delete [] psi_exy[i][j];
      delete [] psi_ezy[i][j];
      delete [] psi_hxy[i][j];
      delete [] psi_hzy[i][j];

      delete [] psi_exz[i][j];
      delete [] psi_eyz[i][j];
      delete [] psi_hxz[i][j];
      delete [] psi_hyz[i][j];
    } 
    delete [] psi_eyx[i];
    delete [] psi_ezx[i];
    delete [] psi_hyx[i];
    delete [] psi_hzx[i];

    delete [] psi_exy[i];
    delete [] psi_ezy[i];
    delete [] psi_hxy[i];
    delete [] psi_hzy[i];

    delete [] psi_exz[i];
    delete [] psi_eyz[i];
    delete [] psi_hxz[i];
    delete [] psi_hyz[i];
  } 

  delete [] psi_eyx;
  delete [] psi_ezx;
  delete [] psi_hyx;
  delete [] psi_hzx;

  delete [] psi_exy;
  delete [] psi_ezy;
  delete [] psi_hxy;
  delete [] psi_hzy;

  delete [] psi_exz;
  delete [] psi_eyz;
  delete [] psi_hxz;
  delete [] psi_hyz;

  //=========================================================
  // AUXILLARY GRIDS
  //=========================================================
  delete [] ex_aux;
  delete [] ey_aux;
  delete [] hx_aux;
  delete [] hy_aux;

  //=========================================================
  // POST-PROCESSING
  //=========================================================

  //---------------------------------------------------------
  // TRANSMISSION
  //---------------------------------------------------------

  if((itransm==1) || (iscat_calc==1))
  {

    // FOURIER-TRANSFORM ARRAYS
    delete [] spect_wavenumbers;

    for(m = 0; m <= spect_npts; ++m)
    {
      delete [] dft_multiplier_t[m];
      delete [] dft_multiplier_ht[m];
    }

    delete [] dft_multiplier_t;
    delete [] dft_multiplier_ht;


    // AUXILLARY GRID TRANSMISSION/REFLECTION
    delete [] aux_ftex_transm;
    delete [] aux_ftey_transm;
    delete [] aux_ftez_transm;
    delete [] aux_fthx_transm;
    delete [] aux_fthy_transm;
    delete [] aux_fthz_transm;
  }

  if(itransm==1)
  {
    for(m=0;m<=spect_npts;++m)
    {
      for(i=0;i<=my_grid_nx+1;++i)
      {
        delete [] transm_ftex[m][i];
        delete [] transm_ftey[m][i];
        delete [] transm_fthx[m][i];
        delete [] transm_fthy[m][i];
        delete [] refl_ftex[m][i];
        delete [] refl_ftey[m][i];
        delete [] refl_fthx[m][i];
        delete [] refl_fthy[m][i];
      }
      delete [] transm_ftex[m];
      delete [] transm_ftey[m];
      delete [] transm_fthx[m];
      delete [] transm_fthy[m];
      delete [] refl_ftex[m];
      delete [] refl_ftey[m];
      delete [] refl_fthx[m];
      delete [] refl_fthy[m];
    }
    delete [] transm_ftex;
    delete [] transm_ftey;
    delete [] transm_fthx;
    delete [] transm_fthy;
    delete [] refl_ftex;
    delete [] refl_ftey;
    delete [] refl_fthx;
    delete [] refl_fthy;

    // AUXILLARY GRID TRANSMISSION/REFLECTION
    //delete [] aux_ftex_transm;
    //delete [] aux_ftey_transm;
    //delete [] aux_ftez_transm;
    //delete [] aux_fthx_transm;
    //delete [] aux_fthy_transm;
    //delete [] aux_fthz_transm;
    delete [] aux_ftex_refl;
    delete [] aux_ftey_refl;
    delete [] aux_ftez_refl;
    delete [] aux_fthx_refl;
    delete [] aux_fthy_refl;
    delete [] aux_fthz_refl;
  }


  //---------------------------------------------------------
  // SCATTERING CROSS SECTION
  //---------------------------------------------------------

  if(iscat_calc==1)
  {

    for(m = 0; m <= spect_npts; ++m)
    {
      delete [] scat_aux_ex[m];
      delete [] scat_aux_ey[m];
      delete [] scat_aux_hx[m];
      delete [] scat_aux_hy[m];
// j5
      for(i = 0; i <= my_grid_nx; ++i)
      {
        delete [] scat_k1_ftex[m][i];
        delete [] scat_k1_ftey[m][i];
        delete [] scat_k1_fthx[m][i];
        delete [] scat_k1_fthy[m][i];
        delete [] scat_k2_ftex[m][i];
        delete [] scat_k2_ftey[m][i];
        delete [] scat_k2_fthx[m][i];
        delete [] scat_k2_fthy[m][i];

        delete [] scat_j1_ftex[m][i];
        delete [] scat_j1_ftez[m][i];
        delete [] scat_j1_fthx[m][i];
        delete [] scat_j1_fthz[m][i];
        delete [] scat_j2_ftex[m][i];
        delete [] scat_j2_ftez[m][i];
        delete [] scat_j2_fthx[m][i];
        delete [] scat_j2_fthz[m][i];
      } 

      delete [] scat_k1_ftex[m];
      delete [] scat_k1_ftey[m];
      delete [] scat_k1_fthx[m];
      delete [] scat_k1_fthy[m];
      delete [] scat_k2_ftex[m];
      delete [] scat_k2_ftey[m];
      delete [] scat_k2_fthx[m];
      delete [] scat_k2_fthy[m];

      delete [] scat_j1_ftex[m];
      delete [] scat_j1_ftez[m];
      delete [] scat_j1_fthx[m];
      delete [] scat_j1_fthz[m];
      delete [] scat_j2_ftex[m];
      delete [] scat_j2_ftez[m];
      delete [] scat_j2_fthx[m];
      delete [] scat_j2_fthz[m];

      for(j = 0; j <= my_grid_ny; ++j)
      {
        delete [] scat_i1_ftez[m][j];
        delete [] scat_i1_ftey[m][j];
        delete [] scat_i1_fthz[m][j];
        delete [] scat_i1_fthy[m][j];
        delete [] scat_i2_ftez[m][j];
        delete [] scat_i2_ftey[m][j];
        delete [] scat_i2_fthz[m][j];
        delete [] scat_i2_fthy[m][j];
      }

      delete [] scat_i1_ftez[m];
      delete [] scat_i1_ftey[m];
      delete [] scat_i1_fthz[m];
      delete [] scat_i1_fthy[m];
      delete [] scat_i2_ftez[m];
      delete [] scat_i2_ftey[m];
      delete [] scat_i2_fthz[m];
      delete [] scat_i2_fthy[m];

    } // ++m


    delete [] scat_aux_ex;
    delete [] scat_aux_ey;
    delete [] scat_aux_hx;
    delete [] scat_aux_hy;

    delete [] scat_k1_ftex;
    delete [] scat_k1_ftey;
    delete [] scat_k1_fthx;
    delete [] scat_k1_fthy;
    delete [] scat_k2_ftex;
    delete [] scat_k2_ftey;
    delete [] scat_k2_fthx;
    delete [] scat_k2_fthy;

    delete [] scat_j1_ftex;
    delete [] scat_j1_ftez;
    delete [] scat_j1_fthx;
    delete [] scat_j1_fthz;
    delete [] scat_j2_ftex;
    delete [] scat_j2_ftez;
    delete [] scat_j2_fthx;
    delete [] scat_j2_fthz;

    delete [] scat_i1_ftez;
    delete [] scat_i1_ftey;
    delete [] scat_i1_fthz;
    delete [] scat_i1_fthy;
    delete [] scat_i2_ftez;
    delete [] scat_i2_ftey;
    delete [] scat_i2_fthz;
    delete [] scat_i2_fthy;

  } // end (if iscatt==1)

  //---------------------------------------------------------
  // NSOM
  //---------------------------------------------------------
  if(insom == 1)
  {
    delete [] transm_nsom_ptoi;
    delete [] transm_nsom_ptoj;
    delete [] transm_nsom_ex;
    delete [] transm_nsom_ey;
    delete [] transm_nsom_hx;
    delete [] transm_nsom_hy;
  }

  cout << "b." << endl;
  // CLOSE THE SIMULATION FILE
  if(isimfile==1)
  {
    simfile << " *Done cleaning up program." << endl << endl;
    simfile.close();
  }
  cout << "a." << endl;

  return;
}


//========================================================================
//========================================================================
//
//	NAME:	void read_parameters()
//	DESC:	Allows the user to set the simulation parameters
//
//	NOTES: 	i.
//
//
//========================================================================
//========================================================================
void read_parameters()
{

  int i;

  ifstream fin("parameters");

 
  // MPI
  get_comments(fin);
  fin >> mpi_nxprocs;
  get_comments(fin);
  fin >> mpi_nyprocs;
  get_comments(fin);
  fin >> mpi_nzprocs;

  get_comments(fin);
  fin >> isendtype;

  // PERIODICITY
  get_comments(fin);
  fin >> ixperiodic;
  get_comments(fin);
  fin >> iyperiodic;
  get_comments(fin);
  fin >> izperiodic;

  // COMPUTATIONAL GRID
  get_comments(fin);
  fin >> grid_xsize;
  get_comments(fin);
  fin >> grid_ysize;
  get_comments(fin);
  fin >> grid_zsize;

  get_comments(fin);
  fin >> grid_dx;
  get_comments(fin);
  fin >> grid_dy;
  get_comments(fin);
  fin >> grid_dz;

  // CPML
  get_comments(fin);
  fin >> cpml_layers;

  get_comments(fin);
  fin >> cpml_epsr;

  get_comments(fin);
  fin >> cpml_mur;

  get_comments(fin);
  fin >> cpml_kappamax;

  get_comments(fin);
  fin >> cpml_sigmamax_coeff;

  get_comments(fin);
  fin >> cpml_alphamax;

  get_comments(fin);
  fin >> cpml_m;

  get_comments(fin);
  fin >> cpml_ma;

  // TIME
  get_comments(fin);
  fin >> time_total;

  get_comments(fin);
  fin >> time_courant_factor;

  // TF/SF & SOURCE
  get_comments(fin);
  fin >> tfsf_zpos;

  get_comments(fin);
  fin >> source_intensity;

  get_comments(fin);
  fin >> src_ieypol;
  get_comments(fin);
  fin >> src_iexpol;
  get_comments(fin);
  fin >> src_i45pol;
  get_comments(fin);
  fin >> src_icircpol;

  get_comments(fin);
  fin >> source_wavelength;

  get_comments(fin);
  fin >> source_gauss_width;

  get_comments(fin);
  fin >> source_gauss_center;

  // TF/SF & SOURCE
  get_comments(fin);
  fin >> aux_epsr;

  get_comments(fin);
  fin >> aux_mur;

  // OUTPUT
  get_comments(fin);
  fin >> isimfile;

  get_comments(fin);
  fin >> noutput;

  get_comments(fin);
  fin >> output_nplanes;

  // !!! bad mem
  output_format = new int [output_nplanes+1];

  output_field = new int [output_nplanes+1];
  output_plane = new int [output_nplanes+1];
  output_plane_midpos = new double [output_nplanes+1];

  ft_wavelengths = new double [output_nplanes+1];

  for(i = 1; i <= output_nplanes; ++i)
  {
    get_comments(fin);
    fin >> output_format[i];

    get_comments(fin);
    fin >> output_field[i];

    get_comments(fin);
    fin >> ft_wavelengths[i];

    get_comments(fin);
    fin >> output_plane[i];

    get_comments(fin);
    fin >> output_plane_midpos[i];
  }

  // SCATTERING CROSS SECTION
  get_comments(fin);
  fin >> iscat_calc;

  get_comments(fin);
  fin >> scat_xcenter;

  get_comments(fin);
  fin >> scat_ycenter;

  get_comments(fin);
  fin >> scat_zcenter;

  get_comments(fin);
  fin >> scat_radius;

  // SPECTRA
  get_comments(fin);
  fin >> itransm;

  get_comments(fin);
  fin >> itransm_nf;

  get_comments(fin);
  fin >> spect_transm_pos;

  get_comments(fin);
  fin >> spect_refl_pos;

  get_comments(fin);
  fin >> nf_aperture_size;

  get_comments(fin);
  fin >> nf_aperture_xcenter;

  get_comments(fin);
  fin >> nf_aperture_ycenter;

  get_comments(fin);
  fin >> spect_npts;

  get_comments(fin);
  fin >> spect_minwave;

  get_comments(fin);
  fin >> spect_maxwave;

  // CLOSE THE FILE
  fin.close();


// *************************************
// *************************************
// *************************************
// *************************************
// *************************************


  //---------------------------------------------------------
  // NSOM
  //---------------------------------------------------------

  // FLAG TO MODEL NSOM
  insom = 0;

  // NSOM DIMENSIONS
  // i. the normal NSOM probe has an aperture width of 80nm
  // ii. the normal NSOM probe has a cylinder width of 240nm
  // iii. the normal NSOM probe has a tip height of 280nm
  // iv. the normal NSOM probe has a cylinder height of 320nm
  nsom_aperture_width = 80.0e-9;
  nsom_cylinder_width = 80.0e-9;
  nsom_tip_height = 0.0e-9;
  nsom_cylinder_height = 500.0e-9;

  // NSOM MATERIALS
  // i. the normal NSOM probe has an interior n of 1.5
  // ii. the normal NSOM probe has a cladidng width of 100nm
  // iii. the normal NSOM probe has a cladding eps of a PEC (infinity)
  nsom_interior_eps = 1.0*1.0;
  nsom_cladding_width = 0.0e-9;
  nsom_cladding_eps = 1.0*1.0;

  // NSOM POSITION (TIP)
  // i. the total grid size in z will be adjusted for this zcenter
  nsom_xcenter = grid_xsize/2.0;
  nsom_ycenter = grid_ysize/2.0;
  nsom_zcenter = 800.0e-9;

  // NSOM MEASURMENT WAVELENGTH
  // i. for 400nm arrays on n=1.5 substrate and n=1.0 superstrate the MAX of the 
  // SPP peak (for glass) is at 673nm in the zero order spectra (calculated)
  nsom_wavelength = 673.0e-9;


  return;
}



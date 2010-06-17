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
// NAME: JFDTD2D_CONFIG: CONFIGURATION FILE FOR JFDTD2D
//
// VERSION: Check www.thecomputationalphysicist.com for updated versions
//
// FILE: jfdtd2d_config.h
//
// AUTHOR: 	Jeffrey M. McMahon
//		jeffrey-mcmahon@northwestern.edu
//
//------------------------------------------------------------------------
//************************************************************************
//
// *last updated on 12/20/08
//
//------------------------------------------------------------------------
//************************************************************************
//
// NOTES:	i.
//
//
// TODO:	i. 
//
//************************************************************************
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//************************************************************************

//************************************************************************
//			INCLUDE FILES
//************************************************************************

#include "global.h"

//************************************************************************
//			CLASSES / STRUCTS
//************************************************************************

//************************************************************************
//			FULL CODE GLOBAL VARIABLES
//************************************************************************

int grid_nx, grid_ny, grid_nz;
int my_grid_nx, my_grid_ny, my_grid_nz, mpi_my_grid_nx1, mpi_my_grid_nx2, mpi_my_grid_ny1, mpi_my_grid_ny2, mpi_my_grid_nz1, mpi_my_grid_nz2;
double grid_xsize, grid_ysize, grid_zsize, grid_dx, grid_dy, grid_dz;

int ***grid_material_ex, ***grid_material_ey, ***grid_material_ez;
double ***permittivity_ex, ***permittivity_ey, ***permittivity_ez, ***permeability_hx, ***permeability_hy, ***permeability_hz;

// !!! temp
int slayer_material[2];
double slayer_eps[2], slayer_mur[2];

// NSOM PROBE
// !!! experimental
int insom;
int nsom_ift;

// DIMENSIONS
double nsom_aperture_width, nsom_cylinder_width, nsom_tip_height, nsom_cylinder_height, nsom_cladding_width, nsom_interior_eps, nsom_cladding_eps;
double nsom_xcenter, nsom_ycenter, nsom_zcenter;
// TRANSMISSION CHARACTERISTICS
int transm_nsom_npts, transm_nsom_kpt, *transm_nsom_ptoi, *transm_nsom_ptoj, transm_nsom_ptok;
double transm_nsom_zpos, nsom_wavelength, nsom_wavenumber;


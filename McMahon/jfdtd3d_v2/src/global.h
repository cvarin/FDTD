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
// FILE: global.h
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

#ifndef GLOBAL_H	
#define GLOBAL_H


//************************************************************************
//			INCLUDE FILES
//************************************************************************

//************************************************************************
//			CLASSES / STRUCTS
//************************************************************************

//************************************************************************
//			FULL CODE GLOBAL VARIABLES
//************************************************************************

extern int grid_nx, grid_ny, grid_nz;
extern int my_grid_nx, my_grid_ny, my_grid_nz, mpi_my_grid_nx1, mpi_my_grid_nx2, mpi_my_grid_ny1, mpi_my_grid_ny2, mpi_my_grid_nz1, mpi_my_grid_nz2;
extern double grid_xsize, grid_ysize, grid_zsize, grid_dx, grid_dy, grid_dz;

extern int ***grid_material_ex, ***grid_material_ey, ***grid_material_ez;
extern double ***permittivity_ex, ***permittivity_ey, ***permittivity_ez, ***permeability_hx, ***permeability_hy, ***permeability_hz;

// !!! temp
extern int slayer_material[2];
extern double slayer_eps[2], slayer_mur[2];

// NSOM PROBE
// !!! experimental
extern int insom;
extern int nsom_ift;

// DIMENSIONS
extern double nsom_aperture_width, nsom_cylinder_width, nsom_tip_height, nsom_cylinder_height, nsom_cladding_width, nsom_interior_eps, nsom_cladding_eps;
extern double nsom_xcenter, nsom_ycenter, nsom_zcenter;
// TRANSMISSION CHARACTERISTICS
extern int transm_nsom_npts, transm_nsom_kpt, *transm_nsom_ptoi, *transm_nsom_ptoj, transm_nsom_ptok;
extern double transm_nsom_zpos, nsom_wavelength, nsom_wavenumber;

#endif	

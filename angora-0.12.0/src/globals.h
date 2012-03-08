/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GLOBALS_H
#define GLOBALS_H

#include "headers.h"

#include "material_id.h"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//global grid variables

/*********************/
/*	ANGORA VERSION **/
/*********************/
//Full version number: MAJOR.MINOR.REVISION
extern int angora_version_major;		//major version number of the Angora package
extern int angora_version_minor;		//minor version number of the Angora package
extern int angora_version_revision;		//revision number of the Angora package


/*******************************/
/*	FILE AND DIRECTORY NAMES ***/
/*******************************/
//current working directory
extern string currentworkdir;
//base path for all other paths  (defaults to working directory)
extern string angora_basepath;
//output path
extern string OutputDir;
//input path
extern string InputDir;
//log output path
extern string LogOutputDir;
//recorder output path
extern string RecorderOutputDir;
//movie-recorder output path
extern string MovieRecorderOutputDir;
//line-recorder output path
extern string LineRecorderOutputDir;
//field-value-recorder output path
extern string FieldValueRecorderOutputDir;
//time-domain near-field-to-far-field transform output path
extern string TimeDomainNFFFTOutputDir;
//phasor-domain near-field-to-far-field transform output path
extern string PhasorDomainNFFFTOutputDir;


/***************************************/
/*	DEFAULT FILE AND DIRECTORY NAMES ***/
/***************************************/
//default output path
extern const string default_OutputDir;
//default input path
extern const string default_InputDir;
//default log output path
extern const string default_LogOutputDir;
//default log file name
extern const string default_LogFileName;
//default recorder output path
extern const string default_RecorderOutputDir;
//default movie-recorder output path
extern const string default_MovieRecorderOutputDir;
//default line-recorder output path
extern const string default_LineRecorderOutputDir;
//default field-value-recorder output path
extern const string default_FieldValueRecorderOutputDir;
//default movie filename
extern const string default_movie_filename;
//default movie file extension
extern const string default_movie_fileextension;
//default line filename
extern const string default_line_filename;
//default line file extension
extern const string default_line_fileextension;
//default field-value filename
extern const string default_fieldvalue_filename;
//default field-value file extension
extern const string default_fieldvalue_fileextension;
//default time-domain near-field-to-far-field transform output path
extern const string default_TimeDomainNFFFTOutputDir;
//default phasor-domain near-field-to-far-field transform output path
extern const string default_PhasorDomainNFFFTOutputDir;
//default time-domain output filename
extern const string default_TimeDomainNFFFT_filename;
//default phasor-domain output filename
extern const string default_PhasorDomainNFFFT_filename;
//default optical image output path
extern const string default_OpticalImagingOutputDir;


/********************************/
/*	GRID DIMENSION VARIABLES	*/
/********************************/
extern double courant,dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
extern int NSTEPS;
//	Indices of the origin cell
extern int OriginX,OriginY,OriginZ;

/********************/
/*	FIELD ARRAYS	*/
/********************/
//	E-field
extern Array<double,3> Ex,Ey,Ez;
//	H-field
extern Array<double,3> Hx,Hy,Hz;


/*********************/
/*	MATERIAL ARRAYS	 */
/*********************/
//Material Pointer Arrays
//	Material pointer arrays for the E-field
extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
//	Material pointer arrays for the H-field
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;


/****************************************/
/*	LAYERING AND MATERIAL INFORMATION	*/
/****************************************/
//Material type information
//number of distinct (eps_x,cond_e_x) combinations present in the simulation space (must be less than sizeof(ElectricMaterialIndexType))
extern int NumOfElectricMaterials_X;
//number of distinct (eps_y,cond_e_y) combinations present in the simulation space (must be less than sizeof(ElectricMaterialIndexType))
extern int NumOfElectricMaterials_Y;
//number of distinct (eps_z,cond_e_z) combinations present in the simulation space (must be less than sizeof(ElectricMaterialIndexType))
extern int NumOfElectricMaterials_Z;
//number of distinct (mu_x,cond_h_x) combinations present in the simulation space (must be less than sizeof(MagneticMaterialIndexType))
extern int NumOfMagneticMaterials_X;
//number of distinct (mu_y,cond_h_y) combinations present in the simulation space (must be less than sizeof(MagneticMaterialIndexType))
extern int NumOfMagneticMaterials_Y;
//number of distinct (mu_z,cond_h_z) combinations present in the simulation space (must be less than sizeof(MagneticMaterialIndexType))
extern int NumOfMagneticMaterials_Z;
extern const int PEC;
extern const int vacuum;
//Layering information
//According to electric properties:
extern Array<ElectricMaterialIndexType_X,1> Layering_e_x;
extern Array<ElectricMaterialIndexType_Y,1> Layering_e_y;
extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
//According to magnetic properties:
extern Array<MagneticMaterialIndexType_X,1> Layering_h_x;
extern Array<MagneticMaterialIndexType_Y,1> Layering_h_y;
extern Array<MagneticMaterialIndexType_Z,1> Layering_h_z;
// Layers are counted from the bottom up, i.e., the lowermost layer is layer #0, and the uppermost layer is layer #N-1 (where N is the number of layers)
//Number of different layers in the grid
extern int number_of_layers;
//Array of material indices for each layer
//electric properties:
extern Array<ElectricMaterialIndexType_X,1> LayerElectricMaterial_X;
extern Array<ElectricMaterialIndexType_Y,1> LayerElectricMaterial_Y;
extern Array<ElectricMaterialIndexType_Z,1> LayerElectricMaterial_Z;
//magnetic properties:
extern Array<MagneticMaterialIndexType_X,1> LayerMagneticMaterial_X;
extern Array<MagneticMaterialIndexType_Y,1> LayerMagneticMaterial_Y;
extern Array<MagneticMaterialIndexType_Z,1> LayerMagneticMaterial_Z;
//Array of cell indices that define the lower boundary of each layer (index of the lowest cell that is included in each layer)
extern Array<int,1> LayerLowerZIndices;
//thicknesses (in grid cells) of each layer
extern Array<int,1> LayerThicknesses;
//boolean array denoting whether there is a perfect-electric-conductor (PEC) sheet below each layer
extern Array<bool,1> IsLayerGrounded;
//Update coefficient arrays
//	Update coefficients for the E-field
extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
//	Update coefficients for the H-field
extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;
//Constitutive parameter arrays
extern Array<double,1> eps_x,eps_y,eps_z;
extern Array<double,1> mu_x,mu_y,mu_z;
extern Array<double,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<double,1> cond_h_x,cond_h_y,cond_h_z;
//Cordinate stretching (kappa) arrays
extern Array<double,1> kappa_e_x,kappa_e_y,kappa_e_z,kappa_h_x,kappa_h_y,kappa_h_z;

//maximum and minimum constitutive parameters in the grid
extern double epsilon_r_max_x,epsilon_r_min_x,epsilon_r_max_y,epsilon_r_min_y,epsilon_r_max_z,epsilon_r_min_z;
extern double mu_r_max_x,mu_r_min_x,mu_r_max_y,mu_r_min_y,mu_r_max_z,mu_r_min_z;
extern double epsilon_r_max,epsilon_r_min;
extern double mu_r_max,mu_r_min;

extern double epsilon_r_upper,mu_r_upper; 	//relative permittivity and permeability of the uppermost layer (changed only in load_geometry())
extern double epsilon_r_lower,mu_r_lower; 	//relative permittivity and permeability of the lowermost layer (changed only in load_geometry())
extern double c_upper;	//velocity of propagation in the uppermost layer (changed only through PlaceSlab in geometry.cpp)
extern double c_lower;	//velocity of propagation in the lowermost layer (changed only through PlaceSlab in geometry.cpp)


/******************************/
/*	MULTIPLE-GRID VARIABLES   */
/******************************/
// int number_of_runs;	//Number of independent grids (i.e. subgroups) in the communicator
extern const int default_number_of_runs;	//Number of FDTD simulation runs is 1 by default
extern int number_of_runs;	//Number of FDTD simulation runs
// int GridIndex;	//index of the grid (i.e. subgroup) that the node belongs to (=0...number_of_runs-1)
extern int GridIndex;	//index of the current simulation run
//with the abandonment of the idea of multiple grids in one communicator, the slowest grid idea is also abandoned
// int SlowestGrid = 0;	//tentatively, the index of the slowest grid (by default, 0)
extern Array<bool,1> grid_is_enabled;	//indices of the grids that are enabled (default: all grids from 0 to NumOfGrids-1)

/********************/
/*	MPI VARIABLES   */
/********************/
#ifndef MPI_DISABLE
//MPI Communicators
extern MPI_Comm MPI_SubComm;		//subcommunicator for the grid
extern MPI_Comm MPI_CartSubComm;	//cartesian subcommunicator for the grid
extern MPI_Status Status;
#endif
//Basic MPI variables
extern int rank, nodes;
//Ranks of adjacent nodes
extern int rank_behind,rank_front,rank_left,rank_right,rank_below,rank_above;
//Grid limits determined by node_sectioning() and node_limits()
extern int rank_x, rank_y, rank_z;		//section rank (or coordinate) in the x, y, z directions
extern int nodes_x, nodes_y, nodes_z;	//number of sections in the x, y, z directions
extern int iback,ifront;	//x-indices of rearmost and foremost cells in the grid
					// that is maintained in this node
extern int jleft,jright;	//y-indices of leftmost and rightmost cells in the grid
					// that is maintained in this node
extern int klower,kupper;	//z-indices of lowermost and uppermost cells in the grid
					// that is maintained in this node
extern int FullIntPosMin_x,FullIntPosMax_x;	//min. and max. indices of field components at full-integer grid positions that are updated at the node
extern int FullIntPosMin_y,FullIntPosMax_y;
extern int FullIntPosMin_z,FullIntPosMax_z;

extern int gridwidth_x, gridwidth_y, gridwidth_z;		//width of the x, y, z sections
//MPI send and receive buffer arrays:
extern Array<double,2> SendBuf_Hx_y,SendBuf_Hz_y,RecvBuf_Hx_y,RecvBuf_Hz_y,
 SendBuf_Hy_x,SendBuf_Hz_x,RecvBuf_Hy_x,RecvBuf_Hz_x,
 SendBuf_Hx_z,SendBuf_Hy_z,RecvBuf_Hx_z,RecvBuf_Hy_z;


/********************************************************/
/*	MAXIMUM FIELD VALUE IN GRID, AND USEFUL ACCURACY	*/
/********************************************************/
extern double max_field_value;	//maximum electric field value encountered in grid
extern double accuracy;		//useful accuracy range of the grid (in -dB)

/********************************/
/*	GLOBAL ITERATION INDICES	*/
/********************************/
extern int i,j,k;

#endif

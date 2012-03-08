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

#ifndef CPW_FS_H
#define CPW_FS_H

//Declaration of the class "Cpw_fs" for a TF/SF plane-wave source in free space

#include "Cpw.h"

extern double epsilon_r_upper;


class Cpw_fs : public Cpw
{
 public:
	 Cpw_fs(const PWDataType& MyData, const double& epsilon_r_space = epsilon_r_upper);	//constructor

	 void WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Cpw_fs object
	 void WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Cpw_fs object
	 void WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Cpw_fs object

 private:
 	 //Given the coordinates of the field component on the TF/SF boundary, the functions below calculate the incident field by projection onto the 1-D grid
	 inline double IncidentHx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHz(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEz(const int& i, const int& j, const int& k, const int& n);

 protected:
 	 //Auxiliary grid updates
	 void UpdateIncidentE(const int& n);
	 void UpdateIncidentH(const int& n);
	 //These methods apply field corrections in the main grid (by projection onto the auxiliary grid or any other method)
	 void ApplyCorrectionE(const int& n);
	 void ApplyCorrectionH(const int& n);

	 void UpdateE(const int& n);	//update incident E field at time index n
	 void UpdateH(const int& n);	//update incident H field at time index n

//	 //incident field waveform
//	 void CreateIncidentWaveform();
	 //E-field waveform of the incident homogeneous plane wave in the uppermost layer at the n'th time step
	 double waveformE_value(const double& n);

 	 //velocity of light in the homogeneous space
	 double c_space;
 	 //wave impedance in the homogeneous space
	 double Z_space;
 	 //grid velocity in the homogeneous space
	 double vp_space;

	 //grid spacing for the 1-D auxiliary grid (different from the 3-D grid)
	 //found using the matched numerical dispersion (MND) method (approximate version - 5.69 in Taflove 2005)
	 double dx_g;
	 double dx_over_dx_g;	//ratio of dx to dx_g

	 TinyVector<double,3> k_E,k_H;	//unit vectors that point in the direction of incident E and H
	 double k_Ex,k_Ey,k_Ez,k_Hx,k_Hy,k_Hz;	//components of the above unit vectors

	 //extra cells at the top and bottom of the 1-D grid (since H-field positions are 1/2 cell outside the TF/SF box)
	 int extra_top,extra_bottom;	//in 1-D grid cells (dx_g)
	 TinyVector<double,3> extra_distance;	//vectorial position of the final E-field position on the 1-D grid,
											//measured from the corner of the TF/SF box in 3-D grid cells (dx)
											//note the conversion factor dx_over_dx_g to convert from 1-D grid cells to 3-D grid cells

	 int AuxGridLength;		//length of the 1-D grid (in grid cells)
	 int AuxGridLower,AuxGridUpper;	//limits of the 1-D grid (indices of the outermost grid cells)

	 //1-D grids for incident field calculation
	 Array<double,1> AuxGrid_E;		//E-field
	 Array<double,1> AuxGrid_H;		//H-field
	 Array<double,1> AuxGrid_Ex,AuxGrid_Ey,AuxGrid_Ez;		//components of the E-field
	 Array<double,1> AuxGrid_Hx,AuxGrid_Hy,AuxGrid_Hz;		//components of the H-field


	 //necessary number of time steps to cover the incident waveform
	 int INCSTEPS;

	 //update coefficients for E and H updates
	 double Ca,Cb,Da,Db;

//	 //E-field waveform of the incident homogeneous plane wave in the uppermost layer
//	 Array<double,1> waveformE;

	 double PositionOnAuxGrid,Einc_interp,Hinc_interp;	//intermediate values during incident field calculation
	 int IntegerPositionOnAuxGrid;	//integer position on the 1-D grid = int(PositionOnAuxGrid). Always use this instead of using int(PositionOnAuxGrid) everywhere, because if PositionOnAuxGrid is infinitesimally smaller than an integer, int(PositionOnAuxGrid) might yield DIFFERENT results depending on where it is used (for example, if int(PositionOnAuxGrid) is passed as an index parameter to a Blitz++ array, the result is sometimes 1 larger than the correct answer - this has been observed in Linux builds, not in VS2005)

	 int HardSourcePoint,AbsorbPointLower;	//hard-source point and lower absorbing point for the plane wave in 1-D grid
	 TinyVector<double,3> Coordinate;	//coordinate of incident field (measured from the plane-wave origin)
	 //aux-grid-coordinate contributions at different x,y,and z positions
	 //these are pre-computed arrays for the efficient projection on the 1-D auxiliary grid
	 Array<double,1> CoordinateOnAuxGrid_x_fullint,CoordinateOnAuxGrid_y_fullint,CoordinateOnAuxGrid_z_fullint,
	 				CoordinateOnAuxGrid_x_halfint,CoordinateOnAuxGrid_y_halfint,CoordinateOnAuxGrid_z_halfint;

	 Array<double,2> PreviousELower;	//lower storage values for the 2nd-order ABC used in the 1-D grid (1st dimension: upward on the auxiliary grid, 2nd dimension: backward in time)
	 double Sc_1D_upper,Sc_1D_lower;	//effective courant numbers for the 1-D grid

	 double a;	//interpolation parameter
};

#endif

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

#ifndef CPW_ML_H
#define CPW_ML_H

//Declaration of the class "Cpw_ml" for a TF/SF plane-wave source in a multilayered medium

#include "Cpw.h"


class Cpw_ml : public Cpw
{
 public:
	 Cpw_ml(const PWDataType& MyData);	//constructor

	 void WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Cpw_ml object
	 void WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Cpw_ml object
	 void WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Cpw_ml object

 private:
	 //Auxiliary grid updates
	 void UpdateIncidentE(const int& n);
	 void UpdateIncidentH(const int& n);
	 //These methods apply field corrections in the main grid (by projection onto the auxiliary grid or any other method)
	 void ApplyCorrectionE(const int& n);
	 void ApplyCorrectionH(const int& n);

	 //Given the coordinates of the field component on the TF/SF boundary, the functions below calculate the incident field by projection onto the time history of the 1-D grid
	 inline double IncidentHx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentHz(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEx(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEy(const int& i, const int& j, const int& k, const int& n);
	 inline double IncidentEz(const int& i, const int& j, const int& k, const int& n);

	 //routine to automatically determine the evanescent regions
	 void DetermineEvanescentRegions();
	 void FractionalUpdate(const int& n);		//update incident fields by a fraction of dt (=dt_g)
	 void UpdateE(const int& n);	//update incident E field at time index n
	 void UpdateH(const int& n);	//update incident H field at time index n
	 void UpdateEvanescentE();	//update E-field in regions with inhomogeneous waves
	 void UpdateEvanescentH();	//update H-field in regions with inhomogeneous waves

	 //calculates the auxiliary-grid time step
	 void calculate_auxgrid_time_step();

	 //TM and TE components of the incident E and H waveforms at the m'th auxiliary-grid time step. The H-values are located 1/2 space cells below the E-values.
	 complex<double> Incident_Ee(const double& m);
	 complex<double> Incident_Eh(const double& m);
	 complex<double> Incident_He(const double& m);
	 complex<double> Incident_Hh(const double& m);

	 //boost shared pointer Hilbert transform of the incident E-field waveform
	 Cwf_shared_ptr waveformHilbert;

	 //basic waveform shape
	 double waveform(double t, double tau);

 	 //grid velocity at the uppermost layer
	 double vp_upper;

	 //time expansion factor (to avoid instability):
	 int timefactor;	// dt will be divided by this amount in the 1-D grid updates.
						// This is necessary because the phase velocity on the grid is greater than c. In order to keep the Courant factor the same, either the grid spacing has to be made larger (which is unacceptable due to increased interpolation error), or the time step has to be reduced.
	 double dt_g;	// time factor used in 1-D grid : dt_g = dt/timefactor

	 TinyVector<double,3> k_Ee,k_Eh,k_He,k_Hh;	//unit vectors that point in the direction of Ee,Eh,He,Hh

	 //1-D grids for incident wave calculation (closest grid to initial point of contact)
	 int AuxGridLower,AuxGridUpper;	//limits of the 1-D grid

	 //1-D auxiliary grids for incident field calculation
	 Array<complex<double>,1> AuxGrid_Ee;	//TM component of E
	 Array<complex<double>,1> AuxGrid_Eh;	//TE component of E
	 Array<complex<double>,1> AuxGrid_Ez;	//z-component of E

	 Array<complex<double>,1> AuxGrid_He;	//TM component of H
	 Array<complex<double>,1> AuxGrid_Hh;	//TE component of H
	 Array<complex<double>,1> AuxGrid_Hz;	//z-component of H
	 Array<complex<double>,1> AuxGrid_He_pr;	//Auxiliary variable He'

	 //time histories of field components on the 1-D grid (composed of time instants above)
	 Array<double,2> TimeHistory_Ee;
	 Array<double,2> TimeHistory_Eh;
	 Array<double,2> TimeHistory_Ez;

	 Array<double,2> TimeHistory_He;
	 Array<double,2> TimeHistory_Hh;
	 Array<double,2> TimeHistory_Hz;

	 int timeindex;	//holds the position of the 1-D grid update within the main update cycle
					// 0<=timeindex<timefactor

	 //marker arrays for indicating the regions where the wave is inhomogeneous
	 Array<bool,1> IsEvanescent_xy,IsEvanescent_z;	//is the incident wave evanescent at this z-position? (first array for integer space locations, second for half-integer space locations)
	 bool EvanescentWavePresent;	//is there a region in the grid where the incident field is evanescent?

	 //structure representing an evanescent layer
	 struct SingleEvanescentLayer
	 {
		 int lowerE, upperE, lowerH, upperH, N_e, N_h, nrhs;
		 Array <double,1> diagE_r, diagE_i, lowdiagE_r, lowdiagE_i, updiagE_r, updiagE_i, righthandE_r, righthandE_i;
		 Array <double,1> diagH_r, diagH_i, lowdiagH_r, lowdiagH_i, updiagH_r, updiagH_i, righthandH_r, righthandH_i;
		 bool homogeneous_at_interfaces;	//is the plane wave homogeneous at the interfaces?
	 };
	 int info;	//temporary status variable used in zgtsv_
	 //array that holds all evanescent layers in the layering structure
	 Array<SingleEvanescentLayer,1> EvanescentLayers;
	 int NumOfEvanescentLayers;	//number of evanescent blocks

	 //relative permittivity (real part) arrays (for E and H updates in 1-D grid)
	 Array<double,1> epsilon_E,epsilon_H;
	 //conductivity arrays (for E and H updates in 1-D grid)
	 Array<double,1> cond_E,cond_H;

	 //update coefficients for E and H updates
	 Array<complex<double>,1> Ca_e,Cb_e,Ca_h,Cb_h;
	 Array<complex<double>,1> Da_e,Db_e,Da_e_pr,Db_e_pr,Da_h,Db_h;
	 //update coefficients for converting from auxiliary varible He' to He
	 Array<complex<double>,1> A_pr,B_pr;

	 //conversions coefficients for Ez and Hz
	 Array<complex<double>,1> Ez_e_pr,Ez_e;	//Ez is obtained from either He' or directly from He
	 Array<complex<double>,1> Hz_h;	//Hz is obtained from Eh

//	 Array<complex<double>,1> waveformE_analytic,waveformH_analytic;	//analytic waveforms of the incident homogeneous plane wave in the uppermost layer

	 double PositionOnPrincipal,TimeDelay,DelayedTime,Einc_interp,Hinc_interp;	//intermediate values during incident field calculation
	 int DelayedIntegerTime;	//integer time position in the time history = int(DelayedTime). Always use this instead of using int(DelayedTime) everywhere, because if DelayedTime is infinitesimally smaller than an integer, int(DelayedTime) might yield DIFFERENT results depending on where it is used (for example, if int(DelayedTime) is passed as an index parameter to a Blitz++ array, the result is sometimes 1 larger than the correct answer - this has been observed in Linux builds, not in VS2005)

	 double PositionInTime,Einc_e,Einc_h,Einc_z,Hinc_e,Hinc_h,Hinc_z;
	 int FeedPoint,AbsorbPointLower,AbsorbPointUpper;	//feed point and lower-upper absorbing points for the plane wave in 1-D grid
	 TinyVector<double,3> Coordinate;	//coordinate of incident field w.r.t. to the plane-wave origin
	 Array<complex<double>,2> PreviousEeLower,PreviousEeUpper,PreviousEhLower,PreviousEhUpper;	//lower and upper storage values for the 2nd-order ABC used in the 1-D grid (1st dimension: upward or downward the auxiliary grid, 2nd dimension: backward in time)
	 double Sc_1D_upper,Sc_1D_lower;	//effective courant numbers for the 1-D grid

	 complex<double> fieldstore;	//storage variable for storing He' values while converting to He

	 double a;	//interpolation parameter
};

#endif

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

//Declaration of the class "Cpw_fs" for a TF/SF plane-wave source in free space

#include "headers.h"

#include "Cpw_fs.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//Uses TinyVector operations
#include <blitz/tinyvec-et.h>

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern double dx,dt,courant;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML,NSTEPS;

extern int number_of_layers;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;

extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;


Cpw_fs::Cpw_fs(const PWDataType& MyData, const double& epsilon_r_space)
		: Cpw(MyData)	//calls the base class constructor first
{
	//Don't do this, because the 2-layer injector consists of 3 free-space injectors
// 	if (number_of_layers!=1)
// 	{//if the # of layers is not 1, there is a developing error, so display error message and exit (better error-reporting scheme later!?)
// 		if (rank==0)
// 		{
// 			cout << "Error: Cpw_fs class cannot be used if number of layers is not 1." << endl;
// 		}
// 		MPI_Barrier(MPI_CartSubComm);
// 		exit(-1);
// 	}

	c_space = c/sqrt(epsilon_r_space); /** Assuming mu_r=1 **/
	Z_space = sqrt(1/epsilon_r_space)*eta_0; /** Assuming mu_r=1 **/

	//Unit vector that points in the direction of the E-field vector
	k_E = CosPsi*cross(k_inc_lateral,unit_z)+SinPsi*cross(cross(k_inc_lateral,unit_z),k_inc);
	k_Ex = k_E(firstDim);
	k_Ey = k_E(secondDim);
	k_Ez = k_E(thirdDim);
	//Unit vector that points in the direction of the H-field vector
	k_H = cross(k_inc,k_E);
	k_Hx = k_H(firstDim);
	k_Hy = k_H(secondDim);
	k_Hz = k_H(thirdDim);

	//find the grid spacing for the 1-D auxiliary grid (approximate MND method - 5.69 in Taflove 2005)
	dx_over_dx_g = 1/sqrt(pow(sin(THETA_INC),4)*(pow(cos(PHI_INC),4)+pow(sin(PHI_INC),4))+pow(cos(THETA_INC),4));

	dx_g = dx/dx_over_dx_g;

	//Length of the 1-D grid (in grid cells)
	//adjusted for worst case: incident direction coincident with diagonal
	AuxGridLength = (int)ceil(sqrt(pow((double)(PWFrontX-PWBackX+1),2)+pow((double)(PWRightY-PWLeftY+1),2)+pow((double)(PWUpperZ-PWLowerZ+1),2))
								*dx_over_dx_g);	//convert to 1-D auxiliary-grid cells
	//extra cells at the top and bottom of the 1-D grid (since H-field positions are 1/2 cell outside the TF/SF box)
	extra_top = 2;	//in 1-D grid cells (dx_g)
	extra_bottom = 2;	//in 1-D grid cells (dx_g)
	extra_distance = -extra_top*k_inc/dx_over_dx_g;

	//Lower and upper limits of the 1-D grid are fixed
	AuxGridLower = 1;
	AuxGridUpper = AuxGridLength+extra_top+extra_bottom;
	//hard-source and absorption points
	AbsorbPointLower = AuxGridLower;	//E-field absorption point (lower)
	HardSourcePoint = AuxGridUpper+1;		//E-field hard-source point is right at uppermost E-field position

	//the coordinates of the plane-wave origin
	//(measured from the back-left-bottom point of the grid)
	PWOrigin = (Data.PWOriginX-1),(Data.PWOriginY-1),(Data.PWOriginZ-1);

	//Determine the coordinates of the plane-wave origination point
	//(measured from the back-left-bottom point of the grid)
	//if incident from upper hemisphere:
	if ((k_inc(firstDim)<=0)&&(k_inc(secondDim)<=0)&&(k_inc(thirdDim)<=0))
	{
		PWOriginationPoint = (PWFrontX),(PWRightY),(PWUpperZ);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)>0)&&(k_inc(secondDim)<=0)&&(k_inc(thirdDim)<=0))
	{
		PWOriginationPoint = (PWBackX-1),(PWRightY),(PWUpperZ);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)>0)&&(k_inc(secondDim)>0)&&(k_inc(thirdDim)<=0))
	{
		PWOriginationPoint = (PWBackX-1),(PWLeftY-1),(PWUpperZ);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)<=0)&&(k_inc(secondDim)>0)&&(k_inc(thirdDim)<=0))
	{
		PWOriginationPoint = (PWFrontX),(PWLeftY-1),(PWUpperZ);
		PWOriginationPoint += extra_distance;
	}
	//if incident from lower hemisphere:
	if ((k_inc(firstDim)<=0)&&(k_inc(secondDim)<=0)&&(k_inc(thirdDim)>0))
	{
		PWOriginationPoint = (PWFrontX),(PWRightY),(PWLowerZ-1);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)>0)&&(k_inc(secondDim)<=0)&&(k_inc(thirdDim)>0))
	{
		PWOriginationPoint = (PWBackX-1),(PWRightY),(PWLowerZ-1);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)>0)&&(k_inc(secondDim)>0)&&(k_inc(thirdDim)>0))
	{
		PWOriginationPoint = (PWBackX-1),(PWLeftY-1),(PWLowerZ-1);
		PWOriginationPoint += extra_distance;
	}
	if ((k_inc(firstDim)<=0)&&(k_inc(secondDim)>0)&&(k_inc(thirdDim)>0))
	{
		PWOriginationPoint = (PWFrontX),(PWLeftY-1),(PWLowerZ-1);
		PWOriginationPoint += extra_distance;
	}

	//components of PWOriginationPoint
	PWOriginationPoint_x = PWOriginationPoint(firstDim);
	PWOriginationPoint_y = PWOriginationPoint(secondDim);
	PWOriginationPoint_z = PWOriginationPoint(thirdDim);

	//aux-grid-coordinate contributions at different x,y,and z positions
	CoordinateOnAuxGrid_x_fullint.resize(Range(1,NCELLS_X+2*NPML+1));
	CoordinateOnAuxGrid_y_fullint.resize(Range(1,NCELLS_Y+2*NPML+1));
	CoordinateOnAuxGrid_z_fullint.resize(Range(1,NCELLS_Z+2*NPML+1));
	CoordinateOnAuxGrid_x_halfint.resize(Range(1,NCELLS_X+2*NPML));
	CoordinateOnAuxGrid_y_halfint.resize(Range(1,NCELLS_Y+2*NPML));
	CoordinateOnAuxGrid_z_halfint.resize(Range(1,NCELLS_Z+2*NPML));

	//aux-grid-coordinate contributions at different x,y,and z positions
	for (int i=1; i<=NCELLS_X+2*NPML+1; i++)
	{
		CoordinateOnAuxGrid_x_fullint(i) = (i-1-PWOriginationPoint_x)*k_inc_x*dx_over_dx_g;
	}
	for (int j=1; j<=NCELLS_Y+2*NPML+1; j++)
	{
		CoordinateOnAuxGrid_y_fullint(j) = (j-1-PWOriginationPoint_y)*k_inc_y*dx_over_dx_g;
	}
	for (int k=1; k<=NCELLS_Z+2*NPML+1; k++)
	{
		CoordinateOnAuxGrid_z_fullint(k) = (k-1-PWOriginationPoint_z)*k_inc_z*dx_over_dx_g;
	}
	for (int i=1; i<=NCELLS_X+2*NPML; i++)
	{
		CoordinateOnAuxGrid_x_halfint(i) = (i-0.5-PWOriginationPoint_x)*k_inc_x*dx_over_dx_g;
	}
	for (int j=1; j<=NCELLS_Y+2*NPML; j++)
	{
		CoordinateOnAuxGrid_y_halfint(j) = (j-0.5-PWOriginationPoint_y)*k_inc_y*dx_over_dx_g;
	}
	for (int k=1; k<=NCELLS_Z+2*NPML; k++)
	{
		CoordinateOnAuxGrid_z_halfint(k) = (k-0.5-PWOriginationPoint_z)*k_inc_z*dx_over_dx_g;
	}

	//calculate grid velocity in the homogeneous space
	vp_space = c_space*GridVelocity(THETA_INC,PHI_INC);
// 	vp_space = c_space;

	//time needed for the plane wave to reach the origin from the plane-wave origination point
	OriginDelay = dot(-k_inc,PWOriginationPoint-PWOrigin)*dx/vp_space;
	//note the use of dx instead of dx_g, because PWOriginationPoint and PWOrigin are measured in dx
	//GridVelocity(THETA_INC,PHI_INC) used because it can cause differences (on the order of dt) for big grids

	//if the there is a sudden jump in the excitation waveform at the plane-wave origin, set the initial time value back
	if ((Data.waveform->starting_time()-OriginDelay)<get_initial_time_value())
	{
		set_initial_time_value(Data.waveform->starting_time()-OriginDelay);
	}

	//accuracy is less than -50 dB with the MND method
	if (!dB_accuracy_set_in_configfile) dB_accuracy = -50;

//	//create the incident waveform
//	CreateIncidentWaveform();

// 	//calculate the velocity anisotropy factor
// 	double vel_anisotropy = GridVelocity(0,0)/GridVelocity(THETA_INC,PHI_INC);
//
// 	//coefficients for the E update
// 	Ca = 1.0;
// 	Cb = dt/dx/(epsilon_r_space*epsilon_0*vel_anisotropy);	//apply the velocity anisotropy factor
//
// 	//coefficients for the H update
// 	Da = 1.0;
// 	Db = dt/dx/(mu_0*vel_anisotropy);	//apply the velocity anisotropy factor

	//coefficients for the E update
	Ca = 1.0;
//	Cb = dt/dx_g/(epsilon_r_space*epsilon_0);	//note the use of dx_g instead of dx
	Cb = (dt/dx_g)*c_space*Z_space;

	//coefficients for the H update
	Da = 1.0;
//	Db = dt/dx_g/mu_0;		//note the use of dx_g instead of dx
	Db = (dt/dx_g)*c_space/Z_space;		//note the use of dx_g instead of dx

	//Auxiliary 1-D grids for incident field computation
	AuxGrid_E.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_H.resize(Range(AuxGridLower,AuxGridUpper));
	AuxGrid_E = 0;
	AuxGrid_H = 0;
	//x,y,z components
	AuxGrid_Ex.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_Ey.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_Ez.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_Hx.resize(Range(AuxGridLower,AuxGridUpper));
	AuxGrid_Hy.resize(Range(AuxGridLower,AuxGridUpper));
	AuxGrid_Hz.resize(Range(AuxGridLower,AuxGridUpper));
	AuxGrid_Ex = 0;
	AuxGrid_Ey = 0;
	AuxGrid_Ez = 0;
	AuxGrid_Hx = 0;
	AuxGrid_Hy = 0;
	AuxGrid_Hz = 0;


	//effective courant number for 1D grid:
	Sc_1D_lower = c_space*dt/dx_g;		//note the use of dx_g instead of dx

	//resize and initialize ABC storage variables PreviousEeLower(k,n) etc.
	//NOTE: k=0,n=0 need not be used, since this is the current value of the E-field at the boundary of the auxiliary grid.
	//However, it will still be used for symmetry (see Schneider course notes on 2nd order ABC).
	PreviousELower.resize(Range(0,2),Range(-1,0));
	PreviousELower = 0;
};

double Cpw_fs::waveformE_value(const double& n)
{
#define ANGORA_CPW_FS_TIME_OFFSET dt
	return Data.E0*Data.waveform->Value(n*dt+OriginDelay+ANGORA_CPW_FS_TIME_OFFSET+get_initial_time_value());
	//ANGORA_CPW_FS_TIME_OFFSET=dt ensures the incident field in the main grid is coincident in time with waveformE (at the position that waveformE is applied). Remember that Einc is dt ahead of the main grid (see time chart in "timing_chart.txt"), and waveformE is coincident with Einc (because of the hard-source in Cpw_fs_upd.cpp), so waveformE is dt ahead of the main grid. Therefore, the first element of the waveformE array effectively corrects the main-grid E-field that is dt ahead in time. The inclusion of this time offset ensures that the E-field in the main grid, and not necessarily waveformE, has the desired shape.
}

void Cpw_fs::WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angle (THETA and PHI) of the scattered PW into PW_THETA and PW_PHI
	PW_THETA.resizeAndPreserve(PW_THETA.size()+1);
	PW_THETA(PW_THETA.size()-1) = M_PI - THETA_INC;		//there is no scattered field, bu the plane wave itself can be counted as the scattered PW
														//NOTE: Different from 2D: THETA is still between 0 and pi
	PW_PHI.resizeAndPreserve(PW_PHI.size()+1);
	PW_PHI(PW_PHI.size()-1) = M_PI + PHI_INC;		//there is no scattered field, bu the plane wave itself can be counted as the scattered PW
													//scattered PW goes into the diagonally-opposite quadrant of the xy plane
}

void Cpw_fs::WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delay (from the origin) of the scattered PW into origindelay_array
	TinyVector<double,3> FFOrigin;
	FFOrigin = (FFOriginX-1),(FFOriginY-1),(FFOriginZ-1);
	double FFdelay = dot(-k_inc,PWOriginationPoint-FFOrigin)*dx/vp_space;	//calculate the delay w.r.t the TRANSFORMER origin
	origindelay_array.resizeAndPreserve(origindelay_array.size()+1);
	origindelay_array(origindelay_array.size()-1) = -FFdelay;	//note the (-) sign, since "E_x_array" and "E_y_array" are actually advanced in time compared to the origin, rather than delayed
// if (rank==0) cout << "In Cpw_fs, PWOriginationPoint is " << PWOriginationPoint << endl;
// if (rank==0) cout << "In Cpw_fs, FFOrigin is " << FFOrigin << endl;
// if (rank==0) cout << "In Cpw_fs, (-k_inc) is " << -k_inc << endl;
// 	cout << origindelay_array(origindelay_array.size()-1) << endl;
}

void Cpw_fs::WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitude of the scattered PWs into field_array
	E_x_array.resizeAndPreserve(E_x_array.size()+1);
	E_y_array.resizeAndPreserve(E_y_array.size()+1);
	//write the TANGENTIAL components (x and y) of the E-field
	E_x_array(E_x_array.size()-1) = AuxGrid_Ex(HardSourcePoint);	//x-component of the E-field
	E_y_array(E_y_array.size()-1) = AuxGrid_Ey(HardSourcePoint);	//y-component of the E-field
}

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

//Definition of the class "Cpw_ml" for a TF/SF plane-wave source in a multilayered medium

#include "headers.h"

#include "Cpw_ml.h"

//for the definition of MaterialId
#include "material_id.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//Uses TinyVector operations
#include <blitz/tinyvec-et.h>

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

// //Uses Gaussian-type waveforms
// #include "waveforms/waveforms.h"

extern double courant,dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
extern int NSTEPS;

extern const int PEC;

extern Array<ElectricMaterialIndexType_X,1> Layering_e_x;
extern Array<ElectricMaterialIndexType_Y,1> Layering_e_y;
extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
extern Array<MagneticMaterialIndexType_X,1> Layering_h_x;
extern Array<MagneticMaterialIndexType_Y,1> Layering_h_y;
extern Array<MagneticMaterialIndexType_Z,1> Layering_h_z;

extern Array<double,1> eps_x,eps_z;
extern Array<double,1> cond_e_x,cond_e_z;

extern int number_of_layers;
extern Array<int,1> LayerLowerZIndices;

extern double epsilon_r_upper,epsilon_r_lower;
extern double c_upper,c_lower;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;
extern int jleft,jright;
extern int klower,kupper;

extern "C" int zgtsv_(int *n, int *nrhs,
	double *dl_r, double *dl_i,
	double *d__r, double *d__i,
	double *du_r, double *du_i,
	double *b_r, double *b_i,
	int *ldb, int *info);

extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;

extern void MPI_exit(const int& exitcode);


Cpw_ml::Cpw_ml(const PWDataType& MyData)
		: Cpw(MyData)	//calls the base class constructor first
{
	/** REMOVE LATER **/
	//check if the plane wave is incident from the upper half space (lower hs incidence is not yet supported)
	if (k_inc_z>=0) //k_inc_z inherited from Cpw
	{
		if (rank==0)
		{
			cout << "Error: Cpw_ml does not yet support plane waves incident from the lower half space." << endl << endl;
		}
		MPI_exit(-1);
	}
	/** REMOVE LATER **/

	//Unit vector that points in the direction of the TE electric field Eh (z x k_inc_lateral), perpendicular to the principal plane
	k_Eh = cross(unit_z,k_inc_lateral);
	//Unit vector that points in the direction of the TE magnetic field Hh (z x k_Eh), contained in the principal plane, perpendicular to z
	k_Hh = -k_inc_lateral;
	//Unit vector that points in the direction of the TM magnetic field He (z x k_inc_lateral = k_Eh), perpendicular to the principal plane
	k_He = k_Eh;
	//Unit vector that points in the direction of the TM electric field Ee (-z x k_He = -k_Hh), contained in the principal plane, perpendicular to z
	k_Ee = -k_Hh;

//********* This portion of the code may need further testing ******************//
	//The auxiliary grid covers enough vertical extent to include all the layers
	//lowermost cell in the auxiliary grid
	AuxGridLower = max(1,min(PWLowerZ-1,LayerLowerZIndices(1)-2));
	//1 cell needed because the H-field just outside the TF/SF box must have a counterpart in the auxiliary grid, and must be in the TF region of the auxiliary grid
	//uppermost cell in the auxiliary grid
	AuxGridUpper = min(NCELLS_Z + 2*NPML,max(PWUpperZ+4,LayerLowerZIndices(number_of_layers-1)+4));
	//a) 1 cell needed because the H-field just outside the TF/SF box must have a counterpart in the auxiliary grid, and must be in the TF region of the auxiliary grid
	//b) 3 more cells are needed because the uppermost 3 E-field components must all be in the SF region of the auxiliary grid, for the 2nd order ABC to work

	AbsorbPointLower = AuxGridLower;	//E-field absorption point (lower)
	AbsorbPointUpper = AuxGridUpper+1;	//E-field absorption point (upper)
	FeedPoint = AbsorbPointUpper-3;		//E-field feed point is 3 cells below the uppermost E-field position
						//(H-field feed point is in the same cell- right above it)
//********* This portion of the code may need further testing ******************//

	//Determine the coordinates of the plane-wave origination point
	//(measured from the back-left-bottom point of the grid)
	if ((k_inc_lateral(firstDim)<=0)&&(k_inc_lateral(secondDim)<=0))
	{
		PWOriginationPoint = (PWFrontX+1),(PWRightY+1),(FeedPoint-1);
	}
	if ((k_inc_lateral(firstDim)>0)&&(k_inc_lateral(secondDim)<=0))
	{
		PWOriginationPoint = (PWBackX-2),(PWRightY+1),(FeedPoint-1);
	}
	if ((k_inc_lateral(firstDim)>0)&&(k_inc_lateral(secondDim)>0))
	{
		PWOriginationPoint = (PWBackX-2),(PWLeftY-2),(FeedPoint-1);
	}
	if ((k_inc_lateral(firstDim)<=0)&&(k_inc_lateral(secondDim)>0))
	{
		PWOriginationPoint = (PWFrontX+1),(PWLeftY-2),(FeedPoint-1);
	}

	//components of PWOriginationPoint
	PWOriginationPoint_x = PWOriginationPoint(firstDim);
	PWOriginationPoint_y = PWOriginationPoint(secondDim);
	PWOriginationPoint_z = PWOriginationPoint(thirdDim);

	//the coordinates of the plane-wave origin
	//(measured from the back-left-bottom point of the grid)
	PWOrigin = (Data.PWOriginX-1),(Data.PWOriginY-1),(Data.PWOriginZ-1);

	//calculate grid velocity in the uppermost layer
	vp_upper = c_upper*GridVelocity(THETA_INC,PHI_INC);
// 	vp_upper = c_upper;

	//time needed for the plane wave to reach the origin from the plane-wave origination point
	OriginDelay = dot(-k_inc,PWOriginationPoint-PWOrigin)*dx/vp_upper;

	//if the there is a sudden jump in the excitation waveform at the plane-wave origin, set the initial time value back
	if ((Data.waveform->starting_time()-OriginDelay)<get_initial_time_value())
	{
		set_initial_time_value(Data.waveform->starting_time()-OriginDelay);
	}

	//fill the relative permittivity and electric_conductivity arrays (for E updates)
	epsilon_E.resize(Range(AuxGridLower,AuxGridUpper+1));
	cond_E.resize(Range(AuxGridLower,AuxGridUpper+1));
	for (int k=AuxGridLower;k<=AuxGridUpper+1;k++)
	{
		epsilon_E(k) = eps_x(Layering_e_x(k)); /** isotropy assumed! **/	//eps_y could also be used
		cond_E(k) = cond_e_x(Layering_e_x(k)); /** isotropy assumed! **/	//cond_y could also be used
	}
	//fill the relative permittivity and electric_conductivity arrays (for H updates)
	epsilon_H.resize(Range(AuxGridLower,AuxGridUpper));
	cond_H.resize(Range(AuxGridLower,AuxGridUpper));
	for (int k=AuxGridLower;k<=AuxGridUpper;k++)
	{
		epsilon_H(k) = eps_z(Layering_e_z(k)); /** isotropy assumed! **/
		cond_H(k) = cond_e_z(Layering_e_z(k)); /** isotropy assumed! **/
	}

	//fill the evanescent wave marker array (is the wave evanescent at this z-position?)
	EvanescentWavePresent = false;
	//first, for horizontal field values at integer space locations
	IsEvanescent_xy.resize(Range(AuxGridLower,AuxGridUpper+1));
	IsEvanescent_xy = false;
	for (int k=AuxGridLower;k<=AuxGridUpper+1;k++)
	{
		if ((epsilon_E(k)-epsilon_r_upper*pow(SinT,2))<=0)
		{
			IsEvanescent_xy(k)=true;	//the incident wave is evanescent at this z-position
			EvanescentWavePresent = true;
		}
	}
	//second, for vertical field values at half-integer space locations
	IsEvanescent_z.resize(Range(AuxGridLower,AuxGridUpper));
	IsEvanescent_z = false;
	for (int k=AuxGridLower;k<=AuxGridUpper;k++)
	{
		if ((epsilon_H(k)-epsilon_r_upper*pow(SinT,2))<=0)
		{
			IsEvanescent_z(k)=true;	//the incident wave is evanescent at this z-position
			EvanescentWavePresent = true;
		}
	}

	if (!dB_accuracy_set_in_configfile)
	{
		if (EvanescentWavePresent)
		{
			if (dB_accuracy<-40)
			{
				dB_accuracy = -40;		//accuracy is -40dB if evanescent waves are present
			}
		}
		else
		{
			if (dB_accuracy<-60)
			{
				dB_accuracy = -60;		//accuracy is -60dB if plane wave is propagating everywhere
			}
		}
	}

	//determine the evanescent regions
	DetermineEvanescentRegions();

	//Calculate the time step dt_g used in the auxiliary grid
	calculate_auxgrid_time_step();

	//Hilbert transform of the incident E-field waveform
	waveformHilbert = Cwf_shared_ptr(Data.waveform->HilbertTransform());

	//fill the coefficient arrays for the Ee,Eh,Hz updates (since Hz is collocated with Ee,Eh)
	Ca_e.resize(Range(AuxGridLower,AuxGridUpper+1));
	Cb_e.resize(Range(AuxGridLower,AuxGridUpper+1));
	Ca_h.resize(Range(AuxGridLower,AuxGridUpper+1));
	Cb_h.resize(Range(AuxGridLower,AuxGridUpper+1));
	Hz_h.resize(Range(AuxGridLower,AuxGridUpper+1));
	for (int k=AuxGridLower;k<=AuxGridUpper+1;k++)
	{
		if (Layering_e_x(k)!=PEC) /** isotropy assumed! **/
		{
			Ca_e(k) = (1.0-cond_E(k)*dt_g/(2*epsilon_0*epsilon_E(k)))
				/(1.0+cond_E(k)*dt_g/(2*epsilon_0*epsilon_E(k)));
			Cb_e(k) = (dt_g/dx/epsilon_0/epsilon_E(k))/(1.0+cond_E(k)*dt_g/(2*epsilon_0*epsilon_E(k)));
			Ca_h(k) = (1.0-cond_E(k)*dt_g/(2*epsilon_0*(epsilon_E(k)-epsilon_r_upper*pow(SinT,2))))
				/(1.0+cond_E(k)*dt_g/(2*epsilon_0*(epsilon_E(k)-epsilon_r_upper*pow(SinT,2))));
			Cb_h(k) = (dt_g/dx/epsilon_0)/(epsilon_E(k)-epsilon_r_upper*pow(SinT,2))
				/(1.0+cond_E(k)*dt_g/(2*epsilon_0*(epsilon_E(k)-epsilon_r_upper*pow(SinT,2))));
			Hz_h(k) = (1/eta_0)*sqrt(epsilon_r_upper)*abs(SinT);	//the coefficient Hz/Eh
		}
		else
		{
			Ca_e(k) = -1.0;
			Cb_e(k) = 0;
			Ca_h(k) = -1.0;
			Cb_h(k) = 0;
		}
	}

	//fill the coefficient arrays for the He,Hh,Ez updates (since Ez is collocated with He,Hh)
	Da_e.resize(Range(AuxGridLower,AuxGridUpper));
	Db_e.resize(Range(AuxGridLower,AuxGridUpper));
	Da_e_pr.resize(Range(AuxGridLower,AuxGridUpper));
	Db_e_pr.resize(Range(AuxGridLower,AuxGridUpper));
	Da_h.resize(Range(AuxGridLower,AuxGridUpper));
	Db_h.resize(Range(AuxGridLower,AuxGridUpper));
	Ez_e.resize(Range(AuxGridLower,AuxGridUpper));
	Ez_e_pr.resize(Range(AuxGridLower,AuxGridUpper));
	//coefficients for the auxiliary variable He'
	A_pr.resize(Range(AuxGridLower,AuxGridUpper));
	B_pr.resize(Range(AuxGridLower,AuxGridUpper));
	for (int k=AuxGridLower;k<=AuxGridUpper;k++)
	{
		Da_e(k) = 1.0;
		Db_e(k) = (dt_g/dx/mu_0*epsilon_H(k))/(epsilon_H(k)-epsilon_r_upper*pow(SinT,2));
		Da_e_pr(k) = (1.0-cond_H(k)*dt_g/(2*epsilon_0*(epsilon_H(k)-epsilon_r_upper*pow(SinT,2))))
				/(1.0+cond_H(k)*dt_g/(2*epsilon_0*(epsilon_H(k)-epsilon_r_upper*pow(SinT,2))));
		Db_e_pr(k) = (dt_g/dx/mu_0)/(epsilon_H(k)-epsilon_r_upper*pow(SinT,2))
				/(1.0+cond_H(k)*dt_g/(2*epsilon_0*(epsilon_H(k)-epsilon_r_upper*pow(SinT,2))));
		Da_h(k) = 1.0;
		Db_h(k) = dt_g/dx/mu_0;
		A_pr(k) = epsilon_H(k)+cond_H(k)*dt_g/(2*epsilon_0);
		B_pr(k) = epsilon_H(k)-cond_H(k)*dt_g/(2*epsilon_0);
		Ez_e(k) = -eta_0*sqrt(epsilon_r_upper)*abs(SinT)/epsilon_H(k);	//the coefficient Ez/He
		Ez_e_pr(k) = -eta_0*sqrt(epsilon_r_upper)*abs(SinT);	//the coefficient Ez/He'
	}

	//Auxiliary 1-D grids for incident field computation
	AuxGrid_Ee.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_Eh.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_Ez.resize(Range(AuxGridLower,AuxGridUpper));
	AuxGrid_He.resize(Range(AuxGridLower-1,AuxGridUpper+1));	//consistent with Hx,Hy, extend 1/2 cell outside last E
	AuxGrid_Hh.resize(Range(AuxGridLower-1,AuxGridUpper+1));	//consistent with Hx,Hy, extend 1/2 cell outside last E
	AuxGrid_Hz.resize(Range(AuxGridLower,AuxGridUpper+1));
	AuxGrid_He_pr.resize(Range(AuxGridLower-1,AuxGridUpper+1));
	AuxGrid_Ee = 0;
	AuxGrid_Eh = 0;
	AuxGrid_Ez = 0;
	AuxGrid_He = 0;
	AuxGrid_Hh = 0;
	AuxGrid_Hz = 0;
	AuxGrid_He_pr = 0;

	//effective courant number for 1D grid:
	Sc_1D_upper = c/sqrt(epsilon_E(AbsorbPointUpper)-epsilon_r_upper*pow(SinT,2))*dt_g/dx;
	Sc_1D_lower = c/sqrt(epsilon_E(AbsorbPointLower)-epsilon_r_upper*pow(SinT,2))*dt_g/dx;

	//resize and initialize ABC storage variables PreviousEeLower(k,n) etc.
	//NOTE: k=0,n=0 need not be used, since this is the current value of the E-field at the boundary of the auxiliary grid.
	//However, it will still be used for symmetry (see Schneider course notes on 2nd order ABC).
	PreviousEeLower.resize(Range(0,2),Range(-1,0));
	PreviousEeUpper.resize(Range(0,2),Range(-1,0));
	PreviousEhLower.resize(Range(0,2),Range(-1,0));
	PreviousEhUpper.resize(Range(0,2),Range(-1,0));
	PreviousEeLower = 0;
	PreviousEeUpper = 0;
	PreviousEhLower = 0;
	PreviousEhUpper = 0;

	//Time histories of field components in the 1-D grid
	TimeHistory_Ee.resize(Range(0,NSTEPS),Range(klower,kupper+1));
	TimeHistory_Eh.resize(Range(0,NSTEPS),Range(klower,kupper+1));
	TimeHistory_Ez.resize(Range(0,NSTEPS),Range(klower,kupper));
	TimeHistory_He.resize(Range(0,NSTEPS),Range(klower-1,kupper+1));	//see note for AuxGrid_He range
	TimeHistory_Hh.resize(Range(0,NSTEPS),Range(klower-1,kupper+1));	//see note for AuxGrid_He range
	TimeHistory_Hz.resize(Range(0,NSTEPS),Range(klower,kupper+1));
	TimeHistory_Ee = 0;
	TimeHistory_Eh = 0;
	TimeHistory_Ez = 0;
	TimeHistory_He = 0;
	TimeHistory_Hh = 0;
	TimeHistory_Hz = 0;
};

void Cpw_ml::calculate_auxgrid_time_step()
//calculates the reduced time step to be used in the auxiliary time grid
{
	//First, find the maximum and minimum velocity of propagation (c_min, c_max) in the auxiliary grid
	//c_min limits the maximum frequency in the incident wave (for dispersion)
	//c_max maximum allowable real time step in the auxiliary grid (for stability)
	double c_min=c;
	double c_max=0;
	//go through horizontal field values at integer space locations,
	//and look at the half-integer locations above and below that position.
	double c1,c2,c3,c4;	//c1,c2 for TM, c3 for TE, c4 for the main grid velocity
	for (int k=AuxGridLower+1; k<=AuxGridUpper;k++)
	{
		c1 = sqrt(1/epsilon_0/epsilon_E(k)
			*epsilon_H(k)/mu_0/abs(epsilon_H(k)-epsilon_r_upper*pow(SinT,2)));
		c2 = sqrt(1/epsilon_0/epsilon_E(k)
			*epsilon_H(k-1)/mu_0/abs(epsilon_H(k-1)-epsilon_r_upper*pow(SinT,2)));
		c3 = sqrt(1/epsilon_0/abs(epsilon_E(k)-epsilon_r_upper*pow(SinT,2))
			*1/mu_0);
		c4 = c/sqrt(epsilon_E(k));
		if (c1>c_max)
		{
			c_max = c1;
		} else if (c2>c_max)
		{
			c_max = c2;
		} else if (c3>c_max)
		{
			c_max = c3;
		}
		else if (c4<c_min)
		{
			c_min = c4;
		}
	}

	//******************************//
	//		FIND TIME STEP dt_g     //
	//******************************//
	//find maximum allowable dt_g:
	if (EvanescentWavePresent)	//evanescent waves introduce more stringent stability conditions
	{
		double courant_aux = 0.5;	//more accurate upper bounds can be found later
		dt_g = courant_aux*dx/c_max;
		if (dt_g<1/Data.waveform->W_0())
		{
			dt_g = pow(dt_g,2)*Data.waveform->W_0();		//if the new stability criterion is more strict, reduce the time step
		}
	}
	else
	{
		double courant_aux = courant;	//if there are no evanescent waves, courant condition is sufficient
		dt_g = courant_aux*dx/c_max;
	}
	//Now that the time step dt_g is determined,
	//make sure dt_g is also dt/n, where n is an odd integer
	timefactor = (int)(ceil(dt/dt_g)+1e-7);
	if (timefactor==2*(timefactor/2))	//if not odd
	{
		timefactor++;	//make it odd
	}
	//dt must be an odd factor of dt_g:
	dt_g = dt/timefactor;	//smaller time step used in auxiliary grid
	//******************************//
	//		TIME STEP dt_g FOUND    //
	//******************************//

//	int aux_grid_length = (NSTEPS+1)*timefactor+1;
//
//	waveformE_analytic.resize(aux_grid_length);
//	waveformH_analytic.resize(aux_grid_length);
//	waveformE_analytic = 0;
//	waveformH_analytic = 0;
//
//	double distance = dx/2*CosT;	//distance traveled from the E-feed point to the H-feed point
//	double time_offset = dt_g;	//time offset that ensures the incident field in the main grid is coincident in time with waveformE (at the position that waveformE is applied). Remember that Einc is dt ahead of the main grid (see time chart in "timing_chart.txt"), and waveformE is (dt-dt_g) behind Einc
////(see time chart in "timing_chart.txt"), so waveformE is dt_g ahead of the main grid. Therefore, the first element of the waveformE array effectively corrects the main-grid E-field that is dt_g ahead in time. The inclusion of this time offset ensures that the E-field in the main grid, and not necessarily waveformE, has the desired shape.
//
//	Cwf_shared_ptr waveformHilbert(Data.waveform->HilbertTransform()); //cleanup done by waveform automatically
//	for (int n=0;n<aux_grid_length;n++)
//	{
//		waveformE_analytic(n)=Data.E0*
//			(Data.waveform->Value(n*dt_g+OriginDelay+time_offset)
//			 +ii*(waveformHilbert->Value(n*dt_g+OriginDelay+time_offset))); //imag. part is Hilb. transform
//		waveformH_analytic(n)=Data.E0*1/(eta_0/sqrt(epsilon_r_upper))*
//			(Data.waveform->Value((n-0.5)*dt_g+OriginDelay+distance/vp_upper+time_offset)
//			+ii*(waveformHilbert->Value((n-0.5)*dt_g+OriginDelay+distance/vp_upper+time_offset))); //imag. part is Hilb. transform
//	}
//
//	//Now, the incident time arrays Incident_Ee, Incident_Eh, Incident_He, Incident_Hh can be defined
//	//+1 extra real time step required for interpolation in some instances
//	Incident_Ee.resize(aux_grid_length);
//	Incident_Eh.resize(aux_grid_length);
//	Incident_He.resize(aux_grid_length);
//	Incident_Hh.resize(aux_grid_length);
//	Incident_Ee=0;
//	Incident_Eh=0;
//	Incident_He=0;
//	Incident_Hh=0;
//
//	Range CommonRange = Range(0,aux_grid_length-1);
//	Incident_Ee(CommonRange)= SinPsi*CosT*waveformE_analytic(CommonRange);
//	Incident_Eh(CommonRange)= -CosPsi*waveformE_analytic(CommonRange);
//	Incident_He(CommonRange)= -SinPsi*waveformH_analytic(CommonRange);
//	Incident_Hh(CommonRange)= CosPsi*CosT*waveformH_analytic(CommonRange);
}

//ANGORA_CPW_ML_TIME_OFFSET ensures the incident field in the main grid is coincident in time with waveformE (at the position that waveformE is applied). Remember that Einc is dt ahead of the main grid (see time chart in "timing_chart.txt"), and waveformE is (dt-dt_g) behind Einc
////(see time chart in "timing_chart.txt"), so waveformE is dt_g ahead of the main grid. Therefore, the first element of the waveformE array effectively corrects the main-grid E-field that is dt_g ahead in time. The inclusion of this time offset ensures that the E-field in the main grid, and not necessarily waveformE, has the desired shape.
#define ANGORA_CPW_ML_TIME_OFFSET dt_g

//distance traveled from the E-feed point to the H-feed point
#define ANGORA_CPW_ML_SPATIAL_OFFSET dx/2*CosT

complex<double> Cpw_ml::Incident_Ee(const double& m)
{
	return SinPsi*CosT*Data.E0*(Data.waveform->Value(m*dt_g+OriginDelay+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())
								+ii*(waveformHilbert->Value(m*dt_g+OriginDelay+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())));
}

complex<double> Cpw_ml::Incident_Eh(const double& m)
{
	return -CosPsi*Data.E0*(Data.waveform->Value(m*dt_g+OriginDelay+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())
								+ii*(waveformHilbert->Value(m*dt_g+OriginDelay+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())));
}

complex<double> Cpw_ml::Incident_He(const double& m)
{
	return -SinPsi*Data.E0*1/(eta_0/sqrt(epsilon_r_upper))*
	(Data.waveform->Value((m-0.5)*dt_g+OriginDelay+ANGORA_CPW_ML_SPATIAL_OFFSET/vp_upper+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())
 +ii*(waveformHilbert->Value((m-0.5)*dt_g+OriginDelay+ANGORA_CPW_ML_SPATIAL_OFFSET/vp_upper+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())));
}

complex<double> Cpw_ml::Incident_Hh(const double& m)
{
	return CosPsi*CosT*Data.E0*1/(eta_0/sqrt(epsilon_r_upper))*
	(Data.waveform->Value((m-0.5)*dt_g+OriginDelay+ANGORA_CPW_ML_SPATIAL_OFFSET/vp_upper+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())
 +ii*(waveformHilbert->Value((m-0.5)*dt_g+OriginDelay+ANGORA_CPW_ML_SPATIAL_OFFSET/vp_upper+ANGORA_CPW_ML_TIME_OFFSET+get_initial_time_value())));
}

void Cpw_ml::DetermineEvanescentRegions()
{
	NumOfEvanescentLayers = 0;
	EvanescentLayers.resize(0);
	int k = AuxGridLower;
	while (k<=AuxGridUpper)
	{
		if (IsEvanescent_z(k))	//start block
		{
			NumOfEvanescentLayers++;	//increment number of evanescent layers
			EvanescentLayers.resizeAndPreserve(EvanescentLayers.size()+1);
			EvanescentLayers(NumOfEvanescentLayers-1).lowerE = k+1;
			EvanescentLayers(NumOfEvanescentLayers-1).lowerH = k;
			if (!IsEvanescent_xy(k))
			{
				EvanescentLayers(NumOfEvanescentLayers-1).homogeneous_at_interfaces = true;
			}
			else
			{
				EvanescentLayers(NumOfEvanescentLayers-1).homogeneous_at_interfaces = false;
			}
			for(;IsEvanescent_z(k+1);k++);	//increment k until the end of the evanescent layer
			EvanescentLayers(NumOfEvanescentLayers-1).upperE = k;
			EvanescentLayers(NumOfEvanescentLayers-1).upperH = k;
		}
		k++;
	}

	for (int i=0; i<NumOfEvanescentLayers; i++)
	{

		//dimension of the tridiagonal matrix (for E updates)
		EvanescentLayers(i).N_e = EvanescentLayers(i).upperE-EvanescentLayers(i).lowerE+1;
		//dimension of the tridiagonal matrix (for H updates)
		EvanescentLayers(i).N_h = EvanescentLayers(i).upperH-EvanescentLayers(i).lowerH+1;
		EvanescentLayers(i).nrhs = 1;

		// The tridiagonal matrices used in the inhomogeneous region are allocated here.
		EvanescentLayers(i).diagE_r.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE));
		EvanescentLayers(i).diagE_i.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE));
		EvanescentLayers(i).lowdiagE_r.resize(Range(EvanescentLayers(i).lowerE+1,EvanescentLayers(i).upperE));
		EvanescentLayers(i).lowdiagE_i.resize(Range(EvanescentLayers(i).lowerE+1,EvanescentLayers(i).upperE));
		EvanescentLayers(i).updiagE_r.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE-1));
		EvanescentLayers(i).updiagE_i.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE-1));
		EvanescentLayers(i).righthandE_r.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE));
		EvanescentLayers(i).righthandE_i.resize(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE));

		EvanescentLayers(i).diagH_r.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH));
		EvanescentLayers(i).diagH_i.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH));
		EvanescentLayers(i).lowdiagH_r.resize(Range(EvanescentLayers(i).lowerH+1,EvanescentLayers(i).upperH));
		EvanescentLayers(i).lowdiagH_i.resize(Range(EvanescentLayers(i).lowerH+1,EvanescentLayers(i).upperH));
		EvanescentLayers(i).updiagH_r.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH-1));
		EvanescentLayers(i).updiagH_i.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH-1));
		EvanescentLayers(i).righthandH_r.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH));
		EvanescentLayers(i).righthandH_i.resize(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH));
		// Tridiagonal matrices initialized.
	}
}

void Cpw_ml::WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angle (THETA and PHI) of the scattered PW into PW_THETA and PW_PHI
	//first, record the specular reflection angles
	PW_THETA.resizeAndPreserve(PW_THETA.size()+1);
	PW_THETA(PW_THETA.size()-1) = -THETA_INC;
	PW_PHI.resizeAndPreserve(PW_PHI.size()+1);
	PW_PHI(PW_PHI.size()-1) = PHI_INC;

	//then, record the transmission angle at the lower half space
	if (abs(sin(THETA_INC)/sqrt(epsilon_r_lower/epsilon_r_upper))<1)
	{//record if the plane wave is not evanescent in the lowermost layer
		//theta angle
		PW_THETA.resizeAndPreserve(PW_THETA.size()+1);
		// calculate the transmission theta angle using Snell's law
		double theta_trans = asin(sin(THETA_INC)/sqrt(epsilon_r_lower/epsilon_r_upper));
		PW_THETA(PW_THETA.size()-1) = M_PI+theta_trans;	//transmission angle is measured clockwise around the -z axis
		//phi angle
		PW_PHI.resizeAndPreserve(PW_PHI.size()+1);
		PW_PHI(PW_PHI.size()-1) = PHI_INC;	//transmission phi angle is the same as the incident phi angle
	}
}

void Cpw_ml::WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delay (from the origin) of the scattered PW into origindelay_array
	//first, for the specular reflection angle
	TinyVector<double,3> FFOrigin,PWScatteredFieldPoint;
	FFOrigin = (FFOriginX-1),(FFOriginY-1),(FFOriginZ-1);
	PWScatteredFieldPoint = (PWOriginationPoint(firstDim)),(PWOriginationPoint(secondDim)),(FeedPoint+1);	//go into the scattered-field (SF) region
	TinyVector<double,3> k_spec(k_inc);	//specular reflection
	k_spec(thirdDim) *= -1;	//reverse the z-component
	double FFdelay = dot(k_spec,PWScatteredFieldPoint-FFOrigin)*dx/vp_upper;	//calculate the delay w.r.t the TRANSFORMER origin
// 	cout << k_inc << endl;
// 	cout << PWScatteredFieldPoint << endl;
// 	cout << FFOrigin << endl;
// 	cout << dot(-k_inc,PWScatteredFieldPoint-FFOrigin) << endl;
// 	cout << FFdelay << endl;
	origindelay_array.resizeAndPreserve(origindelay_array.size()+1);
	origindelay_array(origindelay_array.size()-1) = FFdelay;

	//then, for the transmission angle at the lower half space
	if (abs(sin(THETA_INC)/sqrt(epsilon_r_lower/epsilon_r_upper))<1)
	{//record if the plane wave is not evanescent in the lowermost layer
		PWScatteredFieldPoint = (PWOriginationPoint(firstDim)),(PWOriginationPoint(secondDim)),(AbsorbPointLower-1);	//take the field value at the absorption point
		//transmission angle is measured clockwise around the -z axis
		double theta_trans = asin(sin(THETA_INC)/sqrt(epsilon_r_lower/epsilon_r_upper));		//get the transmission angle
		TinyVector<double,3> k_trans;
		k_trans(firstDim) = -sin(theta_trans)*cos(PHI_INC);
		k_trans(secondDim) = -sin(theta_trans)*sin(PHI_INC);
		k_trans(thirdDim) = -cos(theta_trans);
		FFdelay = dot(k_trans,PWScatteredFieldPoint-FFOrigin)*dx/c_lower;	//calculate the delay w.r.t the TRANSFORMER origin
		origindelay_array.resizeAndPreserve(origindelay_array.size()+1);
		origindelay_array(origindelay_array.size()-1) = FFdelay;	//note the (-) sign, since "field_array" is actually advanced in time compared to the origin, rather than delayed
	}
}

void Cpw_ml::WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitude of the scattered PWs into field_array
	//first, for the specular reflection
	E_x_array.resizeAndPreserve(E_x_array.size()+1);
	E_y_array.resizeAndPreserve(E_y_array.size()+1);
	//write the TANGENTIAL component (either x or y) of the E-field
	//x-component of the E-field
	E_x_array(E_x_array.size()-1) = k_Ee(firstDim)*real(AuxGrid_Ee(FeedPoint+2))
								+ k_Eh(firstDim)*real(AuxGrid_Eh(FeedPoint+2));
	//y-component of the E-field
	E_y_array(E_y_array.size()-1) = k_Ee(secondDim)*real(AuxGrid_Ee(FeedPoint+2))
								+ k_Eh(secondDim)*real(AuxGrid_Eh(FeedPoint+2));

	//then, for the transmission angle at the lower half space
	if (abs(sin(THETA_INC)/sqrt(epsilon_r_lower/epsilon_r_upper))<1)
	{//record if the plane wave is not evanescent in the lowermost layer
		E_x_array.resizeAndPreserve(E_x_array.size()+1);
		E_y_array.resizeAndPreserve(E_y_array.size()+1);
		//write the TANGENTIAL component (either x or y) of the E-field
		//x-component of the E-field
		E_x_array(E_x_array.size()-1) = k_Ee(firstDim)*real(AuxGrid_Ee(AbsorbPointLower))
									+ k_Eh(firstDim)*real(AuxGrid_Eh(AbsorbPointLower));
		//y-component of the E-field
		E_y_array(E_y_array.size()-1) = k_Ee(secondDim)*real(AuxGrid_Ee(AbsorbPointLower))
									+ k_Eh(secondDim)*real(AuxGrid_Eh(AbsorbPointLower));
	}
}

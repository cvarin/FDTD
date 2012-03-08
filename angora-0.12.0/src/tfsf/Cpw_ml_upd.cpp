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

// Includes the methods that update the auxiliary grid associated with a single plane wave in a multilayered medium.

#include "headers.h"

#include "Cpw_ml.h"

extern int rank;
extern int klower,kupper;

extern "C" int zgtsv_(int *n, int *nrhs,
	double *dl_r, double *dl_i,
	double *d__r, double *d__i,
	double *du_r, double *du_i,
	double *b_r, double *b_i,
	int *ldb, int *info);

extern int i,j,k;


void Cpw_ml::UpdateIncidentH(const int& n)
{
	//update incident E and H until H-field is at time n+1/2
	for (timeindex=0;timeindex<timefactor/2;timeindex++)
	{
		FractionalUpdate(n); //update incident fields by dt_g
	}
	//record the H field at n+1/2 (see time chart in "timing_chart.txt")
	for (k=max(AuxGridLower-1,klower-1);k<=min(AuxGridUpper+1,kupper+1);k++)
	{
		TimeHistory_He(n,k) = real(AuxGrid_He(k));
		TimeHistory_Hh(n,k) = real(AuxGrid_Hh(k));
	}
	for (k=max(AuxGridLower,klower);k<=min(AuxGridUpper+1,kupper+1);k++)
	{
		TimeHistory_Hz(n,k) = real(AuxGrid_Hz(k));
	}
}

void Cpw_ml::UpdateIncidentE(const int& n)
{
	//update incident E and H until E-field is at time n+1
	for (timeindex=timefactor/2;timeindex<timefactor;timeindex++)
	{
		FractionalUpdate(n);	//update incident fields by dt_g
	}
	//record the E field at n+1 (see time chart in "timing_chart.txt")
	for (k=max(AuxGridLower,klower);k<=min(AuxGridUpper+1,kupper+1);k++)
	{
		TimeHistory_Ee(n,k) = real(AuxGrid_Ee(k));
		TimeHistory_Eh(n,k) = real(AuxGrid_Eh(k));
	}
	for (k=max(AuxGridLower,klower);k<=min(AuxGridUpper,kupper);k++)
	{
		TimeHistory_Ez(n,k) = real(AuxGrid_Ez(k));
	}
	//The recorded E field is dt/2 ahead of the recorded H field, and they have the same time index (for ex. n=0).
	//This is desirable since E is updated first in the main grid.
	//(see the time chart in "timing_chart.txt")
}

void Cpw_ml::FractionalUpdate(const int& n)
{
	//makes ONE grid update for both E and H
	//(a full main grid update includes (timefactor)xUpdateIncident)
	//Updating sequence:
	//1)E-field components are updated.
	//2)H-field components are updated.

	//Update the E-field in propagating regions
	UpdateE(n);
	//Update the E-field in evanescent regions
	UpdateEvanescentE();

	//Update the H-field in propagating regions
	UpdateH(n);
	//Update the H-field in evanescent regions
	UpdateEvanescentH();
}

void Cpw_ml::UpdateE(const int& n)
{//updates the E-field components at time index n
	//Update Ee,Eh
	for (k=AbsorbPointLower+1;k<=AbsorbPointUpper-1;k++)
	{
		if (!IsEvanescent_xy(k))	//if wave is propagating, perform leap-frog update
		{
			//Ee
			AuxGrid_Ee(k) = Ca_e(k)*AuxGrid_Ee(k) - Cb_e(k)*(AuxGrid_He(k)-AuxGrid_He(k-1));
			//Eh
			AuxGrid_Eh(k) = Ca_h(k)*AuxGrid_Eh(k) - Cb_h(k)*(AuxGrid_Hh(k)-AuxGrid_Hh(k-1));
		}
	}
	//Update Ez
	for (k=AbsorbPointLower;k<=AbsorbPointUpper-1;k++)
	{
		if (!IsEvanescent_z(k))	//if wave is propagating, perform leap-frog update
		{
			//Ez
			AuxGrid_Ez(k) = -AuxGrid_Ez(k) + 2.0*Ez_e_pr(k)*AuxGrid_He_pr(k);
		}
	}

	//TF/SF corrections and upper ABC are applied in the uppermost layer, where the plane wave is assumed to be PROPAGATING (non-evanescent)

	//Apply (TF/SF) correction for Einc: (NOTE: Uppermost layer is assumed LOSSLESS)
	//correct Ee
	AuxGrid_Ee(FeedPoint) += -Cb_e(FeedPoint)*Incident_He(n*timefactor+timeindex);
	//correct Eh
	AuxGrid_Eh(FeedPoint) += -Cb_h(FeedPoint)*Incident_Hh(n*timefactor+timeindex);

	//Apply 2nd order ABC	(NOTE: Uppermost and lowermost layers are assumed LOSSLESS)
	//NOTE: PreviousE{e,h}{Upper,Lower}(0,0) need not be used, since it is equal to AuxGrid_E{e,h}(AbsorbPoint{Upper,Lower}) [before the update]
	//, but it is used for symmetry.
	//Upper absorbing point
	//absorb Ee
	AuxGrid_Ee(AbsorbPointUpper) = (-1/(1/Sc_1D_upper+2+Sc_1D_upper))*(
		(1/Sc_1D_upper-2+Sc_1D_upper)*(AuxGrid_Ee(AbsorbPointUpper-2)+PreviousEeUpper(0,-1))
		+ 2*(Sc_1D_upper-1/Sc_1D_upper)*(PreviousEeUpper(0,0)+PreviousEeUpper(2,0)-AuxGrid_Ee(AbsorbPointUpper-1)-PreviousEeUpper(1,-1))
		- 4*(1/Sc_1D_upper+Sc_1D_upper)*PreviousEeUpper(1,0)) - PreviousEeUpper(2,-1);
	//absorb Eh
	AuxGrid_Eh(AbsorbPointUpper) = (-1/(1/Sc_1D_upper+2+Sc_1D_upper))*(
		(1/Sc_1D_upper-2+Sc_1D_upper)*(AuxGrid_Eh(AbsorbPointUpper-2)+PreviousEhUpper(0,-1))
		+ 2*(Sc_1D_upper-1/Sc_1D_upper)*(PreviousEhUpper(0,0)+PreviousEhUpper(2,0)-AuxGrid_Eh(AbsorbPointUpper-1)-PreviousEhUpper(1,-1))
		- 4*(1/Sc_1D_upper+Sc_1D_upper)*PreviousEhUpper(1,0)) - PreviousEhUpper(2,-1);
	//update the field history values
	PreviousEeUpper(Range(0,2),-1) = PreviousEeUpper(Range(0,2),0);	//shift toward the future by one step
	PreviousEeUpper(Range(0,2),0) = AuxGrid_Ee(Range(AbsorbPointUpper,AbsorbPointUpper-2,-1));	//record the current values as history
	PreviousEhUpper(Range(0,2),-1) = PreviousEhUpper(Range(0,2),0);	//shift toward the future by one step
	PreviousEhUpper(Range(0,2),0) = AuxGrid_Eh(Range(AbsorbPointUpper,AbsorbPointUpper-2,-1));	//record the current values as history

	//Lower absorbing point	(assuming the wave is PROPAGATING at this point !!!! May generalize to evanescent waves later)
	if (!IsEvanescent_xy(AbsorbPointLower))	//if wave is not evanescent at either below or above this position, perform absorption
	{
		//absorb Ee
		AuxGrid_Ee(AbsorbPointLower) = (-1/(1/Sc_1D_lower+2+Sc_1D_lower))*(
			(1/Sc_1D_lower-2+Sc_1D_lower)*(AuxGrid_Ee(AbsorbPointLower+2)+PreviousEeLower(0,-1))
			+ 2*(Sc_1D_lower-1/Sc_1D_lower)*(PreviousEeLower(0,0)+PreviousEeLower(2,0)-AuxGrid_Ee(AbsorbPointLower+1)-PreviousEeLower(1,-1))
			- 4*(1/Sc_1D_lower+Sc_1D_lower)*PreviousEeLower(1,0)) - PreviousEeLower(2,-1);
		//absorb Eh
		AuxGrid_Eh(AbsorbPointLower) = (-1/(1/Sc_1D_lower+2+Sc_1D_lower))*(
			(1/Sc_1D_lower-2+Sc_1D_lower)*(AuxGrid_Eh(AbsorbPointLower+2)+PreviousEhLower(0,-1))
			+ 2*(Sc_1D_lower-1/Sc_1D_lower)*(PreviousEhLower(0,0)+PreviousEhLower(2,0)-AuxGrid_Eh(AbsorbPointLower+1)-PreviousEhLower(1,-1))
			- 4*(1/Sc_1D_lower+Sc_1D_lower)*PreviousEhLower(1,0)) - PreviousEhLower(2,-1);
		//update the field history values
		PreviousEeLower(Range(0,2),-1) = PreviousEeLower(Range(0,2),0);	//shift toward the future by one step
		PreviousEeLower(Range(0,2),0) = AuxGrid_Ee(Range(AbsorbPointLower,AbsorbPointLower+2));	//record the current values as history
		PreviousEhLower(Range(0,2),-1) = PreviousEhLower(Range(0,2),0);	//shift toward the future by one step
		PreviousEhLower(Range(0,2),0) = AuxGrid_Eh(Range(AbsorbPointLower,AbsorbPointLower+2));	//record the current values as history
	}
}

void Cpw_ml::UpdateH(const int& n)
{//updates the H-field components at time index n
	//Update He,Hh
	for (k=AbsorbPointLower;k<=AbsorbPointUpper-1;k++)
	{
		if (!IsEvanescent_z(k))		//if wave is propagating, perform leap-frog update
		{
			//update He' first
			fieldstore = AuxGrid_He_pr(k);
			AuxGrid_He_pr(k) = Da_e_pr(k)*AuxGrid_He_pr(k) - Db_e_pr(k)*(AuxGrid_Ee(k+1)-AuxGrid_Ee(k));
			//then obtain He from He'
			AuxGrid_He(k) = AuxGrid_He(k) + (A_pr(k)*AuxGrid_He_pr(k) - B_pr(k)*fieldstore);
			//finally update Hh
			AuxGrid_Hh(k) = Da_h(k)*AuxGrid_Hh(k) - Db_h(k)*(AuxGrid_Eh(k+1)-AuxGrid_Eh(k));
		}
	}
	//Update Hz
	for (k=AbsorbPointLower+1;k<=AbsorbPointUpper-1;k++)
	{
		if (!IsEvanescent_xy(k))
		{
			//Hz
			AuxGrid_Hz(k) = -AuxGrid_Hz(k) + 2.0*Hz_h(k)*AuxGrid_Eh(k);
		}
	}

	//TF/SF correction is made in the uppermost layer, where the plane wave is assumed to be PROPAGATING (non-evanescent)

	//Apply (TF/SF) correction for Hinc: (NOTE: Uppermost layer is assumed LOSSLESS)
	//correct He' first
	AuxGrid_He_pr(FeedPoint) += -Db_e_pr(FeedPoint)*Incident_Ee(n*timefactor+timeindex);
	//then obtain He from He' (linear dependence, since uppermost medium is lossless => er=er')
	AuxGrid_He(FeedPoint) = epsilon_H(FeedPoint)*AuxGrid_He_pr(FeedPoint);
	//finally correct Hh
	AuxGrid_Hh(FeedPoint) += -Db_h(FeedPoint)*Incident_Eh(n*timefactor+timeindex);
}

void Cpw_ml::UpdateEvanescentE()
{
	for (i=0; i<NumOfEvanescentLayers; i++)
	{
		if (EvanescentLayers(i).homogeneous_at_interfaces)
		{
			//Update Ee
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE; k++)
			{
				EvanescentLayers(i).diagE_r(k) = real(1.0-(Cb_e(k)*(Db_e(k)+Db_e(k-1)))/pow(w_0*dt_g,2));
				EvanescentLayers(i).diagE_i(k) = imag(1.0-(Cb_e(k)*(Db_e(k)+Db_e(k-1)))/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE-1; k++)
			{
				EvanescentLayers(i).updiagE_r(k) = real(Cb_e(k)*Db_e(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).updiagE_i(k) = imag(Cb_e(k)*Db_e(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagE_r(k+1) = real(Cb_e(k+1)*Db_e(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagE_i(k+1) = imag(Cb_e(k+1)*Db_e(k)/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE; k++)
			{
				EvanescentLayers(i).righthandE_r(k) = 0.0;
				EvanescentLayers(i).righthandE_i(k) = 0.0;
			}
			EvanescentLayers(i).righthandE_r(EvanescentLayers(i).upperE) = -real(Cb_e(EvanescentLayers(i).upperE)
				*Db_e(EvanescentLayers(i).upperE)/pow(w_0*dt_g,2)*AuxGrid_Ee(EvanescentLayers(i).upperE+1));
			EvanescentLayers(i).righthandE_i(EvanescentLayers(i).upperE) = -imag(Cb_e(EvanescentLayers(i).upperE)
				*Db_e(EvanescentLayers(i).upperE)/pow(w_0*dt_g,2)*AuxGrid_Ee(EvanescentLayers(i).upperE+1));
			EvanescentLayers(i).righthandE_r(EvanescentLayers(i).lowerE) = -real(Cb_e(EvanescentLayers(i).lowerE)
				*Db_e(EvanescentLayers(i).lowerE-1)/pow(w_0*dt_g,2)*AuxGrid_Ee(EvanescentLayers(i).lowerE-1));
			EvanescentLayers(i).righthandE_i(EvanescentLayers(i).lowerE) = -imag(Cb_e(EvanescentLayers(i).lowerE)
				*Db_e(EvanescentLayers(i).lowerE-1)/pow(w_0*dt_g,2)*AuxGrid_Ee(EvanescentLayers(i).lowerE-1));
			//solve for the E-field
			zgtsv_(&EvanescentLayers(i).N_e, &EvanescentLayers(i).nrhs, EvanescentLayers(i).lowdiagE_r.data(), EvanescentLayers(i).lowdiagE_i.data(), EvanescentLayers(i).diagE_r.data(), EvanescentLayers(i).diagE_i.data(), EvanescentLayers(i).updiagE_r.data(), EvanescentLayers(i).updiagE_i.data(), EvanescentLayers(i).righthandE_r.data(), EvanescentLayers(i).righthandE_i.data(), &EvanescentLayers(i).N_e, &info);
			if (info!=0)
			{
				cout << endl << "Matrix inversion failed!!" << endl;
				exit(-1);
			}

			//Place the computed E-field values into the auxiliary grid
			AuxGrid_Ee(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE)) = EvanescentLayers(i).righthandE_r + ii*EvanescentLayers(i).righthandE_i;

			//Update Eh
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE; k++)
			{
				EvanescentLayers(i).diagE_r(k) = real(1.0-(Cb_h(k)*(Db_h(k)+Db_h(k-1)))/pow(w_0*dt_g,2));
				EvanescentLayers(i).diagE_i(k) = imag(1.0-(Cb_h(k)*(Db_h(k)+Db_h(k-1)))/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE-1; k++)
			{
				EvanescentLayers(i).updiagE_r(k) = real(Cb_h(k)*Db_h(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).updiagE_i(k) = imag(Cb_h(k)*Db_h(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagE_r(k+1) = real(Cb_h(k+1)*Db_h(k)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagE_i(k+1) = imag(Cb_h(k+1)*Db_h(k)/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerE; k<=EvanescentLayers(i).upperE; k++)
			{
				EvanescentLayers(i).righthandE_r(k) = 0.0;
				EvanescentLayers(i).righthandE_i(k) = 0.0;
			}
			EvanescentLayers(i).righthandE_r(EvanescentLayers(i).upperE) = -real(Cb_h(EvanescentLayers(i).upperE)
				*Db_h(EvanescentLayers(i).upperE)/pow(w_0*dt_g,2)*AuxGrid_Eh(EvanescentLayers(i).upperE+1));
			EvanescentLayers(i).righthandE_i(EvanescentLayers(i).upperE) = -imag(Cb_h(EvanescentLayers(i).upperE)
				*Db_h(EvanescentLayers(i).upperE)/pow(w_0*dt_g,2)*AuxGrid_Eh(EvanescentLayers(i).upperE+1));
			EvanescentLayers(i).righthandE_r(EvanescentLayers(i).lowerE) = -real(Cb_h(EvanescentLayers(i).lowerE)
				*Db_h(EvanescentLayers(i).lowerE-1)/pow(w_0*dt_g,2)*AuxGrid_Eh(EvanescentLayers(i).lowerE-1));
			EvanescentLayers(i).righthandE_i(EvanescentLayers(i).lowerE) = -imag(Cb_h(EvanescentLayers(i).lowerE)
				*Db_h(EvanescentLayers(i).lowerE-1)/pow(w_0*dt_g,2)*AuxGrid_Eh(EvanescentLayers(i).lowerE-1));
			//solve for the E-field
			zgtsv_(&EvanescentLayers(i).N_e, &EvanescentLayers(i).nrhs, EvanescentLayers(i).lowdiagE_r.data(), EvanescentLayers(i).lowdiagE_i.data(), EvanescentLayers(i).diagE_r.data(), EvanescentLayers(i).diagE_i.data(), EvanescentLayers(i).updiagE_r.data(), EvanescentLayers(i).updiagE_i.data(), EvanescentLayers(i).righthandE_r.data(), EvanescentLayers(i).righthandE_i.data(), &EvanescentLayers(i).N_e, &info);
			if (info!=0)
			{
				cout << endl << "Matrix inversion failed!!" << endl;
				exit(-1);
			}

			//Place the computed E-field values into the auxiliary grid
			AuxGrid_Eh(Range(EvanescentLayers(i).lowerE,EvanescentLayers(i).upperE))
				= EvanescentLayers(i).righthandE_r + ii*EvanescentLayers(i).righthandE_i;
		}
		else
		{
			//Update Ee, Eh
			for (k=EvanescentLayers(i).lowerE-1; k<=EvanescentLayers(i).upperE+1; k++)
			{
				//He
				AuxGrid_Ee(k) = -AuxGrid_Ee(k) + 2.0*ii*Cb_e(k)/(w_0*dt_g)*(AuxGrid_He(k)-AuxGrid_He(k-1));
				//Hh
				AuxGrid_Eh(k) = -AuxGrid_Eh(k) + 2.0*ii*Cb_h(k)/(w_0*dt_g)*(AuxGrid_Hh(k)-AuxGrid_Hh(k-1));
			}
		}
	}

	//Update Ez
	for (k=AbsorbPointLower;k<=AbsorbPointUpper-1;k++)
	{
		if (IsEvanescent_z(k))
		{
			AuxGrid_Ez(k) = -AuxGrid_Ez(k) + 2.0*Ez_e(k)*AuxGrid_He(k);
		}
	}
}

void Cpw_ml::UpdateEvanescentH()
{
	for (i=0; i<NumOfEvanescentLayers; i++)
	{
		if (!EvanescentLayers(i).homogeneous_at_interfaces)
		{
			//Update He
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH; k++)
			{
				EvanescentLayers(i).diagH_r(k) = real(1.0-(Db_e(k)*(Cb_e(k+1)+Cb_e(k)))/pow(w_0*dt_g,2));
				EvanescentLayers(i).diagH_i(k) = imag(1.0-(Db_e(k)*(Cb_e(k+1)+Cb_e(k)))/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH-1; k++)
			{
				EvanescentLayers(i).updiagH_r(k) = real(Db_e(k)*Cb_e(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).updiagH_i(k) = imag(Db_e(k)*Cb_e(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagH_r(k+1) = real(Db_e(k+1)*Cb_e(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagH_i(k+1) = imag(Db_e(k+1)*Cb_e(k+1)/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH; k++)
			{
				EvanescentLayers(i).righthandH_r(k) = 0.0;
				EvanescentLayers(i).righthandH_i(k) = 0.0;
			}
			EvanescentLayers(i).righthandH_r(EvanescentLayers(i).upperH) = -real(Db_e(EvanescentLayers(i).upperH)
				*Cb_e(EvanescentLayers(i).upperH+1)/pow(w_0*dt_g,2)*AuxGrid_He(EvanescentLayers(i).upperH+1));
			EvanescentLayers(i).righthandH_i(EvanescentLayers(i).upperH) = -imag(Db_e(EvanescentLayers(i).upperH)
				*Cb_e(EvanescentLayers(i).upperH+1)/pow(w_0*dt_g,2)*AuxGrid_He(EvanescentLayers(i).upperH+1));
			EvanescentLayers(i).righthandH_r(EvanescentLayers(i).lowerH) = -real(Db_e(EvanescentLayers(i).lowerH)
				*Cb_e(EvanescentLayers(i).lowerH)/pow(w_0*dt_g,2)*AuxGrid_He(EvanescentLayers(i).lowerH-1));
			EvanescentLayers(i).righthandH_i(EvanescentLayers(i).lowerH) = -imag(Db_e(EvanescentLayers(i).lowerH)
				*Cb_e(EvanescentLayers(i).lowerH)/pow(w_0*dt_g,2)*AuxGrid_He(EvanescentLayers(i).lowerH-1));
			//solve for the H-field
			zgtsv_(&EvanescentLayers(i).N_h, &EvanescentLayers(i).nrhs, EvanescentLayers(i).lowdiagH_r.data(), EvanescentLayers(i).lowdiagH_i.data(), EvanescentLayers(i).diagH_r.data(), EvanescentLayers(i).diagH_i.data(), EvanescentLayers(i).updiagH_r.data(), EvanescentLayers(i).updiagH_i.data(), EvanescentLayers(i).righthandH_r.data(), EvanescentLayers(i).righthandH_i.data(), &EvanescentLayers(i).N_h, &info);
			if (info!=0)
			{
				cout << endl << "Matrix inversion failed!!" << endl;
				exit(-1);
			}

			//Place the computed H-field values into the auxiliary grid
			AuxGrid_He(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH))
				= EvanescentLayers(i).righthandH_r + ii*EvanescentLayers(i).righthandH_i;

			//Update Hh
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH; k++)
			{
				EvanescentLayers(i).diagH_r(k) = real(1.0-(Db_h(k)*(Cb_h(k+1)+Cb_h(k)))/pow(w_0*dt_g,2));
				EvanescentLayers(i).diagH_i(k) = imag(1.0-(Db_h(k)*(Cb_h(k+1)+Cb_h(k)))/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH-1; k++)
			{
				EvanescentLayers(i).updiagH_r(k) = real(Db_h(k)*Cb_h(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).updiagH_i(k) = imag(Db_h(k)*Cb_h(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagH_r(k+1) = real(Db_h(k+1)*Cb_h(k+1)/pow(w_0*dt_g,2));
				EvanescentLayers(i).lowdiagH_i(k+1) = imag(Db_h(k+1)*Cb_h(k+1)/pow(w_0*dt_g,2));
			}
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH; k++)
			{
				EvanescentLayers(i).righthandH_r(k) = 0.0;
				EvanescentLayers(i).righthandH_i(k) = 0.0;
			}
			EvanescentLayers(i).righthandH_r(EvanescentLayers(i).upperH) = -real(Db_h(EvanescentLayers(i).upperH)
				*Cb_h(EvanescentLayers(i).upperH+1)/pow(w_0*dt_g,2)*AuxGrid_Hh(EvanescentLayers(i).upperH+1));
			EvanescentLayers(i).righthandH_i(EvanescentLayers(i).upperH) = -imag(Db_h(EvanescentLayers(i).upperH)
				*Cb_h(EvanescentLayers(i).upperH+1)/pow(w_0*dt_g,2)*AuxGrid_Hh(EvanescentLayers(i).upperH+1));
			EvanescentLayers(i).righthandH_r(EvanescentLayers(i).lowerH) = -real(Db_h(EvanescentLayers(i).lowerH)
				*Cb_h(EvanescentLayers(i).lowerH)/pow(w_0*dt_g,2)*AuxGrid_Hh(EvanescentLayers(i).lowerH-1));
			EvanescentLayers(i).righthandH_i(EvanescentLayers(i).lowerH) = -imag(Db_h(EvanescentLayers(i).lowerH)
				*Cb_h(EvanescentLayers(i).lowerH)/pow(w_0*dt_g,2)*AuxGrid_Hh(EvanescentLayers(i).lowerH-1));
			//solve for the H-field
			zgtsv_(&EvanescentLayers(i).N_h, &EvanescentLayers(i).nrhs, EvanescentLayers(i).lowdiagH_r.data(), EvanescentLayers(i).lowdiagH_i.data(), EvanescentLayers(i).diagH_r.data(), EvanescentLayers(i).diagH_i.data(), EvanescentLayers(i).updiagH_r.data(), EvanescentLayers(i).updiagH_i.data(), EvanescentLayers(i).righthandH_r.data(), EvanescentLayers(i).righthandH_i.data(), &EvanescentLayers(i).N_h, &info);
			if (info!=0)
			{
				cout << endl << "Matrix inversion failed!!" << endl;
				exit(-1);
			}

			//Place the computed H-field values into the auxiliary grid
			AuxGrid_Hh(Range(EvanescentLayers(i).lowerH,EvanescentLayers(i).upperH))
				= EvanescentLayers(i).righthandH_r + ii*EvanescentLayers(i).righthandH_i;
		}
		else
		{
			//Update He, Hh
			for (k=EvanescentLayers(i).lowerH; k<=EvanescentLayers(i).upperH; k++)
			{
				//He
				AuxGrid_He(k) = -AuxGrid_He(k) + 2.0*ii*Db_e(k)/(w_0*dt_g)*(AuxGrid_Ee(k+1)-AuxGrid_Ee(k));
				//Hh
				AuxGrid_Hh(k) = -AuxGrid_Hh(k) + 2.0*ii*Db_h(k)/(w_0*dt_g)*(AuxGrid_Eh(k+1)-AuxGrid_Eh(k));
			}
		}
	}

	//Update Hz
	for (k=AbsorbPointLower+1;k<=AbsorbPointUpper-1;k++)
	{
		if (IsEvanescent_xy(k))
		{
			//Hz
			AuxGrid_Hz(k) = -AuxGrid_Hz(k) + 2.0*Hz_h(k)*AuxGrid_Eh(k);
		}
	}
}

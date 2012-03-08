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

//Declaration of the class "Cpw_2l" for a TF/SF plane-wave source in a lossless 2-layered medium

#include "headers.h"

#include "Cpw_2l.h"

//base Angora exception class
#include "angora_excp.h"

//for the definition of MaterialId
#include "material_id.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//Uses TinyVector operations
#include <blitz/tinyvec-et.h>

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern Array<bool,1> IsLayerGrounded;

extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern int number_of_layers;
extern Array<int,1> LayerLowerZIndices;

extern double epsilon_r_upper,epsilon_r_lower;
extern double c_upper,c_lower;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int GridIndex;
extern int rank;

extern int i,j,k;


//*************************************************//
//***************  Cpw_2l_hs class  ***************//
//*************************************************//

Cpw_2l_hs::Cpw_2l_hs(const PWDataType& MyData, const string& myhalfspace, const int& myHighestIndexInLowerHalfSpace, const double& epsilon_r_space)
//DON'T FORGET TO TURN OFF THE "RECORD AUXILIARY GRID" OPTION IN "DATA" BEFORE FEEDING IT INTO CPW_FS
	 	: Cpw_fs(MyData,epsilon_r_space), halfspace(myhalfspace), HighestIndexInLowerHalfSpace(myHighestIndexInLowerHalfSpace)
{
	 if ((halfspace!="upper")&&(halfspace!="lower"))
	 {
	 	ostringstream errstream;
	 	errstream << "Invalid half space (" << halfspace << ") for plane wave";
		throw AngoraDeveloperException(errstream.str());
	 }

	//resize and construct the boolean arrays that determine if the field component is in the desired half plane
	IsFullIntegerComponentInHalfSpace.resize(Range(1,NCELLS_Z+2*NPML+1));
	IsHalfIntegerComponentInHalfSpace.resize(Range(1,NCELLS_Z+2*NPML));

	for (int k=1; k<=NCELLS_Z+2*NPML+1; k++)
	{
		if (halfspace=="upper")
		{
			IsFullIntegerComponentInHalfSpace(k) = (k>=HighestIndexInLowerHalfSpace+1);	//is k above HighestIndexInLowerHalfSpace+1?
		}
		else if (halfspace=="lower")
		{
			IsFullIntegerComponentInHalfSpace(k) = (k<HighestIndexInLowerHalfSpace+1);		//is k below HighestIndexInLowerHalfSpace+1?
		}
	}
	for (int k=1; k<=NCELLS_Z+2*NPML; k++)
	{
		if (halfspace=="upper")
		{
			IsHalfIntegerComponentInHalfSpace(k) = (k>HighestIndexInLowerHalfSpace);	//is k above HighestIndexInLowerHalfSpace?
		}
		else if (halfspace=="lower")
		{
			IsHalfIntegerComponentInHalfSpace(k) = (k<=HighestIndexInLowerHalfSpace);	//is k below HighestIndexInLowerHalfSpace?
		}
	}
};


//**********************************************//
//***************  Cpw_2l class  ***************//
//**********************************************//

Cpw_2l::Cpw_2l(const PWDataType& MyData)
		: Cpw(MyData)
{
	if (number_of_layers!=2)
	{//if the # of layers is not 2, there is a developing error, so throw a developer exception
		throw AngoraDeveloperException("Cpw_2l class cannot be used if number of layers is not 2.");
	}

	//highest cell index in the lower half space is taken from the layering info
	HighestIndexInLowerHalfSpace = LayerLowerZIndices(1)-1;	//layer #1 is the upper half space

	//initialize with 3 plane waves: incident, reflected, transmitted
	PlaneWaves.resize(3);

	//plane-wave data objects for incident, reflected, transmitted waves
	PWDataType PWData_inc(Data),PWData_refl(Data),PWData_trans(Data);		//initialize all plane waves with the incident wave data

	//turn off some features in the child plane waves
// 	DisplayWarnings and DisplayCellsPerLambda may be turned off later

	//Determine the incidence layer
	if (CosT>=0)
	{//upper half space
		incidence_layer = "upper";
	}
	else
	{//lower half space
		incidence_layer = "lower";
	}

	//determine the limits of the three TF/SF boxes belonging to the incident, reflected and transmitted PWs
//	if (incidence_layer=="upper")
//	{//PW is incident from the upper half space
//		//remember that these are the *margins* from the PML interface, not the cell indices
//		PWData_inc.PWMarginLowerZ = HighestIndexInLowerHalfSpace-NPML-1;	//the lower face of the TF/SF box is 1 cell below the interface
//		PWData_refl.PWMarginLowerZ = HighestIndexInLowerHalfSpace-NPML-1; //the lower face of the TF/SF box is 1 cell below the interface
//		PWData_trans.PWMarginUpperZ = NCELLS_Z+NPML-HighestIndexInLowerHalfSpace-1; //the upper face of the TF/SF box is 1 cell above the interface
//	}
//	else if (incidence_layer=="lower")
//	{//PW is incident from the upper half space
//		//remember that these are the *margins* from the PML interface, not the cell indices
//		PWData_inc.PWMarginUpperZ = NCELLS_Z+NPML-HighestIndexInLowerHalfSpace-1;	//the upper face of the TF/SF box is 1 cell above the interface
//		PWData_refl.PWMarginUpperZ = NCELLS_Z+NPML-HighestIndexInLowerHalfSpace-1; //the upper face of the TF/SF box is 1 cell above the interface
//		PWData_trans.PWMarginLowerZ = HighestIndexInLowerHalfSpace-NPML-1; //the lower face of the TF/SF box is 1 cell below the interface
//	}

	//Normalized relative permittivities in the incidence and transmission layers
	eps_r_i = 1;
	if (incidence_layer=="upper")
	{//upper half space
		eps_r_t = epsilon_r_lower/epsilon_r_upper;
	}
	else if (incidence_layer=="lower")
	{//lower half space
		eps_r_t = epsilon_r_upper/epsilon_r_lower;
	}

	//Cosines of propagation angles
	if ((eps_r_i<pow(SinT,2))||(eps_r_t<pow(SinT,2)))
	{//is the wave evanescent in any of the layers?
		EvanescentPWException evanpw;
		throw evanpw; //throw the exception (caught from Ctfsf)
	}
	else
	{
		CosT_i = sqrt(1-pow2(SinT)/eps_r_i);
		CosT_t = sqrt(1-pow2(SinT)/eps_r_t);
	}

	//Fresnel transmission and reflection coefficients
	if (IsLayerGrounded(1)) //is the top layer grounded?
	{
		//Reflection coefficients are -1, transmission coefficients are 0
		Refl_TE = -1;
		Trans_TE = 0;
		Refl_TM = -1;
		Trans_TM = 0;
	}
	else
	{
		//TE transmission and reflection coefficients  (Balanis pg. 187)
		Refl_TE = ((CosT_i/sqrt(eps_r_t))-(CosT_t/sqrt(eps_r_i)))/((CosT_i/sqrt(eps_r_t))+(CosT_t/sqrt(eps_r_i)));
		Trans_TE = 1 + Refl_TE;
		//TM transmission and reflection coefficients  (Balanis pg. 191)
		Refl_TM = (-(CosT_i/sqrt(eps_r_i))+(CosT_t/sqrt(eps_r_t)))/((CosT_i/sqrt(eps_r_i))+(CosT_t/sqrt(eps_r_t)));
		Trans_TM = 2*(CosT_i/sqrt(eps_r_t))/((CosT_i/sqrt(eps_r_i))+(CosT_t/sqrt(eps_r_t)));
	}

	//The definitions of the Fresnel coefficients agree with those in Balanis,"Adv. Eng. EM", pg.185-191
	//First of all, remember that the polarization angle PSI of the E-field is measured clockwise w.r.t. (k_inc x ^z), with k_inc pointing away from the clock surface.
	//For TE incidence, Balanis' book assumes the direction of the reflected E-field to be the same as that of the incident E-field. Therefore, as long as (k_inc x ^z) stays the same upon reflection, the formulas agree.
	//For TM incidence, however, Balanis' book assumes the direction of the reflected H-field to be opposite of that of the incident H-field. Therefore, if (k_inc x ^z) stays the same, the TM reflection coefficient should be inverted in sign.
	//The direction convention for the transmitted fields always agree with the incident field; therefore, no change in the formulas is necessary.

	//reflection and transmission angles
	//reflection
	THETA_refl = M_PI - THETA_INC;
	//The sines of the incidence and reflection angles should have the same sign, so that (k_inc x ^z) stays the same upon reflection. Remember that PHI_refl=PHI, so it is the sign of sin(THETA_refl) that determines the direction of k_inc_lateral (see Cpw.cpp), therefore the sign of (k_inc x ^z).
	//If SinT is not very close to 0, the signs of THETA_INC and THETA_refl should normally have the same sign. However, if this is not the case, THETA_refl might end up very close to the z-axis, but on the wrong side of it. The following check corrects this exceptional situation.
	if (((SinT>=0)&&(sin(THETA_refl)<0))||((SinT<0)&&(sin(THETA_refl)>=0)))	//are the signs of the incidence and reflection angles the same?
	{
		THETA_refl = -THETA_refl;
	}
	PHI_refl = PHI_INC;		//same principal plane, same PHI
	//transmission
	THETA_trans = asin(sin(THETA_INC)/sqrt(eps_r_t/eps_r_i));
	if (incidence_layer=="lower")
	{//result of "asin" is between -pi and pi, so convert to lower-half space incidence direction
		THETA_trans = M_PI-THETA_trans;
	}
	PHI_trans = PHI_INC;		//same principal plane, same PHI

	//reflected and transmitted electric field amplitudes
	E_TE = Data.E0*CosPsi;
	E_TM = Data.E0*SinPsi;
	E_TE_refl = Refl_TE*E_TE;
	E_TM_refl = -Refl_TM*E_TM;	///Note the (-) sign: See the explanation above!
	E_TE_trans = Trans_TE*E_TE;
	E_TM_trans = Trans_TM*E_TM;

	E_refl = sqrt(pow2(E_TE_refl)+pow2(E_TM_refl));
	E_trans = sqrt(pow2(E_TE_trans)+pow2(E_TM_trans));
	//polarization angles of reflected and transmitted waves
	PSI_refl = atan2(E_TM_refl,E_TE_refl);
	PSI_trans = atan2(E_TM_trans,E_TE_trans);

	//determine the plane-wave origins of the reflected and scattered waves
	//This is an issue only when PWOriginZ!=OriginZ
	//The x and y coordinates of the PW origins of the reflected and transmitted waves are always the same as the incident wave.
	//However, the z coordinates must be adjusted to match the waves at the interface.
	// Origins of the incident and reflected waves are symmetrical w.r.t. the interface
	PWData_refl.PWOriginZ = (HighestIndexInLowerHalfSpace+1)-(Data.PWOriginZ-(HighestIndexInLowerHalfSpace+1));
	// Origin of the transmitted wave
	PWData_trans.PWOriginZ = (HighestIndexInLowerHalfSpace+1)-((HighestIndexInLowerHalfSpace+1)-Data.PWOriginZ)*sqrt(eps_r_i-pow2(SinT))/sqrt(eps_r_t-pow2(SinT));

	//INCIDENT WAVE:
	//data object for the incident wave stays the same as the original data object (since the original data object specifies the incident plane wave in the first place)
	//if the incident wave is incident from the upper hemisphere:
	if (incidence_layer=="upper")
	{//incident wave is in the upper half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_inc,"upper",HighestIndexInLowerHalfSpace,epsilon_r_upper));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[0] = new_pw_ptr;
	}
	//if the incident wave is incident from the lower hemisphere:
	else if (incidence_layer=="lower")
	{//incident wave is in the lower half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_inc,"lower",HighestIndexInLowerHalfSpace,epsilon_r_lower));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[0] = new_pw_ptr;
	}

	//REFLECTED WAVE:
	PWData_refl.E0 = E_refl;
	PWData_refl.THETA = THETA_refl;
	PWData_refl.PHI = PHI_refl;
	PWData_refl.PSI = PSI_refl;
	//if the incident wave is incident from the upper hemisphere:
	if (incidence_layer=="upper")
	{//reflected wave is in the upper half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_refl,"upper",HighestIndexInLowerHalfSpace,epsilon_r_upper));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[1] = new_pw_ptr;
	}
	//if the incident wave is incident from the lower hemisphere:
	else if (incidence_layer=="lower")
	{//reflected wave is in the lower half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_refl,"lower",HighestIndexInLowerHalfSpace,epsilon_r_lower));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[1] = new_pw_ptr;
	}

	//TRANSMITTED WAVE:
	PWData_trans.E0 = E_trans;
	PWData_trans.THETA = THETA_trans;
	PWData_trans.PHI = PHI_trans;
	PWData_trans.PSI = PSI_trans;
	//if the incident wave is incident from the upper hemisphere:
	if (incidence_layer=="upper")
	{//transmitted wave is in the lower half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_trans,"lower",HighestIndexInLowerHalfSpace,epsilon_r_lower));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[2] = new_pw_ptr;
	}
	//if the incident wave is incident from the lower hemisphere:
	else if (incidence_layer=="lower")
	{//transmitted wave is in the lower half space
		boost::shared_ptr<Cpw_2l_hs> new_pw_ptr(new Cpw_2l_hs(PWData_trans,"upper",HighestIndexInLowerHalfSpace,epsilon_r_upper));//all PWs have the same PW index (that of the Cpw_2l object)
		PlaneWaves[2] = new_pw_ptr;
	}
}

void Cpw_2l::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the incident, reflected and transmitted plane waves.
	for (int planewaveindex=0; planewaveindex<3; planewaveindex++)		//there are 3 plane waves (incident, reflected, transmitted)
	{
		PlaneWaves[planewaveindex]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Cpw_2l::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the incident, reflected and transmitted plane waves.
	for (int planewaveindex=0; planewaveindex<3; planewaveindex++)		//there are 3 plane waves (incident, reflected, transmitted)
	{
		PlaneWaves[planewaveindex]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

//**********************************************************************************************************************//
//These are defined (since they are pure virtual in Cpw), but left blank                                                //
// What is done by these functions is instead done by CorrectE, CorrectH (virtual Cpw functions, redefined in Cpw_2l)   //
//**********************************************************************************************************************//
void Cpw_2l::UpdateIncidentH(const int& n)
{
}
void Cpw_2l::UpdateIncidentE(const int& n)
{
}
void Cpw_2l::ApplyCorrectionE(const int& n)
{
}
void Cpw_2l::ApplyCorrectionH(const int& n)
{
}
//**********************************************************************************************************************//
//**********************************************************************************************************************//
//**********************************************************************************************************************//


void Cpw_2l::WriteScatteredPWDirection(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angle (THETA and PHI) of the scattered PW into PW_THETA and PW_PHI
	//NOTE: "WriteScatteredPWDirection(*,*) is inherited by Cpw_2l_hs from Cpw_fs
	//reflected plane wave has index 1
	PlaneWaves[1]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	//transmitted plane wave has index 2
	PlaneWaves[2]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
}

void Cpw_2l::WriteScatteredPWDelayFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delay (from the origin) of the scattered PW into origindelay_array
	//NOTE: "WriteScatteredPWDelayFromOrigin(*,*) is inherited by Cpw_2l_hs from Cpw_fs
	//reflected plane wave has index 1
	PlaneWaves[1]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	//transmitted plane wave has index 2
	PlaneWaves[2]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
}

void Cpw_2l::WriteScatteredPWFieldAmplitude(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitude of the scattered PWs into field_array
	//NOTE: "WriteScatteredPWFieldAmplitude(*,*) is inherited by Cpw_2l_hs from Cpw_fs
	//reflected plane wave has index 1
	PlaneWaves[1]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	//transmitted plane wave has index 2
	PlaneWaves[2]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
}

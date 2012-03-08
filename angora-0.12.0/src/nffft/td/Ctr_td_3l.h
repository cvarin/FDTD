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

#ifndef CTR_TD_3L_H
#define CTR_TD_3L_H

//Declaration of the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_3l" for 3-layered lossless media
//Derived from the abstract base class "Ctr_td".

#include "Ctr_td.h"


class Ctr_td_3l : public Ctr_td
{
 public:
	 Ctr_td_3l(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index);		//constructor

	 void UpdateFarField(const int& n);		//Update auxiliary far-field waveforms using the current E,H on the virtual surface
	 void ConstructFarField();	//Construct far-field waveforms using the final auxiliary waveforms

 private:
	 void UpdateElectric(const int& n);		//Update Electric_{X,Y,Z}_{e,h}, using the H values on the virtual surface
	 void UpdateMagnetic(const int& n);		//Update Magnetic_{X,Y,Z}_{e,h}, using the E values on the virtual surface

	 void Update_Electric_X(const int& i, const int& j, const int& k, const int& n);	//Update Electric_X_{e,h}
	 void Update_Electric_Y(const int& i, const int& j, const int& k, const int& n);	//Update Electric_Y_{e,h}
	 void Update_Electric_Z(const int& i, const int& j, const int& k, const int& n);	//Update Electric_Z_{e,h}
	 void Update_Magnetic_X(const int& i, const int& j, const int& k, const int& n);	//Update Magnetic_X_{e,h}
	 void Update_Magnetic_Y(const int& i, const int& j, const int& k, const int& n);	//Update Magnetic_Y_{e,h}
	 void Update_Magnetic_Z(const int& i, const int& j, const int& k, const int& n);	//Update Magnetic_Z_{e,h}

	 double TheoreticalFarFieldTheta(const int& n);		//the theoretical theta-component of the far field
	 double TheoreticalFarFieldPhi(const int& n);		//the theoretical phi-component of the far field

	 //Index of the highest cell in dielectric slab
	 int SlabPos;
	 //Thickness of dielectric slab
	 int SlabThickness;
	 //Dielectric permittivity of the slab
	 double epsilon_r_slab;
	 //is the dielectric slab grounded?
	 bool grounded;

// 	 double TheoreticalFarFieldWaveform(double t, double tau);	//Return the appropriate Gaussian waveform, depending on the source current waveform

	 void Calculate_TL_GreenFunctions();	//calculates the TL Green functions
	 void Calculate_TL_GreenFunctions_UpperHalfSpace();	//calculates the TL Green functions for the upper half space
	 void Calculate_TL_GreenFunctions_LowerHalfSpace();	//calculates the TL Green functions for the lower half space

	 double epsilon_o;	//relative permittivity in the observation half space s(either uppermost or lowermost halfspace)
	 double c_o;	//velocity of light in the observation half space
	 double Z_o;	//wave impedance in the observation half space

	 double eps_r_0,eps_r_1,eps_r_2;	//relative permittivities of the uppermost, middle, and lowermost layers, normalized by epsilon_o
	 double vp_0,vp_1,vp_2;	//transmission-line velocities in the uppermost, middle, and lowermost layers

	 double Z_e_0,Z_e_1,Z_e_2;	//TM transmission-line impedances in the uppermost, middle, and lowermost layers
	 double Z_h_0,Z_h_1,Z_h_2;	//TE transmission-line impedances in the uppermost, middle, and lowermost layers

	 double Refl_e_10,Trans_e_10;	//TM refl./trans. coefficients passing from the middle to uppermost layer
	 double Refl_e_01,Trans_e_01;	//TM refl./trans. coefficients passing from the uppermost to middle layer
	 double Refl_e_21,Trans_e_21;	//TM refl./trans. coefficients passing from the lowermost to middle layer
	 double Refl_e_12,Trans_e_12;	//TM refl./trans. coefficients passing from the middle to lowermost layer

	 double Refl_h_10,Trans_h_10;	//TE refl./trans. coefficients passing from the middle to uppermost layer
	 double Refl_h_01,Trans_h_01;	//TE refl./trans. coefficients passing from the uppermost to middle layer
	 double Refl_h_21,Trans_h_21;	//TE refl./trans. coefficients passing from the lowermost to middle layer
	 double Refl_h_12,Trans_h_12;	//TE refl./trans. coefficients passing from the middle to lowermost layer

	 struct Impulse
	//Represents an impulse train, each element being the amplitude and delay of one element. Elements do not have to be ordered in increasing delays.
	 {
		 double Amp;
		 double Delay;
	 };

	 typedef Array<Impulse,1> ImpulseTrain;	//impulse train is defined as an array of impulses

	 ImpulseTrain V_e;	//TM voltage impulse train at z=0+ due to a unit impulsive voltage wave at z=0- traveling in the z+ direction
	 ImpulseTrain V_h;	//TE voltage impulse train at z=0+ due to a unit impulsive voltage wave at z=0- traveling in the z+ direction

	 //Transmission-line (TL) Green functions
	 Array<ImpulseTrain,1> V_e_v;		//V^(e)_(v) at different vertical locations
	 Array<ImpulseTrain,1> V_h_v;		//V^(e)_(v) at different vertical locations
	 Array<ImpulseTrain,1> V_e_i;		//V^(e)_(v) at different vertical locations
	 Array<ImpulseTrain,1> V_h_i;		//V^(e)_(v) at different vertical locations

	 //Extra steps required for the far-field arrays
	 double MaxLateralAdvance,MaxVertAdvance,MaxTotalAdvance;	//Maximum time advances in different directions, and maximum total advance
	 double MaxLateralDelay,MaxVertDelay,MaxTotalDelay;	//Maximum delays in different directions, and maximum total delay
	 double MaxDelay;	//the maximum delay due to reflections
	 int ExtraSteps;
	 int TotalSteps;				//Length of the far-field waveforms
	 double r_offset;

	 //Far-field arrays
	 Array<double,1> Electric_X_e;	//TM contribution of the x-component of the electric current (J) to the far field (differentiation, integration and constant factor omitted, done later in post-processing--ConstructFarField)
	 Array<double,1> Electric_X_h;	//TE contribution of the x-component of the electric current (J) to the far field(differentiation, integration and constant factor omitted, done later in post-processing--ConstructFarField)
	 Array<double,1> Electric_Y_e;	//TM contribution of the y-component of the electric current (J) to the far field(differentiation, integration and constant factor omitted, done later in post-processing--ConstructFarField)
	 Array<double,1> Electric_Y_h;	//TE contribution of the y-component of the electric current (J) to the far field(differentiation, integration and constant factor omitted, done later in post-processing--ConstructFarField)
	 Array<double,1> Electric_Z_e;	//TM contribution of the z-component of the electric current (J) to the far field(differentiation, integration and constant factor omitted, done later in post-processing--ConstructFarField)
	 //(No TE contribution from Jz)

	 Array<double,1> Magnetic_X_e;	//TM contribution of the x-component of the magnetic current (M) to the far field
	 Array<double,1> Magnetic_X_h;	//TE contribution of the x-component of the magnetic current (M) to the far field
	 Array<double,1> Magnetic_Y_e;	//TM contribution of the y-component of the magnetic current (M) to the far field
	 Array<double,1> Magnetic_Y_h;	//TE contribution of the y-component of the magnetic current (M) to the far field
	 Array<double,1> Magnetic_Z_h;	//TE contribution of the z-component of the magnetic current (M) to the far field
	 //(No TM contribution from Mz)

	 //Magnetic_(X,Y,Z)_(e,h) lags Electric_(X,Y,Z)_(e,h) by "0.5dt" (since E lags H by the same amount after an update)

	 //Post-processing variables
	 Array<double,1> Theta_J;	//Contribution of Js to the Theta component
	 Array<double,1> Phi_J;		//Contribution of Js to the Phi component
	 Array<double,1> Phi_M;		//Contribution of Ms to the Theta component
	 Array<double,1> Theta_M;	//Contribution of Ms to the Phi component

	 double delay_total,delay_lateral,delay_vertical,Jx_e,Jx_h,Jy_e,Jy_h,Jz_e,Mx_e,Mx_h,My_e,My_h,Mz_h;
	 int impulse;	//impulse index (for looping through the impulses in an impulse train)
};

#endif

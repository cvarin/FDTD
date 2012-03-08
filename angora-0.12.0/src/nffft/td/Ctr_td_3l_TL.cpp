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

//Defines the transmission-line Green functions in the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_3l" for 3-layered lossless media

#include "headers.h"

#include "Ctr_td_3l.h"

extern double dx,dt;
extern int NCELLS_Z,NPML;

extern double c_o;


void Ctr_td_3l::Calculate_TL_GreenFunctions()
{//calculates the TL Green functions

	//transmission-line velocities in different layers
	vp_0 = c_o/sqrt(eps_r_0-pow(SinT,2));
	vp_1 = c_o/sqrt(eps_r_1-pow(SinT,2));
	vp_2 = c_o/sqrt(eps_r_2-pow(SinT,2));

	//Transmission-line impedances
	//TM impedances
	Z_e_0 = Z_o*sqrt(eps_r_0-pow(SinT,2))/eps_r_0;
	Z_e_1 = Z_o*sqrt(eps_r_1-pow(SinT,2))/eps_r_1;
	Z_e_2 = Z_o*sqrt(eps_r_2-pow(SinT,2))/eps_r_2;
	//TE impedances
	Z_h_0 = Z_o/sqrt(eps_r_0-pow(SinT,2));
	Z_h_1 = Z_o/sqrt(eps_r_1-pow(SinT,2));
	Z_h_2 = Z_o/sqrt(eps_r_2-pow(SinT,2));

	//Transmission and reflection coefficients
	//TM transmission and reflection coefficients
	Refl_e_10 = ((Z_e_0-Z_e_1)/(Z_e_0+Z_e_1));
	Trans_e_10 = 1 + Refl_e_10;
	Refl_e_01 = -Refl_e_10;
	Trans_e_01 = 1 + Refl_e_01;
	if (grounded)
	{
		Refl_e_21 = -1;
		Trans_e_21 = 0;
		Refl_e_12 = -1;
		Trans_e_12 = 0;
	}
	else
	{
		Refl_e_21 = ((Z_e_1-Z_e_2)/(Z_e_1+Z_e_2));
		Trans_e_21 = 1 + Refl_e_21;
		Refl_e_12 = -Refl_e_21;
		Trans_e_12 = 1 + Refl_e_12;
	}

	//TE transmission and reflection coefficients
	Refl_h_10 = ((Z_h_0-Z_h_1)/(Z_h_0+Z_h_1));
	Trans_h_10 = 1 + Refl_h_10;
	Refl_h_01 = -Refl_h_10;
	Trans_h_01 = 1 + Refl_h_01;
	if (grounded)
	{
		Refl_h_21 = -1;
		Trans_h_21 = 0;
		Refl_h_12 = -1;
		Trans_h_12 = 0;
	}
	else
	{
		Refl_h_21 = ((Z_h_1-Z_h_2)/(Z_h_1+Z_h_2));
		Trans_h_21 = 1 + Refl_h_21;
		Refl_h_12 = -Refl_h_21;
		Trans_h_12 = 1 + Refl_h_12;
	}

	if (CosT>=0)
	{
		Calculate_TL_GreenFunctions_UpperHalfSpace();
	}
	else
	{
		Calculate_TL_GreenFunctions_LowerHalfSpace();
	}
}

void Ctr_td_3l::Calculate_TL_GreenFunctions_UpperHalfSpace()
{//calculates the TL Green functions for the upper half space
	//V_e and V_h
	V_e.resize(1);
	V_e(0).Amp = Trans_e_10;
	V_e(0).Delay = 0;
	double refl = Refl_e_10*Refl_e_12;
	double DecayThreshold = 1e-2;
	while (abs(refl)>DecayThreshold)
	{
		V_e.resizeAndPreserve(V_e.size()+1);
		V_e(V_e.size()-1).Amp = refl*Trans_e_10;
		V_e(V_e.size()-1).Delay = V_e(V_e.size()-2).Delay + 2*SlabThickness*dx/vp_1/dt;
		refl *= Refl_e_10*Refl_e_12;
	}

	V_h.resize(1);
	V_h(0).Amp = Trans_h_10;
	V_h(0).Delay = 0;
	refl = Refl_h_10*Refl_h_12;
	while (abs(refl)>DecayThreshold)
	{
		V_h.resizeAndPreserve(V_h.size()+1);
		V_h(V_h.size()-1).Amp = refl*Trans_h_10;
		V_h(V_h.size()-1).Delay = V_h(V_h.size()-2).Delay + 2*SlabThickness*dx/vp_1/dt;
		refl *= Refl_h_10*Refl_h_12;
	}

	//Define V_e_v,V_h_v
	//V_e_v,V_h_v is defined at half-integer positions
	V_e_v.resize(Range(1,NCELLS_Z+2*NPML));
	V_h_v.resize(Range(1,NCELLS_Z+2*NPML));
	//below the slab
	for (int k=1; k<=SlabPos-SlabThickness; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_2/dt;
		//V_e_v
		V_e_v(k).resize(V_e.size());
		//only upward-traveling impulse contributes
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse).Amp = Trans_e_21*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse).Delay = V_e(impulse).Delay - delay_vertical
				- dx*SlabThickness/vp_2/dt
				+ dx*SlabThickness/vp_1/dt;
		}
		//V_h_v
		V_h_v(k).resize(V_h.size());
		//only upward-traveling impulse contributes
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse).Amp = Trans_h_21*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse).Delay = V_h(impulse).Delay - delay_vertical
				- dx*SlabThickness/vp_2/dt
				+ dx*SlabThickness/vp_1/dt;
		}
	}
	//inside the slab
	for (int k=SlabPos-SlabThickness+1; k<=SlabPos; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_1/dt;
		//V_e_v
		V_e_v(k).resize(2*V_e.size());
		//upward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse).Amp = 0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse).Delay = V_e(impulse).Delay - delay_vertical;
		}
		//downward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse+V_e.size()).Amp = -Refl_e_12*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse+V_e.size()).Delay = V_e(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
		//V_h_v
		V_h_v(k).resize(2*V_h.size());
		//upward-traveling impulse
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse).Amp = 0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse).Delay = V_h(impulse).Delay - delay_vertical;
		}
		//downward-traveling impulse
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse+V_h.size()).Amp = -Refl_h_12*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse+V_h.size()).Delay = V_h(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
	}
	//above the slab
	for (int k=SlabPos+1; k<=NCELLS_Z+2*NPML; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_0/dt;
		V_e_v(k).resize(V_e.size()+2);
		V_h_v(k).resize(V_h.size()+2);
		//upward-traveling impulse
		V_e_v(k)(0).Amp = 0.5;
		V_e_v(k)(0).Delay = -delay_vertical;
		V_h_v(k)(0).Amp = 0.5;
		V_h_v(k)(0).Delay = -delay_vertical;
		//downward-traveling impulse
		V_e_v(k)(1).Amp = -0.5*Refl_e_01;
		V_e_v(k)(1).Delay = delay_vertical;
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse+2).Amp = -Refl_e_12*Trans_e_01*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse+2).Delay = V_e(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
		V_h_v(k)(1).Amp = -0.5*Refl_h_01;
		V_h_v(k)(1).Delay = delay_vertical;
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse+2).Amp = -Refl_h_12*Trans_h_01*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse+2).Delay = V_h(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
	}

	//Define V_e_i,V_h_i
	//V_e_i,V_h_i is defined at full-integer positions
	V_e_i.resize(Range(1,NCELLS_Z+2*NPML+1));
	V_h_i.resize(Range(1,NCELLS_Z+2*NPML+1));
	//below the slab and at the lower interface
	//(it can be shown that at the interface, inside and outside Green functions are the same)
	for (int k=1; k<=SlabPos-SlabThickness+1; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_2/dt;
		//V_e_i
		V_e_i(k).resize(V_e.size());
		//only upward-traveling impulse contributes
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse).Amp = Z_e_2*Trans_e_21*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse).Delay = V_e(impulse).Delay - delay_vertical
				- dx*SlabThickness/vp_2/dt
				+ dx*SlabThickness/vp_1/dt;
		}
		//V_h_i
		V_h_i(k).resize(V_h.size());
		//only upward-traveling impulse contributes
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse).Amp = Z_h_2*Trans_h_21*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse).Delay = V_h(impulse).Delay - delay_vertical
				- dx*SlabThickness/vp_2/dt
				+ dx*SlabThickness/vp_1/dt;
		}
	}
	//inside the slab
	for (int k=SlabPos-SlabThickness+2; k<=SlabPos; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_1/dt;
		//V_e_i
		V_e_i(k).resize(2*V_e.size());
		//upward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse).Amp = Z_e_1*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse).Delay = V_e(impulse).Delay - delay_vertical;
		}
		//downward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse+V_e.size()).Amp = Z_e_1*Refl_e_12*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse+V_e.size()).Delay = V_e(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
		//V_h_i
		V_h_i(k).resize(2*V_h.size());
		//upward-traveling impulse
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse).Amp = Z_h_1*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse).Delay = V_h(impulse).Delay - delay_vertical;
		}
		//downward-traveling impulse
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse+V_h.size()).Amp = Z_h_1*Refl_h_12*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse+V_h.size()).Delay = V_h(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
	}
	//above the slab and at the upper interface
	//(it can be shown that at the interface, inside and outside Green functions are the same)
	for (int k=SlabPos+1; k<=NCELLS_Z+2*NPML+1; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_0/dt;
		V_e_i(k).resize(V_e.size()+2);
		V_h_i(k).resize(V_h.size()+2);
		//upward-traveling impulse
		V_e_i(k)(0).Amp = Z_e_0*0.5;
		V_h_i(k)(0).Amp = Z_h_0*0.5;
		V_e_i(k)(0).Delay = -delay_vertical;
		V_h_i(k)(0).Delay = -delay_vertical;
		//downward-traveling impulse
		V_e_i(k)(1).Amp = Z_e_0*0.5*Refl_e_01;
		V_h_i(k)(1).Amp = Z_h_0*0.5*Refl_h_01;
		V_e_i(k)(1).Delay = delay_vertical;
		V_h_i(k)(1).Delay = delay_vertical;
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse+2).Amp = Z_e_0*Refl_e_12*0.5*Trans_e_01*V_e(impulse).Amp;
			V_e_i(k)(impulse+2).Delay = V_e(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse+2).Amp = Z_h_0*Refl_h_12*0.5*Trans_h_01*V_h(impulse).Amp;
			V_h_i(k)(impulse+2).Delay = V_h(impulse).Delay + delay_vertical
															+ 2*dx*SlabThickness/vp_1/dt;
		}
	}
}

void Ctr_td_3l::Calculate_TL_GreenFunctions_LowerHalfSpace()
{//calculates the TL Green functions for the lower half space
	//V_e and V_h
	V_e.resize(1);
	V_e(0).Amp = Trans_e_12;
	V_e(0).Delay = 0;
	double refl = Refl_e_12*Refl_e_10;
	double DecayThreshold = 1e-2;
	while (abs(refl)>DecayThreshold)
	{
		V_e.resizeAndPreserve(V_e.size()+1);
		V_e(V_e.size()-1).Amp = refl*Trans_e_12;
		V_e(V_e.size()-1).Delay = V_e(V_e.size()-2).Delay + 2*SlabThickness*dx/vp_1/dt;
		refl *= Refl_e_12*Refl_e_10;
	}

	V_h.resize(1);
	V_h(0).Amp = Trans_h_12;
	V_h(0).Delay = 0;
	refl = Refl_h_12*Refl_h_10;
	while (abs(refl)>DecayThreshold)
	{
		V_h.resizeAndPreserve(V_h.size()+1);
		V_h(V_h.size()-1).Amp = refl*Trans_h_12;
		V_h(V_h.size()-1).Delay = V_h(V_h.size()-2).Delay + 2*SlabThickness*dx/vp_1/dt;
		refl *= Refl_h_12*Refl_h_10;
	}

	//Define V_e_v,V_h_v
	//V_e_v,V_h_v is defined at half-integer positions
	V_e_v.resize(Range(1,NCELLS_Z+2*NPML));
	V_h_v.resize(Range(1,NCELLS_Z+2*NPML));
	//below the slab
	for (int k=1; k<=SlabPos-SlabThickness; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_2/dt;
		V_e_v(k).resize(V_e.size()+2);
		V_h_v(k).resize(V_h.size()+2);
		//downward-traveling impulse
		V_e_v(k)(0).Amp = -0.5;
		V_e_v(k)(0).Delay = delay_vertical;
		V_h_v(k)(0).Amp = -0.5;
		V_h_v(k)(0).Delay = delay_vertical;
		//upward-traveling impulse
		V_e_v(k)(1).Amp = 0.5*Refl_e_21;
		V_h_v(k)(1).Amp = 0.5*Refl_h_21;
		V_e_v(k)(1).Delay = -(delay_vertical + 2*dx*SlabThickness/vp_2/dt);
		V_h_v(k)(1).Delay = -(delay_vertical + 2*dx*SlabThickness/vp_2/dt);
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse+2).Amp = Refl_e_10*Trans_e_21*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse+2).Delay = V_e(impulse).Delay
										- (delay_vertical + 2*dx*SlabThickness/vp_2/dt)
										+ 2*dx*SlabThickness/vp_1/dt;
		}
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse+2).Amp = Refl_h_10*Trans_h_21*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse+2).Delay = V_h(impulse).Delay
										- (delay_vertical + 2*dx*SlabThickness/vp_2/dt)
										+ 2*dx*SlabThickness/vp_1/dt;
		}
	}
	//inside the slab
	for (int k=SlabPos-SlabThickness+1; k<=SlabPos; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_1/dt;
		//V_e_v
		V_e_v(k).resize(2*V_e.size());
		//downward-traveling impulses
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse).Amp = -0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse).Delay = V_e(impulse).Delay
									+ (delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//upward-traveling impulses
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse+V_e.size()).Amp = Refl_e_10*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse+V_e.size()).Delay = V_e(impulse).Delay
									+ (-delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//V_h_v
		V_h_v(k).resize(2*V_h.size());
		//downward-traveling impulses
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse).Amp = -0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse).Delay = V_h(impulse).Delay
									+ (delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//upward-traveling impulses
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse+V_h.size()).Amp = Refl_h_10*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse+V_h.size()).Delay = V_h(impulse).Delay
									+ (-delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
	}
	//above the slab
	for (int k=SlabPos+1; k<=NCELLS_Z+2*NPML; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ+0.5)/vp_0/dt;
		//V_e_v
		V_e_v(k).resize(V_e.size());
		//only downward-traveling impulse contributes
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_v(k)(impulse).Amp = -Trans_e_01*0.5*V_e(impulse).Amp;
			V_e_v(k)(impulse).Delay = V_e(impulse).Delay
									+ delay_vertical
									+ dx*SlabThickness/vp_1/dt
									- dx*SlabThickness/vp_2/dt;
		}
		//V_h_v
		V_h_v(k).resize(V_h.size());
		//only downward-traveling impulse contributes
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_v(k)(impulse).Amp = -Trans_h_01*0.5*V_h(impulse).Amp;
			V_h_v(k)(impulse).Delay = V_h(impulse).Delay
									+ delay_vertical
									+ dx*SlabThickness/vp_1/dt
									- dx*SlabThickness/vp_2/dt;
		}
	}

	//Define V_e_i,V_h_i
	//V_e_i,V_h_i is defined at full-integer positions
	V_e_i.resize(Range(1,NCELLS_Z+2*NPML+1));
	V_h_i.resize(Range(1,NCELLS_Z+2*NPML+1));
	//below the slab and at the lower interface
	//(it can be shown that at the interface, inside and outside Green functions are the same)
	for (int k=1; k<=SlabPos-SlabThickness+1; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_2/dt;
		V_e_i(k).resize(V_e.size()+2);
		V_h_i(k).resize(V_h.size()+2);
		//upward-traveling impulse
		V_e_i(k)(0).Amp = Z_e_2*0.5;
		V_h_i(k)(0).Amp = Z_h_2*0.5;
		V_e_i(k)(0).Delay = delay_vertical;
		V_h_i(k)(0).Delay = delay_vertical;
		//downward-traveling impulse
		V_e_i(k)(1).Amp = Z_e_2*0.5*Refl_e_21;
		V_h_i(k)(1).Amp = Z_h_2*0.5*Refl_h_21;
		V_e_i(k)(1).Delay = -(delay_vertical + 2*dx*SlabThickness/vp_2/dt);
		V_h_i(k)(1).Delay = -(delay_vertical + 2*dx*SlabThickness/vp_2/dt);
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse+2).Amp = Z_e_2*Refl_e_10*Trans_e_21*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse+2).Delay = V_e(impulse).Delay
										- (delay_vertical + 2*dx*SlabThickness/vp_2/dt)
														  + 2*dx*SlabThickness/vp_1/dt;
		}
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse+2).Amp = Z_h_2*Refl_h_10*Trans_h_21*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse+2).Delay = V_h(impulse).Delay
										- (delay_vertical + 2*dx*SlabThickness/vp_2/dt)
														  + 2*dx*SlabThickness/vp_1/dt;
		}
	}
	//inside the slab
	for (int k=SlabPos-SlabThickness+2; k<=SlabPos; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_1/dt;
		//V_e_i
		V_e_i(k).resize(2*V_e.size());
		//downward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse).Amp = Z_e_1*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse).Delay = V_e(impulse).Delay
									+ (delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//upward-traveling impulse
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse+V_e.size()).Amp = Z_e_1*Refl_e_10*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse+V_e.size()).Delay = V_e(impulse).Delay
									+ (-delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//V_h_i
		V_h_i(k).resize(2*V_h.size());
		//downward-traveling impulses
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse).Amp = Z_h_1*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse).Delay = V_h(impulse).Delay
									+ (delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
		//upward-traveling impulse
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse+V_h.size()).Amp = Z_h_1*Refl_h_10*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse+V_h.size()).Delay = V_h(impulse).Delay
									+ (-delay_vertical + dx*SlabThickness/vp_1/dt)
									-  dx*SlabThickness/vp_2/dt;
		}
	}
	//above the slab and at the upper interface
	//(it can be shown that at the interface, inside and outside Green functions are the same)
	for (int k=SlabPos+1; k<=NCELLS_Z+2*NPML+1; k++)
	{
		delay_vertical = dx*(k-FarFieldOriginZ)/vp_0/dt;
		//V_e_i
		V_e_i(k).resize(V_e.size());
		//only downward-traveling impulse contributes
		for (int impulse=0; impulse<V_e.size(); impulse++)
		{
			V_e_i(k)(impulse).Amp = Z_e_0*Trans_e_01*0.5*V_e(impulse).Amp;
			V_e_i(k)(impulse).Delay = V_e(impulse).Delay
									+ delay_vertical
									+ dx*SlabThickness/vp_1/dt
									- dx*SlabThickness/vp_2/dt;
		}
		//V_h_i
		V_h_i(k).resize(V_h.size());
		//only downward-traveling impulse contributes
		for (int impulse=0; impulse<V_h.size(); impulse++)
		{
			V_h_i(k)(impulse).Amp = Z_h_0*Trans_h_01*0.5*V_h(impulse).Amp;
			V_h_i(k)(impulse).Delay = V_h(impulse).Delay
									+ delay_vertical
									+ dx*SlabThickness/vp_1/dt
									- dx*SlabThickness/vp_2/dt;
		}
	}
}

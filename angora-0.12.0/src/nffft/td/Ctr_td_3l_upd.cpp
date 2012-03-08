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

// Defines the updates due to magnetic (M) and electric (J) surface currents in the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_3l" for 3-layered lossless media

#include "headers.h"

#include "Ctr_td_3l.h"

//for the definition of MaterialId
#include "material_id.h"

extern double dx,dt;

extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
extern Array<double,1> eps_z;

extern double c_o;

extern int rank;

void Ctr_td_3l::Update_Magnetic_X(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to Mx
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX)*SinTCosP+(j-FarFieldOriginY+0.5)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_e_v:
	for (impulse=0;impulse<V_e_v(k).size();impulse++)
	{
		delay_total = delay_lateral + V_e_v(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Mx_e = V_e_v(k)(impulse).Amp*Mx;
		Magnetic_X_e(n+int(delay_total))		+= (1-a)*Mx_e;
		Magnetic_X_e(n+int(delay_total)+1)		+= a*Mx_e;
	}
	//For TE waveform V_h_v:
	for (impulse=0;impulse<V_h_v(k).size();impulse++)
	{
		delay_total = delay_lateral + V_h_v(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Mx_h = V_h_v(k)(impulse).Amp*Mx;
		Magnetic_X_h(n+int(delay_total))		+= (1-a)*Mx_h;
		Magnetic_X_h(n+int(delay_total)+1)		+= a*Mx_h;
	}
}

void Ctr_td_3l::Update_Magnetic_Y(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to My
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX+0.5)*SinTCosP+(j-FarFieldOriginY)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_e_v:
	for (impulse=0;impulse<V_e_v(k).size();impulse++)
	{
		delay_total = delay_lateral + V_e_v(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		My_e = V_e_v(k)(impulse).Amp*My;
		Magnetic_Y_e(n+int(delay_total))		+= (1-a)*My_e;
		Magnetic_Y_e(n+int(delay_total)+1)		+= a*My_e;
	}
	//For TE waveform V_h_v:
	for (impulse=0;impulse<V_h_v(k).size();impulse++)
	{
		delay_total = delay_lateral + V_h_v(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		My_h = V_h_v(k)(impulse).Amp*My;
		Magnetic_Y_h(n+int(delay_total))		+= (1-a)*My_h;
		Magnetic_Y_h(n+int(delay_total)+1)		+= a*My_h;
	}
}

void Ctr_td_3l::Update_Magnetic_Z(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to Mz
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX+0.5)*SinTCosP+(j-FarFieldOriginY+0.5)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_h_i:
	for (impulse=0;impulse<V_h_i(k).size();impulse++)
	{
		delay_total = delay_lateral + V_h_i(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Mz_h = V_h_i(k)(impulse).Amp*Mz;
		Magnetic_Z_h(n+int(delay_total))		+= (1-a)*Mz_h;
		Magnetic_Z_h(n+int(delay_total)+1)		+= a*Mz_h;
	}
}

void Ctr_td_3l::Update_Electric_X(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to Jx
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX+0.5)*SinTCosP+(j-FarFieldOriginY)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_e_i:
	for (impulse=0;impulse<V_e_i(k).size();impulse++)
	{
		delay_total = delay_lateral + V_e_i(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Jx_e = V_e_i(k)(impulse).Amp*Jx;
		Electric_X_e(n+int(delay_total))		+= (1-a)*Jx_e;
		Electric_X_e(n+int(delay_total)+1)		+= a*Jx_e;
	}
	//For TE waveform V_h_i:
	for (impulse=0;impulse<V_h_i(k).size();impulse++)
	{
		delay_total = delay_lateral + V_h_i(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Jx_h = V_h_i(k)(impulse).Amp*Jx;
		Electric_X_h(n+int(delay_total))		+= (1-a)*Jx_h;
		Electric_X_h(n+int(delay_total)+1)		+= a*Jx_h;
	}
}

void Ctr_td_3l::Update_Electric_Y(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to Jy
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX)*SinTCosP+(j-FarFieldOriginY+0.5)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_e_i:
	for (impulse=0;impulse<V_e_i(k).size();impulse++)
	{
		delay_total = delay_lateral + V_e_i(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Jy_e = V_e_i(k)(impulse).Amp*Jy;
		Electric_Y_e(n+int(delay_total))		+= (1-a)*Jy_e;
		Electric_Y_e(n+int(delay_total)+1)		+= a*Jy_e;
	}
	//For TE waveform V_h_i:
	for (impulse=0;impulse<V_h_i(k).size();impulse++)
	{
		delay_total = delay_lateral + V_h_i(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Jy_h = V_h_i(k)(impulse).Amp*Jy;
		Electric_Y_h(n+int(delay_total))		+= (1-a)*Jy_h;
		Electric_Y_h(n+int(delay_total)+1)		+= a*Jy_h;
	}
}

void Ctr_td_3l::Update_Electric_Z(const int& i, const int& j, const int& k, const int& n)
{//Update the far-field waveforms due to Jz
	delay_lateral = (r_offset-dx*((i-FarFieldOriginX)*SinTCosP+(j-FarFieldOriginY)*SinTSinP))/c_o/dt;	//Lateral delay
	//For TM waveform V_e_v:
	for (impulse=0;impulse<V_e_v(k).size();impulse++)
	{
		delay_total = delay_lateral + V_e_v(k)(impulse).Delay;
		a = delay_total-int(delay_total);
		Jz_e = epsilon_o/(eps_z(Layering_e_z(k)))*V_e_v(k)(impulse).Amp*Jz;
		Electric_Z_e(n+int(delay_total))		+= (1-a)*Jz_e;
		Electric_Z_e(n+int(delay_total)+1)		+= a*Jz_e;
	}
}

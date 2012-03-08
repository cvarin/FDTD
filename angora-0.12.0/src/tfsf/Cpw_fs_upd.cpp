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

// Includes the methods that update the auxiliary grid associated with a single plane wave in free space.

#include "headers.h"

#include "Cpw_fs.h"

extern int i,j,k;


void Cpw_fs::UpdateIncidentH(const int& n)
{
	//H-field is already at time n+1/2, since 1-D grid has the same dt as the main grid
	//no update necessary
}

void Cpw_fs::UpdateIncidentE(const int& n)
{
	//update incident E and H until E-field is at time n+1
	//Update the E-field
	UpdateE(n);
	//Update the H-field
	UpdateH(n);
}

void Cpw_fs::UpdateE(const int& n)
{//updates the E-field components at time index n
	for (k=AbsorbPointLower+1;k<=HardSourcePoint-1;k++)
	{
		AuxGrid_E(k) = Ca*AuxGrid_E(k) + Cb*(AuxGrid_H(k)-AuxGrid_H(k-1));
		//Note the "+" sign before Cb, instead of "-" in Cpw_fs_update.cpp
		//In the multilayered version, the characteristic impedance of the line did not matter, since the directions of the field components in the 1-D grid were well-defined in the main grid. We only had to make sure that the updates were consistent with the TL equations. The rest was taken care of by the 1-D TF/SF corrections in the 1-D grid.
		//Here, however, we want the characteristic impedance to be positive for a wave that is propagating downward in the 1-D grid (which is our assumption). Otherwise, the H-field in the 1-D grid would have the opposite sign relative to the E-field for a downward-propagating wave, which is inconsistent with the desired plane wave in the  main grid.
	}

	//Apply hard-source for Einc
	AuxGrid_E(HardSourcePoint) = waveformE_value(n);

	//Apply 2nd order ABC
	//NOTE: PreviousELower(0,0) need not be used, since it is equal to AuxGrid_E(AbsorbPointLower) [before the update]
	//, but it is used for symmetry.
	//Lower absorbing point
	AuxGrid_E(AbsorbPointLower) = (-1/(1/Sc_1D_lower+2+Sc_1D_lower))*(
		(1/Sc_1D_lower-2+Sc_1D_lower)*(AuxGrid_E(AbsorbPointLower+2)+PreviousELower(0,-1))
		+ 2*(Sc_1D_lower-1/Sc_1D_lower)*(PreviousELower(0,0)+PreviousELower(2,0)-AuxGrid_E(AbsorbPointLower+1)-PreviousELower(1,-1))
		- 4*(1/Sc_1D_lower+Sc_1D_lower)*PreviousELower(1,0)) - PreviousELower(2,-1);
	//update the field history values
	PreviousELower(Range(0,2),-1) = PreviousELower(Range(0,2),0);	//shift toward the future by one step
	PreviousELower(Range(0,2),0) = AuxGrid_E(Range(AbsorbPointLower,AbsorbPointLower+2));	//record the current values as history

	//calculate the x,y,z components separately
	AuxGrid_Ex = k_Ex*AuxGrid_E;
	AuxGrid_Ey = k_Ey*AuxGrid_E;
	AuxGrid_Ez = k_Ez*AuxGrid_E;
}

void Cpw_fs::UpdateH(const int& n)
{//updates the H-field components at time index n
	for (k=AbsorbPointLower;k<=HardSourcePoint-1;k++)
	{
		AuxGrid_H(k) = Da*AuxGrid_H(k) + Db*(AuxGrid_E(k+1)-AuxGrid_E(k));
		//See note in UpdateE regarding the "+" sign
	}

	//calculate the x,y,z components separately
	AuxGrid_Hx = k_Hx*AuxGrid_H;
	AuxGrid_Hy = k_Hy*AuxGrid_H;
	AuxGrid_Hz = k_Hz*AuxGrid_H;
}

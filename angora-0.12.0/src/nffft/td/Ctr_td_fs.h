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

#ifndef CTR_TD_FS_H
#define CTR_TD_FS_H

//Declaration of the class "Ctr_td_fs" for a TIME_DOMAIN near-field-to-far-field transformer in free space
//Derived from the abstract base class "Ctr_td".

#include "Ctr_td.h"


class Ctr_td_fs : public Ctr_td
{
 public:
	 Ctr_td_fs(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index);		//constructor

	 void UpdateFarField(const int& n);		//Update auxiliary far-field waveforms using the current E,H on the virtual surface
	 void ConstructFarField();				//Construct far-field waveforms using the final auxiliary waveforms

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

	 double TheoreticalFarFieldWaveform(double t, double tau);	//Return the appropriate Gaussian waveform, depending on the source current waveform
};

#endif

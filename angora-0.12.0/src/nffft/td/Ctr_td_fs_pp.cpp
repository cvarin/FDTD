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

//Defines the post-processing steps in the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_fs" in free space

#include "headers.h"

#include "Ctr_td_fs.h"


void Ctr_td_fs::ConstructFarField()
{
}

double Ctr_td_fs::TheoreticalFarFieldTheta(const int& n)
//Theoretical theta-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr)
{
}

double Ctr_td_fs::TheoreticalFarFieldPhi(const int& n)
//Theoretical phi-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr)
{
}

double Ctr_td_fs::TheoreticalFarFieldWaveform(double t, double tau)
{//Return the appropriate far-field waveform shape (normalized), depending on the source current waveform
}

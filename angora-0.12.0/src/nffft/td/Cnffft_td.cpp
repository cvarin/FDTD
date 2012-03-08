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

//Defines the TIME_DOMAIN near-field-to-far-field transformer (NFFFT) object "Cnffft_td", which represents a collection of transformer objects "Ctr_td".

#include "headers.h"

#include "Cnffft_td.h"

#include "Ctr_td_fs.h"
#include "Ctr_td_3l.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern int number_of_layers;

extern int rank;


int Cnffft_td::AddTransformer(const TrDataType_td& MyData, const string& FarFieldFileName)
{
	if (number_of_layers>3)
	{
		if (rank==0)
		{
			cout << endl << "Warning: Time-domain NFFFT cannot be used for more than 3 planar material layers. Skipping... " << endl;
		}
	}
	else
	{
		//attach a near-field-to-far-field transformer for a 3-layered medium
		boost::shared_ptr<Ctr_td> new_tr_ptr(new Ctr_td_3l(MyData,FarFieldFileName,NumberOfTransformers()));
		Transformers.push_back(new_tr_ptr);
		//near-field-to-far-field transformer for free space or general N-layered media can be implemented later.
	}
	return NumberOfTransformers()-1;
}

void Cnffft_td::UpdateFarField(const int& n)
{
	for (int i=0; i<NumberOfTransformers(); i++)
	{
		Transformers[i]->UpdateFarField(n);		//update the far field arrays for each transformer
	}
}

void Cnffft_td::ConstructFarField()
{
	for (int i=0; i<NumberOfTransformers(); i++)
	{
		Transformers[i]->ConstructFarField();		//construct the far field arrays for each transformer
	}
}

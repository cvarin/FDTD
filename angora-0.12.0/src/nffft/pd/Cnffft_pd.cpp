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

//Defines the PHASOR-DOMAIN near-field-to-far-field transformer (NFFFT) object "Cnffft_pd", which represents a collection of transformer objects "Ctr_pd".

#include "headers.h"

#include "Cnffft_pd.h"

#include "Ctr_pd_fs.h"
#include "Ctr_pd_2l.h"
#include "Ctr_pd_3l.h"

extern int number_of_layers;

extern int rank;


int Cnffft_pd::AddTransformer(const TrDataType_pd& MyData, const string& FarFieldFileName)
{
	if (number_of_layers>3)
	{
		boost::shared_ptr<Ctr_pd> new_tr_ptr(new Ctr_pd_ml(MyData,FarFieldFileName,NumberOfTransformers()));
//cout << "Added ml transformer for " << number_of_layers << "-layered medium!!" << endl;
		Transformers.push_back(new_tr_ptr);
	}
	else
	{
		if (number_of_layers==1)
		{//attach a free-space transformer
			boost::shared_ptr<Ctr_pd> new_tr_ptr(new Ctr_pd_fs(MyData,FarFieldFileName,NumberOfTransformers()));
			Transformers.push_back(new_tr_ptr);
		}
		else if (number_of_layers==2)
		{//attach a 2-layered medium transformer
			boost::shared_ptr<Ctr_pd> new_tr_ptr(new Ctr_pd_2l(MyData,FarFieldFileName,NumberOfTransformers()));
			Transformers.push_back(new_tr_ptr);
		}
		else if (number_of_layers==3)
		{//attach a 3-layered medium transformer
			boost::shared_ptr<Ctr_pd> new_tr_ptr(new Ctr_pd_3l(MyData,FarFieldFileName,NumberOfTransformers()));
			Transformers.push_back(new_tr_ptr);
		}
	}
	return NumberOfTransformers()-1;
}

void Cnffft_pd::UpdateFarField(const int& n)
{
	for (int i=0; i<NumberOfTransformers(); i++)
	{
		Transformers[i]->UpdateFarField(n);		//update the far field arrays for each transformer
	}
}

void Cnffft_pd::ConstructFarField()
{
	if ((rank==0)&&(NumberOfTransformers()>0))
	{
		cout << "Constructing far-field arrays..." << endl;
	}
	for (int i=0; i<NumberOfTransformers(); i++)
	{
		Transformers[i]->ConstructFarField();		//construct the far field arrays for each transformer
	}
}

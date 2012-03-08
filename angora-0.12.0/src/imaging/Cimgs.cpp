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

//Definitions for the collector class "Cimgs" for optical image objects Cimg

#include "headers.h"

#include "Cimgs.h"

#include "Cimg.h"

extern int rank;


int Cimgs::AddOpticalImage(const ImgDataType& MyData, const string& ImgFileName)
{
	boost::shared_ptr<Cimg> new_image_ptr(new Cimg(MyData,ImgFileName,Images.size()));
	Images.push_back(new_image_ptr);
}

void Cimgs::UpdateFarField(const int& n)
{
	for (int i=0; i<Images.size(); i++)
	{
		Images[i]->UpdateFarField(n);		//update the far field arrays for each optical image
	}
}

void Cimgs::ConstructImages()
{
	if ((rank==0)&&(Images.size()>0))
	{
		cout << "Constructing optical images..." << endl;
	}
	for (int i=0; i<Images.size(); i++)
	{
		Images[i]->ConstructImage();		//construct each optical image
	}
}

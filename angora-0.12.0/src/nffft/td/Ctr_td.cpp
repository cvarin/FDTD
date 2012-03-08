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

//Definition of the abstract base class "Ctr_td" for a TIME_DOMAIN near-field-to-far-field transformer

#include "headers.h"

#include "Ctr_td.h"

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern int rank;


Ctr_td::Ctr_td(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index)
		: Data(MyData), TransformerIndex(Index)
{
	//Far-field observation angles
	THETA = Data.THETA;	//must be between -M_PI/2 and M_PI/2
	PHI = Data.PHI;		//must be between 0 and 2*M_PI

	SinTCosP = sin(THETA)*cos(PHI);
	SinTSinP = sin(THETA)*sin(PHI);
	SinT = sin(THETA);
	CosT = cos(THETA);
	SinP = sin(PHI);
	CosP = cos(PHI);

	//Virtual surface dimensions
	SurfaceBackX=NPML+Data.NFFFTMarginBackX;					//x-index of the cells neighboring the back face
	SurfaceFrontX=(NCELLS_X+2*NPML)-NPML-Data.NFFFTMarginFrontX;		//x-index of the cells neighboring the front face
	SurfaceLeftY=NPML+Data.NFFFTMarginLeftY;							//y-index of the cells neighboring the right face
	SurfaceRightY=(NCELLS_Y+2*NPML)-NPML-Data.NFFFTMarginRightY;		//y-index of the cells neighboring the right face
	SurfaceLowerZ=NPML+Data.NFFFTMarginLowerZ;						//z-index of the cells neighboring the lower face
	SurfaceUpperZ=(NCELLS_Z+2*NPML)-NPML-Data.NFFFTMarginUpperZ;		//z-index of the cells neighboring the upper face

	if (rank==0)
	{
		FarFieldFile.open(FarFieldFileName.c_str(),ios::binary);
		if (!FarFieldFile)
		{
			cout << "Error opening far-field output file!" << endl << endl;
			exit(-1);
		}
	}
}

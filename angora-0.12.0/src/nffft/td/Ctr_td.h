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

#ifndef CTR_TD_H
#define CTR_TD_H

//Declaration of the abstract base class "Ctr_td" for a TIME_DOMAIN near-field-to-far-field transformer

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

#include <fstream>

//only the declaration of Cpointsources needed: use forward declaration
class Cpointsources;

extern int OriginX,OriginY,OriginZ;


class TrDataType_td
{
 public:
	 TrDataType_td():
	   	THETA(0*M_PI/180),
	   	PHI(0*M_PI/180),
		NFFFTMarginBackX(3),
		NFFFTMarginFrontX(3),
		NFFFTMarginLeftY(3),
		NFFFTMarginRightY(3),
		NFFFTMarginLowerZ(3),
		NFFFTMarginUpperZ(3),
		NFFFTOriginX(OriginX),
		NFFFTOriginY(OriginY),
		NFFFTOriginZ(OriginZ),
	   	PointSourcesPtr(NULL)
	   {};

	 double THETA,PHI;	//observation angles
	 //distances (in cells) of the NFFFT box from the PML boundary
	 int NFFFTMarginBackX,NFFFTMarginFrontX,NFFFTMarginLeftY,NFFFTMarginRightY,NFFFTMarginLowerZ,NFFFTMarginUpperZ;
	 //the index of the cell assigned as the origin of the NFFFT box
	 //(its rear-left-lower corner is the origin)
	 double NFFFTOriginX,NFFFTOriginY,NFFFTOriginZ;		//(can be fractional, therefore defined as "double")
	 const Cpointsources* PointSourcesPtr;	//points to a non-changeable Cpointsources object
};

class Ctr_td
{
 public:
	 Ctr_td(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index);	//constructor
 	 virtual ~Ctr_td(){};	//virtual destructor for cleaning up some objects (e.g. ofstream object)
						//this is needed when deleting derived objects using a base class pointer

	 virtual void UpdateFarField(const int& n)=0;		//Do some auxiliary updates using the current E,H on the virtual surface
	 virtual void ConstructFarField()=0;	//Do the necessary post-processing steps to obtain the final far field

 protected:
	 virtual double TheoreticalFarFieldTheta(const int& n)=0;	//the theoretical theta-component of the far field
	 virtual double TheoreticalFarFieldPhi(const int& n)=0;		//the theoretical phi-component of the far field

	 const int TransformerIndex;	//index of the current transformer in the NFFFT object Cnffft_td (1st, 2nd,.. etc.)

	 const TrDataType_td Data;	//data type that holds the transformer data

	 //NFFFT observation angles
	 double THETA,PHI;
	 //Direction sines and cosines
	 double SinTCosP,SinTSinP,SinT,CosT,SinP,CosP;

	 //Location of the transform box
	 int SurfaceBackX,SurfaceFrontX,SurfaceLeftY,SurfaceRightY,SurfaceLowerZ,SurfaceUpperZ;
	 //the index of the cell assigned as the origin of the NFFFT box
	 //(its rear-left-lower corner is the origin)
	 double FarFieldOriginX,FarFieldOriginY,FarFieldOriginZ;

	 Array<double,1> waveformTheta;	//Theta component of the radiated electric field (normalized by 1/r, advanced by r+rOffset)
	 Array<double,1> waveformPhi;	//Phi component of the radiated electric field (normalized by 1/r, advanced by r+rOffset)

#ifndef MPI_DISABLE
	 MPI_Status Status;
#endif
	 enum FarFieldTag {ThetaTag=0,PhiTag=1};

	 ofstream FarFieldFile;		//far-field file

	 //Theoretical far-field calculations
	 double theta_component,phi_component;
	 double temp;
	 string Orientation;
// 	 string WaveShape;
	 int x,y,z;
	 double SourceX,SourceY;
/*	 double j0;
	 double delay;
	 double tau;
	 virtual double TheoreticalFarFieldWaveform(double t, double tau)=0;	//Return the appropriate Gaussian waveform, depending on the source current waveform*/

	 double a;	//interpolation parameter

	 double Jx,Jy,Jz,Mx,My,Mz;
};

#endif

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

// Includes the methods that perform the TF/SF corrections associated with a single plane wave in a multilayered medium.

#include "headers.h"

#include "Cpw_ml.h"

//for the definition of MaterialId
#include "material_id.h"

extern double dx,dt;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern double c_upper;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int i,j,k;


// ApplyCorrectionE and ApplyCorrectionH are actually the same as Cpw_fs. Instead of inheriting them from Cpw_fs, they are copied here.
// This is done to be able to define IncidentEx,etc.inline, since inline functions have to be in the same file as the functions from which they are called.

void Cpw_ml::ApplyCorrectionE(const int& n)
{
//NOTE: If a TF/SF E-field position is included in the node, it must be updated, whether or not it is at the boundary of the node. However, a TF/SF H-field position at the boundary of the node (immediately outside the outer E-field shell of the node) need not be updated, since it will soon be overwritten by the H-field from the next node.

	//Bottom face of the plane-wave source TF/SF box
	if ((klower<=PWLowerZ)&&(kupper>=PWLowerZ-1))	//is the lower face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		k=PWLowerZ;
		//Correct Ex
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (j=TFSF_min_y;j<=TFSF_max_y+1;j++)
			{
				Hinc = IncidentHy(i,j,k-1,n);		//calculate incident Hy by projection
				Ex(i,j,k) += Cb_X(Media_Ex(i,j,k))*Hinc;
			}
		}
		//Correct Ey
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (j=TFSF_min_y;j<=TFSF_max_y;j++)
			{
				Hinc = IncidentHx(i,j,k-1,n);		//calculate incident Hx by projection
				Ey(i,j,k) += -Cb_Y(Media_Ey(i,j,k))*Hinc;
			}
		}
	}

	//Upper face of the plane-wave source TF/SF box
	if ((klower<=PWUpperZ+1)&&(kupper>=PWUpperZ))	//is the upper face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		k=PWUpperZ+1;
		//Correct Ex
		for (i=TFSF_min_x;i<=TFSF_max_x;i++)
		{
			for (j=TFSF_min_y;j<=TFSF_max_y+1;j++)
			{
				Hinc = IncidentHy(i,j,k,n);		//calculate incident Hy by projection
				Ex(i,j,k) += -Cb_X(Media_Ex(i,j,k))*Hinc;
			}
		}
		//Correct Ey
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++)
		{
			for (j=TFSF_min_y;j<=TFSF_max_y;j++)
			{
				Hinc = IncidentHx(i,j,k,n);		//calculate incident Hx by projection
				Ey(i,j,k) += Cb_Y(Media_Ey(i,j,k))*Hinc;
			}
		}
	}

	//Back face of the plane-wave source TF/SF box
	if ((iback<=PWBackX)&&(ifront>=PWBackX-1))	//is the back face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		i=PWBackX;
		//Correct Ey
		for (j=TFSF_min_y;j<=TFSF_max_y;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Hinc = IncidentHz(i-1,j,k,n);		//calculate incident Hz by projection
				Ey(i,j,k) += Cb_Y(Media_Ey(i,j,k))*Hinc;
			}
		}
		//Correct Ez
		for (j=TFSF_min_y;j<=TFSF_max_y+1;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Hinc = IncidentHy(i-1,j,k,n);		//calculate incident Hy by projection
				Ez(i,j,k) += -Cb_Z(Media_Ez(i,j,k))*Hinc;
			}
		}
	}

	//Front face of the plane-wave source TF/SF box
	if ((iback<=PWFrontX+1)&&(ifront>=PWFrontX))	//is the front face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		i=PWFrontX+1;
		//Correct Ey
		for (j=TFSF_min_y;j<=TFSF_max_y;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Hinc = IncidentHz(i,j,k,n);		//calculate incident Hz by projection
				Ey(i,j,k) += -Cb_Y(Media_Ey(i,j,k))*Hinc;
			}
		}
		//Correct Ez
		for (j=TFSF_min_y;j<=TFSF_max_y+1;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Hinc = IncidentHy(i,j,k,n);		//calculate incident Hy by projection
				Ez(i,j,k) += Cb_Z(Media_Ez(i,j,k))*Hinc;
			}
		}
	}

	//Left face of the plane-wave source TF/SF box
	if ((jleft<=PWLeftY)&&(jright>=PWLeftY-1))	//is the left face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		j=PWLeftY;
		//Correct Ex
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Hinc = IncidentHz(i,j-1,k,n);		//calculate incident Hz by projection
				Ex(i,j,k) += -Cb_X(Media_Ex(i,j,k))*Hinc;
			}
		}
		//Correct Ez
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Hinc = IncidentHx(i,j-1,k,n);		//calculate incident Hx by projection
				Ez(i,j,k) += Cb_Z(Media_Ez(i,j,k))*Hinc;
			}
		}
	}

	//Right face of the plane-wave source TF/SF box
	if ((jleft<=PWRightY+1)&&(jright>=PWRightY))	//is the right face of the box included in the node? (it can be included simultaneously in two neighboring nodes, and updated by both)
	{
		j=PWRightY+1;
		//Correct Ex
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Hinc = IncidentHz(i,j,k,n);		//calculate incident Hz by projection
				Ex(i,j,k) += Cb_X(Media_Ex(i,j,k))*Hinc;
			}
		}
		//Correct Ez
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Hinc = IncidentHx(i,j,k,n);		//calculate incident Hx by projection
				Ez(i,j,k) += -Cb_Z(Media_Ez(i,j,k))*Hinc;
			}
		}
	}
}

void Cpw_ml::ApplyCorrectionH(const int& n)
{

	//Bottom face of the plane-wave source TF/SF box
	if ((klower<=PWLowerZ-1)&&(kupper>=PWLowerZ-1))	//is the lower face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		k=PWLowerZ-1;
		//Correct Hx
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (j=TFSF_min_y;j<=TFSF_max_y;j++)
			{
				Einc = IncidentEy(i,j,k+1,n);		//calculate incident Ey by projection
				Hx(i,j,k) += -Db_X(Media_Hx(i,j,k))*Einc;
			}
		}
		//Correct Hy
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (j=TFSF_min_y;j<=TFSF_max_y+1;j++)
			{
				Einc = IncidentEx(i,j,k+1,n);		//calculate incident Ex by projection
				Hy(i,j,k) += Db_Y(Media_Hy(i,j,k))*Einc;
			}
		}
	}

	//Upper face of the plane-wave source TF/SF box
	if ((klower<=PWUpperZ+1)&&(kupper>=PWUpperZ+1))	//is the upper face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		k=PWUpperZ+1;
		//Correct Hx
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++)
		{
			for (j=TFSF_min_y;j<=TFSF_max_y;j++)
			{
				Einc = IncidentEy(i,j,k,n);		//calculate incident Ey by projection
				Hx(i,j,k) += Db_X(Media_Hx(i,j,k))*Einc;
			}
		}
		//Correct Hy
		for (i=TFSF_min_x;i<=TFSF_max_x;i++)
		{
			for (j=TFSF_min_y;j<=TFSF_max_y+1;j++)
			{
				Einc = IncidentEx(i,j,k,n);		//calculate incident Ex by projection
				Hy(i,j,k) += -Db_Y(Media_Hy(i,j,k))*Einc;
			}
		}
	}

	//Back face of the plane-wave source TF/SF box
	if ((iback<=PWBackX-1)&&(ifront>=PWBackX-1))	//is the back face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		i=PWBackX-1;
		//Correct Hy
		for (j=TFSF_min_y;j<=TFSF_max_y+1;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Einc = IncidentEz(i+1,j,k,n);		//calculate incident Ez by projection
				Hy(i,j,k) += -Db_Y(Media_Hy(i,j,k))*Einc;
			}
		}
		//Correct Hz
		for (j=TFSF_min_y;j<=TFSF_max_y;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Einc = IncidentEy(i+1,j,k,n);		//calculate incident Ey by projection
				Hz(i,j,k) += Db_Z(Media_Hz(i,j,k))*Einc;
			}
		}
	}

	//Front face of the plane-wave source TF/SF box
	if ((iback<=PWFrontX+1)&&(ifront>=PWFrontX+1))	//is the front face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		i=PWFrontX+1;
		//Correct Hy
		for (j=TFSF_min_y;j<=TFSF_max_y+1;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Einc = IncidentEz(i,j,k,n);		//calculate incident Ez by projection
				Hy(i,j,k) += Db_Y(Media_Hy(i,j,k))*Einc;
			}
		}
		//Correct Hz
		for (j=TFSF_min_y;j<=TFSF_max_y;j++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Einc = IncidentEy(i,j,k,n);		//calculate incident Ey by projection
				Hz(i,j,k) += -Db_Z(Media_Hz(i,j,k))*Einc;
			}
		}
	}

	//Left face of the plane-wave source TF/SF box
	if ((jleft<=PWLeftY-1)&&(jright>=PWLeftY-1))	//is the left face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		j=PWLeftY-1;
		//Correct Hx
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Einc = IncidentEz(i,j+1,k,n);		//calculate incident Ez by projection
				Hx(i,j,k) += Db_X(Media_Hx(i,j,k))*Einc;
			}
		}
		//Correct Hz
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Einc = IncidentEx(i,j+1,k,n);		//calculate incident Ex by projection
				Hz(i,j,k) += -Db_Z(Media_Hz(i,j,k))*Einc;
			}
		}
	}

	//Right face of the plane-wave source TF/SF box
	if ((jleft<=PWRightY+1)&&(jright>=PWRightY+1))	//is the right face of the box included in the node? (correct H-field if you are the node who gives this H-field information to the next node, during the MPI H-field exchange)
	{
		j=PWRightY+1;
		//Correct Hx
		for (i=TFSF_min_x;i<=TFSF_max_x+1;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z;k++)
			{
				Einc = IncidentEz(i,j,k,n);		//calculate incident Ez by projection
				Hx(i,j,k) += -Db_X(Media_Hx(i,j,k))*Einc;
			}
		}
		//Correct Hz
		for (i=TFSF_min_x;i<=TFSF_max_x;i++){
			for (k=TFSF_min_z;k<=TFSF_max_z+1;k++)
			{
				Einc = IncidentEx(i,j,k,n);		//calculate incident Ex by projection
				Hz(i,j,k) += Db_Z(Media_Hz(i,j,k))*Einc;
			}
		}
	}
}

//******************************************************************
//  PLANE-WAVE INCIDENCE IN LAYERED MEDIUM (PROJECTION ON 1-D GRID)
//******************************************************************
//These inline functions have to be in the same file as the functions (ApplyCorrectionE,H) from which they are called.
inline double Cpw_ml::IncidentEx(const int& i, const int& j, const int& k, const int& n)
{
	Einc_e = 0;
	Einc_h = 0;
	Coordinate = (i-0.5),(j-1),(k-1);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Einc_e = (1-a)*TimeHistory_Ee(DelayedIntegerTime,k)+a*TimeHistory_Ee(DelayedIntegerTime+1,k);
		Einc_h = (1-a)*TimeHistory_Eh(DelayedIntegerTime,k)+a*TimeHistory_Eh(DelayedIntegerTime+1,k);
	}
	return Einc_e*k_Ee(firstDim)+Einc_h*k_Eh(firstDim);
};

inline double Cpw_ml::IncidentEy(const int& i, const int& j, const int& k, const int& n)
{
	Einc_e = 0;
	Einc_h = 0;
	Coordinate = (i-1),(j-0.5),(k-1);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Einc_e = (1-a)*TimeHistory_Ee(DelayedIntegerTime,k)+a*TimeHistory_Ee(DelayedIntegerTime+1,k);
		Einc_h = (1-a)*TimeHistory_Eh(DelayedIntegerTime,k)+a*TimeHistory_Eh(DelayedIntegerTime+1,k);
	}
	return Einc_e*k_Ee(secondDim)+Einc_h*k_Eh(secondDim);
};

inline double Cpw_ml::IncidentEz(const int& i, const int& j, const int& k, const int& n)
{
	Einc_z = 0;
	Coordinate = (i-1),(j-1),(k-0.5);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Einc_z = (1-a)*TimeHistory_Ez(DelayedIntegerTime,k)+a*TimeHistory_Ez(DelayedIntegerTime+1,k);
	}
	return Einc_z;
};

inline double Cpw_ml::IncidentHx(const int& i, const int& j, const int& k, const int& n)
{
	Hinc_e = 0;
	Hinc_h = 0;
	Coordinate = (i-1),(j-0.5),(k-0.5);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Hinc_e = (1-a)*TimeHistory_He(DelayedIntegerTime,k)+a*TimeHistory_He(DelayedIntegerTime+1,k);
		Hinc_h = (1-a)*TimeHistory_Hh(DelayedIntegerTime,k)+a*TimeHistory_Hh(DelayedIntegerTime+1,k);
	}
	return Hinc_e*k_He(firstDim)+Hinc_h*k_Hh(firstDim);
};

inline double Cpw_ml::IncidentHy(const int& i, const int& j, const int& k, const int& n)
{
	Hinc_e = 0;
	Hinc_h = 0;
	Coordinate = (i-0.5),(j-1),(k-0.5);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Hinc_e = (1-a)*TimeHistory_He(DelayedIntegerTime,k)+a*TimeHistory_He(DelayedIntegerTime+1,k);
		Hinc_h = (1-a)*TimeHistory_Hh(DelayedIntegerTime,k)+a*TimeHistory_Hh(DelayedIntegerTime+1,k);
	}
	return Hinc_e*k_He(secondDim)+Hinc_h*k_Hh(secondDim);
};

inline double Cpw_ml::IncidentHz(const int& i, const int& j, const int& k, const int& n)
{
	Hinc_z = 0;
	Coordinate = (i-0.5),(j-0.5),(k-1);
	Coordinate -= PWOriginationPoint;
	PositionOnPrincipal = dot(Coordinate,k_inc_lateral);
	TimeDelay = PositionOnPrincipal*dx*abs(SinT)/c_upper/dt;	//time delay (in dt) compared to the reference 1-D grid
	DelayedTime = n - TimeDelay;					//corresponding time index in time history
	DelayedIntegerTime = int(DelayedTime);
	if (DelayedTime >= 0)	//if time index is nonnegative, interpolate field value in time history
	{
		a = DelayedTime - DelayedIntegerTime;
		Hinc_z = (1-a)*TimeHistory_Hz(DelayedIntegerTime,k)+a*TimeHistory_Hz(DelayedIntegerTime+1,k);
	}
	return Hinc_z;
};

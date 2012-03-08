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

//Includes routines that carry out the main-grid updates.

#include "headers.h"

#include "update.h"

//for definition of MaterialId
#include "material_id.h"

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
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
extern Array<double,1> kappa_e_x,kappa_e_y,kappa_e_z,kappa_h_x,kappa_h_y,kappa_h_z;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int FullIntPosMin_x,FullIntPosMax_x;
extern int FullIntPosMin_y,FullIntPosMax_y;
extern int FullIntPosMin_z,FullIntPosMax_z;

extern int i,j,k;

void updateEx(const int& n);
void updateEy(const int& n);
void updateEz(const int& n);
void updateHx(const int& n);
void updateHy(const int& n);
void updateHz(const int& n);


//**********************************
//	UPDATE THE ELECTRIC FIELD
//**********************************

void updateE(const int& n)
{
	updateEx(n);
	updateEy(n);
	updateEz(n);
}

//**********************************
//	UPDATE THE MAGNETIC FIELD
//**********************************
void updateH(const int& n)
{
	updateHx(n);
	updateHy(n);
	updateHz(n);
}

void updateEx(const int& n)
{
	//Update Ex
	for (i=iback; i<=ifront; i++){
		for (j=FullIntPosMin_y; j<=FullIntPosMax_y; j++){
			for (k=FullIntPosMin_z; k<=FullIntPosMax_z; k++)
			{
				Ex(i,j,k)=Ca_X(Media_Ex(i,j,k))*Ex(i,j,k)
					+ Cb_X(Media_Ex(i,j,k))*
					((Hz(i,j,k)-Hz(i,j-1,k))/kappa_e_y(j)
					+(Hy(i,j,k-1)-Hy(i,j,k))/kappa_e_z(k));
			}
		}
	}
}

void updateEy(const int& n)
{
	//Update Ey
	for (i=FullIntPosMin_x; i<=FullIntPosMax_x; i++){
		for (j=jleft; j<=jright; j++){
			for (k=FullIntPosMin_z; k<=FullIntPosMax_z; k++)
			{
				Ey(i,j,k)=Ca_Y(Media_Ey(i,j,k))*Ey(i,j,k)
					+ Cb_Y(Media_Ey(i,j,k))*
					((Hx(i,j,k)-Hx(i,j,k-1))/kappa_e_z(k)
					+(Hz(i-1,j,k)-Hz(i,j,k))/kappa_e_x(i));
			}
		}
	}
}

void updateEz(const int& n)
{
	//Update Ez
	for (i=FullIntPosMin_x; i<=FullIntPosMax_x; i++){
		for (j=FullIntPosMin_y; j<=FullIntPosMax_y; j++){
			for (k=klower; k<=kupper; k++)
			{
				Ez(i,j,k)=Ca_Z(Media_Ez(i,j,k))*Ez(i,j,k)
					+ Cb_Z(Media_Ez(i,j,k))*
					((Hy(i,j,k)-Hy(i-1,j,k))/kappa_e_x(i)
						+(Hx(i,j-1,k)-Hx(i,j,k))/kappa_e_y(j));
			}
		}
	}
}

void updateHx(const int& n)
{
	//Update Hx
	for (i=FullIntPosMin_x; i<=FullIntPosMax_x; i++){
		for (j=jleft; j<=jright; j++){
			for (k=klower; k<=kupper; k++)
			{
				Hx(i,j,k)=Da_X(Media_Hx(i,j,k))*Hx(i,j,k)
					+ Db_X(Media_Hx(i,j,k))*
					((Ey(i,j,k+1)-Ey(i,j,k))/kappa_h_z(k)
						+(Ez(i,j,k)-Ez(i,j+1,k))/kappa_h_y(j));
			}
		}
	}
}

void updateHy(const int& n)
{
	//Update Hy
	for (i=iback; i<=ifront; i++){
		for (j=FullIntPosMin_y; j<=FullIntPosMax_y; j++){
			for (k=klower; k<=kupper; k++)
			{
				Hy(i,j,k)=Da_Y(Media_Hy(i,j,k))*Hy(i,j,k)
					+ Db_Y(Media_Hy(i,j,k))*
					((Ez(i+1,j,k)-Ez(i,j,k))/kappa_h_x(i)
						+(Ex(i,j,k)-Ex(i,j,k+1))/kappa_h_z(k));
			}
		}
	}
}

void updateHz(const int& n)
{
    //Update Hz
	for (i=iback; i<=ifront; i++){
		for (j=jleft; j<=jright; j++){
			for (k=FullIntPosMin_z; k<=FullIntPosMax_z; k++)
			{
				Hz(i,j,k)=Da_Z(Media_Hz(i,j,k))*Hz(i,j,k)
					+ Db_Z(Media_Hz(i,j,k))*
					((Ex(i,j+1,k)-Ex(i,j,k))/kappa_h_y(j)
						+(Ey(i,j,k)-Ey(i+1,j,k))/kappa_h_x(i));
			}
		}
	}
}

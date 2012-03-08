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

//Declaration of the class "Ctr_td_fs" for a TIME_DOMAIN near-field-to-far-field transformer in free space

#include "headers.h"

#include "Ctr_td_fs.h"

extern Array<double,3> Ex,Ey,Ez,Hx,Hy,Hz;

extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int i,j,k;


Ctr_td_fs::Ctr_td_fs(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_td(MyData,FarFieldFileName,Index)
{
}

void Ctr_td_fs::UpdateFarField(const int& n)
{
	 UpdateElectric(n);
	 UpdateMagnetic(n);
}

//The following code is purposefully repeated in other classes derived from Ctr_td, instead of using a common definition. The use of virtual functions for Update_Magnetic_Y, etc. would decrease the efficiency.
void Ctr_td_fs::UpdateMagnetic(const int& n)
{
	//Using the equivalence theorem, the E-values on the virtual surface are converted to Meq-values.
	//The Meq values are half-cell away from the virtual surface, and are used to update Magnetic_(X,Y,Z).
	//It must be ensured that each E-value updates Magnetic_(X,Y,Z) in only ONE NODE.

	//Bottom face of the virtual surface - Use Ex, Ey
	if ((klower<=SurfaceLowerZ)&&(kupper>=SurfaceLowerZ)) //is the bottom face of the box included in the node?
		//lowermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		k=SurfaceLowerZ;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
			{
				My = Ex(i,j,k);
				Update_Magnetic_Y(i,j,k-1,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
			{
				Mx = -Ey(i,j,k);
				Update_Magnetic_X(i,j,k-1,n);
			}
		}
	}

	//Upper face of the virtual surface - Use Ex, Ey
	if ((klower<=SurfaceUpperZ)&&(kupper>=SurfaceUpperZ)) //is the upper face of the box included in the node?
		//uppermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		k=SurfaceUpperZ+1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
			{
				My = -Ex(i,j,k);
				Update_Magnetic_Y(i,j,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
			{
				Mx = Ey(i,j,k);
				Update_Magnetic_X(i,j,k,n);
			}
		}
	}

	//Back face of the virtual surface - Use Ey, Ez
	if ((iback<=SurfaceBackX)&&(ifront>=SurfaceBackX))	//is the back face of the box included in the node?
		//lowermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		i=SurfaceBackX;
		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Mz = Ey(i,j,k);
				Update_Magnetic_Z(i-1,j,k,n);
			}
		}

		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				My = -Ez(i,j,k);
				Update_Magnetic_Y(i-1,j,k,n);
			}
		}
	}

	//Front face of the virtual surface - Use Ey, Ez
	if ((iback<=SurfaceFrontX)&&(ifront>=SurfaceFrontX)) //is the front face of the box included in the node?
		//lowermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		i=SurfaceFrontX+1;
		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Mz = -Ey(i,j,k);
				Update_Magnetic_Z(i,j,k,n);
			}
		}

		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				My = Ez(i,j,k);
				Update_Magnetic_Y(i,j,k,n);
			}
		}
	}

	//Left face of the virtual surface - Use Ex, Ez
	if ((jleft<=SurfaceLeftY)&&(jright>=SurfaceLeftY))	//is the left face of the box included in the node?
		//lowermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		j=SurfaceLeftY;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Mz = -Ex(i,j,k);
				Update_Magnetic_Z(i,j-1,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Mx = Ez(i,j,k);
				Update_Magnetic_X(i,j-1,k,n);
			}
		}
	}

	//Right face of the virtual surface - Use Ex, Ez
	if ((jleft<=SurfaceRightY)&&(jright>=SurfaceRightY)) //is the right face of the box included in the node?
		//lowermost cell must be included in the node
		//- which cannot happen with more than one node
	{
		j=SurfaceRightY+1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Mz = Ex(i,j,k);
				Update_Magnetic_Z(i,j,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Mx = -Ez(i,j,k);
				Update_Magnetic_X(i,j,k,n);
			}
		}
	}
}

void Ctr_td_fs::UpdateElectric(const int& n)
{
	//Using the equivalence theorem, the H-values half-cell away from the virtual surface are converted to Jeq-values.
	//The Jeq values coincide with the E-field values on the virtual surface, and are used to update Electric_(X,Y,Z).
	//It must be ensured that each H-value updates Electric_(X,Y,Z) in only ONE NODE.

	//Bottom face of the virtual surface - Use Hx, Hy
	if ((klower<=SurfaceLowerZ-1)&&(kupper>=SurfaceLowerZ-1)) //is the bottom face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		k=SurfaceLowerZ-1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
			{
				Jy = -Hx(i,j,k);
				Update_Electric_Y(i,j,k+1,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
			{
				Jx = Hy(i,j,k);
				Update_Electric_X(i,j,k+1,n);
			}
		}
	}

	//Upper face of the virtual surface - Use Hx, Hy
	if ((klower<=SurfaceUpperZ+1)&&(kupper>=SurfaceUpperZ+1)) //is the upper face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		k=SurfaceUpperZ+1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
			{
				Jy = Hx(i,j,k);
				Update_Electric_Y(i,j,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
			{
				Jx = -Hy(i,j,k);
				Update_Electric_X(i,j,k,n);
			}
		}
	}

	//Back face of the virtual surface - Use Hy, Hz
	if ((iback<=SurfaceBackX-1)&&(ifront>=SurfaceBackX-1))	//is the back face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		i=SurfaceBackX-1;
		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Jz = -Hy(i,j,k);
				Update_Electric_Z(i+1,j,k,n);
			}
		}

		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Jy = Hz(i,j,k);
				Update_Electric_Y(i+1,j,k,n);
			}
		}
	}

	//Front face of the virtual surface - Use Hy, Hz
	if ((iback<=SurfaceFrontX+1)&&(ifront>=SurfaceFrontX+1)) //is the front face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		i=SurfaceFrontX+1;
		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY+1);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Jz = Hy(i,j,k);
				Update_Electric_Z(i,j,k,n);
			}
		}

		for (j=max(jleft,SurfaceLeftY);j<=min(jright,SurfaceRightY);j++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Jy = -Hz(i,j,k);
				Update_Electric_Y(i,j,k,n);
			}
		}
	}

	//Left face of the virtual surface - Use Hx, Hz
	if ((jleft<=SurfaceLeftY-1)&&(jright>=SurfaceLeftY-1))	//is the left face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		j=SurfaceLeftY-1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Jz = Hx(i,j,k);
				Update_Electric_Z(i,j+1,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Jx = -Hz(i,j,k);
				Update_Electric_X(i,j+1,k,n);
			}
		}
	}

	//Right face of the virtual surface - Use Hx, Hz
	if ((jleft<=SurfaceRightY+1)&&(jright>=SurfaceRightY+1)) //is the right face of the box included in the node?
		//(H-values that are half-cell away from the virtual surface must be included in the node - which cannot happen with more than one node)
	{
		j=SurfaceRightY+1;
		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX+1);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ);k++)
			{
				Jz = -Hx(i,j,k);
				Update_Electric_Z(i,j,k,n);
			}
		}

		for (i=max(iback,SurfaceBackX);i<=min(ifront,SurfaceFrontX);i++)
		{
			for (k=max(klower,SurfaceLowerZ);k<=min(kupper,SurfaceUpperZ+1);k++)
			{
				Jx = Hz(i,j,k);
				Update_Electric_X(i,j,k,n);
			}
		}
	}
}

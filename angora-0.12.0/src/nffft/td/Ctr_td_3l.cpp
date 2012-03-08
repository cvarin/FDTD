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

//Defines the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_3l" for 3-layered lossless media

#include "headers.h"

#include "Ctr_td_3l.h"

//for the definition of MaterialId
#include "material_id.h"

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern double dx,dt;
extern int NCELLS_Z,NPML;
extern int NSTEPS;

extern int OriginX,OriginY,OriginZ;

extern Array<double,3> Ex,Ey,Ez,Hx,Hy,Hz;

extern int number_of_layers;

extern Array<ElectricMaterialIndexType_X,1> LayerElectricMaterial_X;
extern Array<ElectricMaterialIndexType_Y,1> LayerElectricMaterial_Y;
extern Array<ElectricMaterialIndexType_Z,1> LayerElectricMaterial_Z;
extern Array<MagneticMaterialIndexType_X,1> LayerMagneticMaterial_X;
extern Array<MagneticMaterialIndexType_Y,1> LayerMagneticMaterial_Y;
extern Array<MagneticMaterialIndexType_Z,1> LayerMagneticMaterial_Z;

extern Array<int,1> LayerLowerZIndices;
extern Array<int,1> LayerThicknesses;
extern Array<bool,1> IsLayerGrounded;

extern Array<double,1> eps_x;	//this transformer works only for isotropic materials, may generalize later

extern double epsilon_r_upper,epsilon_r_lower;
extern double c_o;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int i,j,k;


Ctr_td_3l::Ctr_td_3l(const TrDataType_td& MyData, const string& FarFieldFileName, const int& Index)
		: Ctr_td(MyData,FarFieldFileName,Index)
{
	if (number_of_layers>3)
	{//if the # of layers is above 3, there is a developing error, so display error message and exit (better error-reporting scheme later!?)
		if (rank==0)
		{
			cout << "Error: Ctr_td_3l class cannot be used if number of layers is more than 3." << endl;
		}
#ifndef MPI_DISABLE
		MPI_Barrier(MPI_CartSubComm);
#endif
		exit(-1);
	}

	//take the position, thickness and permittivity of the slab from the layering info
	if (number_of_layers==3)
	{//if the # of layers is 3, the dielectric slab is the middle layer
		if (IsLayerGrounded(2))
		{//if the uppermost layer is grounded, then we basically have a 2-layered medium
			int SlabLayerIndex = 2; // dummy slab has the same properties as the uppermost layer
			SlabPos = LayerLowerZIndices(SlabLayerIndex);	//dummy slab of thickness 1 above the material interface: index of the highest (and lowest) cell in slab is the same as the lowest z-index cell in the upper layer
			SlabThickness = 1;	//dummy slab of thickness 1
			epsilon_r_slab = eps_x(LayerElectricMaterial_X(SlabLayerIndex)); /** isotropy assumed! **/	//permittivity of the upper layer
			grounded = true; //the uppermost layer is grounded
		}
		else
		{//uppermost layer is not grounded (but the middle layer still could be)
			int SlabLayerIndex = 1; // slab is the middle layer, with index 1
			SlabPos = LayerLowerZIndices(SlabLayerIndex+1)-1;	//index of the highest cell in the slab
			SlabThickness = LayerThicknesses(SlabLayerIndex);	//thickness of the middle layer
			epsilon_r_slab = eps_x(LayerElectricMaterial_X(SlabLayerIndex)); /** isotropy assumed! **/	//permittivity of the middle layer
			grounded = IsLayerGrounded(SlabLayerIndex);	//is the middle layer grounded?
		}
	}
	else if (number_of_layers==2)
	{//if the # of layers is 2, assign a "dummy" slab of thickness 1 above the material interface
		int SlabLayerIndex = 1; // dummy slab has the same properties as the upper layer
		SlabPos = LayerLowerZIndices(SlabLayerIndex);	//dummy slab of thickness 1 above the material interface: index of the highest (and lowest) cell in slab is the same as the lowest z-index cell in the upper layer
		SlabThickness = 1;	//dummy slab of thickness 1
		epsilon_r_slab = eps_x(LayerElectricMaterial_X(SlabLayerIndex)); /** isotropy assumed! **/	//permittivity of the upper layer
		grounded = IsLayerGrounded(SlabLayerIndex); //the upper layer may be grounded
	}
	else if (number_of_layers==1)
	{//if the # of layers is 1, assign a "dummy" slab of thickness 1 above the grid origin
		int SlabLayerIndex = 0; // dummy slab has the same properties as the homogeneous space
		SlabPos = OriginZ;	//dummy slab of thickness 1 above the grid origin:  index of the highest (and lowest) cell in slab is the same as OriginZ
		SlabThickness = 1; //dummy slab of thickness 1
		epsilon_r_slab = eps_x(LayerElectricMaterial_X(SlabLayerIndex)); /** isotropy assumed! **/	//permittivity of the homogeneous space
		grounded = false;	//the dummy slab cannot be grounded, since there is only one layer
	}

	//Origin of the far-field graphs (r is measured from this point)
	FarFieldOriginX = Data.NFFFTOriginX;
	FarFieldOriginY = Data.NFFFTOriginY;
	FarFieldOriginZ = SlabPos+1;

	//determine the observation half space
	//epsilon_o is the relative permittivity at that half space
	if (CosT>=0)
	{//observation half space is the uppermost half space
		epsilon_o = epsilon_r_upper;			//relative permittivity of the uppermost half space
	}
	else
	{//observation half space is the lowermost half space
		epsilon_o = epsilon_r_lower;			//relative permittivity of the lowermost half space
	}
	eps_r_0 = epsilon_r_upper/epsilon_o;	//relative permittivity of the uppermost half space, normalized by epsilon_o
	eps_r_1 = epsilon_r_slab/epsilon_o;		//relative permittivity of the slab, normalized by epsilon_o
	eps_r_2 = epsilon_r_lower/epsilon_o;	//relative permittivity of the lower half space, normalized by epsilon_o

	c_o = c/sqrt(epsilon_o);	//velocity of light in the observation half space
	Z_o = eta_0/sqrt(epsilon_o);		//wave impedance in the observation half space

	//******************************************************************************//
	//		BEGIN - TL GREEN FUNCTIONS FOR GROUNDED/UNGROUNDED DIELECTRIC SLAB		//
	//******************************************************************************//

	Calculate_TL_GreenFunctions();

	//******************************************************************************//
	//		END - TL GREEN FUNCTIONS FOR GROUNDED/UNGROUNDED DIELECTRIC SLAB		//
	//******************************************************************************//

	//the maximum delay (in time steps) due to reflections from the slab
	MaxDelay = max(V_e(V_e.size()-1).Delay,V_h(V_h.size()-1).Delay);

	//Max. advance in the xy-plane:
	MaxLateralAdvance = dx*sqrt(double((SurfaceFrontX-SurfaceBackX+1)*(SurfaceFrontX-SurfaceBackX+1)+(SurfaceRightY-SurfaceLeftY+1)*(SurfaceRightY-SurfaceLeftY+1)))/c_o/dt;
	//Max. advance in the z-direction:
	MaxVertAdvance = dx*abs(SurfaceUpperZ-SurfaceLowerZ+1)*sqrt(eps_r_1)/c_o/dt;
	//Max. total advance (may NEVER be attained, but included for all-positive delays):
	MaxTotalAdvance = MaxLateralAdvance+MaxVertAdvance;

	//Offset distance (to make all delays positive):
	r_offset = MaxTotalAdvance*dt*c_o;

	//Max. delay in the xy-plane:
	MaxLateralDelay = MaxLateralAdvance;
	//Max. delay in the z-direction:
	MaxVertDelay = MaxVertAdvance;
	//Max. total delay:
	MaxTotalDelay = MaxLateralDelay+MaxVertDelay+MaxDelay;

	ExtraSteps = int(MaxTotalDelay+MaxTotalAdvance)+1;	//Extra time steps needed for the far-field waveform

	TotalSteps = NSTEPS + ExtraSteps;	//Total number of time steps in the far-field waveform

	waveformTheta.resize(TotalSteps);
	waveformPhi.resize(TotalSteps);

	Electric_X_e.resize(TotalSteps);
	Electric_X_h.resize(TotalSteps);
	Electric_Y_e.resize(TotalSteps);
	Electric_Y_h.resize(TotalSteps);
	Electric_Z_e.resize(TotalSteps);
	Magnetic_X_e.resize(TotalSteps);
	Magnetic_X_h.resize(TotalSteps);
	Magnetic_Y_e.resize(TotalSteps);
	Magnetic_Y_h.resize(TotalSteps);
	Magnetic_Z_h.resize(TotalSteps);
	Theta_J.resize(TotalSteps);
	Phi_J.resize(TotalSteps);
	Phi_M.resize(TotalSteps);
	Theta_M.resize(TotalSteps);

	waveformTheta=0;
	waveformPhi=0;
	Electric_X_e=0;
	Electric_X_h=0;
	Electric_Y_e=0;
	Electric_Y_h=0;
	Electric_Z_e=0;
	Magnetic_X_e=0;
	Magnetic_X_h=0;
	Magnetic_Y_e=0;
	Magnetic_Y_h=0;
	Magnetic_Z_h=0;
	Theta_J=0;
	Phi_J=0;
	Phi_M=0;
	Theta_M=0;

	if (rank==0)
	{
		//first, record the package version
//		FarFieldFile.write((char*)&package_version,sizeof(package_version));
		FarFieldFile.write((char*)&angora_version_major,sizeof(angora_version_major));
		FarFieldFile.write((char*)&angora_version_minor,sizeof(angora_version_minor));
		FarFieldFile.write((char*)&angora_version_revision,sizeof(angora_version_revision));
		FarFieldFile.write((char*)&NSTEPS,sizeof(NSTEPS));
		FarFieldFile.write((char*)&TotalSteps,sizeof(TotalSteps));
		FarFieldFile.write((char*)&dt,sizeof(dt));
	}
}

void Ctr_td_3l::UpdateFarField(const int& n)
{
	 UpdateElectric(n);
	 UpdateMagnetic(n);
}

//The following code is purposefully repeated in other classes derived from Ctr_td, instead of using a common definition. The use of virtual functions for Update_Magnetic_Y, etc. would decrease the efficiency.
void Ctr_td_3l::UpdateMagnetic(const int& n)
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

void Ctr_td_3l::UpdateElectric(const int& n)
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

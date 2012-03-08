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

//Defines the PML object "CPML"

#include "headers.h"

//for definition of MaterialId
#include "material_id.h"

#include "CPML.h"

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern Array<ElectricMaterialIndexType_X,1> Layering_e_x;
extern Array<ElectricMaterialIndexType_Y,1> Layering_e_y;
extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
//According to magnetic properties:
extern Array<MagneticMaterialIndexType_X,1> Layering_h_x;
extern Array<MagneticMaterialIndexType_Y,1> Layering_h_y;
extern Array<MagneticMaterialIndexType_Z,1> Layering_h_z;

extern Array<double,1> Cb_X,Cb_Y,Cb_Z;
extern Array<double,1> Db_X,Db_Y,Db_Z;
extern Array<double,1> eps_x,eps_y,eps_z;

extern const int PEC;

extern int rank;
extern int rank_x, rank_y, rank_z;
extern int nodes_x, nodes_y, nodes_z;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int i,j,k;

extern Array<double,1> kappa_e_x,kappa_e_y,kappa_e_z,kappa_h_x,kappa_h_y,kappa_h_z;


Cpml::Cpml(const int& my_pml_thickness, const double& SizeOfScatterer)
		: pml_thickness(my_pml_thickness)
{
	if (pml_thickness>0)
	{
		//PML Parameters
		m=4;	//Order for polynomial grading
		epsilon_eff_x=0;	//Effective permittivity in the yz interface plane will be found by arithmetic weighting
		int usable_thickness = 0;	//usable thickness of the layered structure
		k = NCELLS_Z+2*pml_thickness;
		do
		{
			usable_thickness++;
			epsilon_eff_x += eps_z(Layering_e_z(k));  /** isotropy assumed **/
			k--;
		} while (k>=1); //terminate at k=1
		epsilon_eff_x = epsilon_eff_x/usable_thickness;

		epsilon_eff_y=epsilon_eff_x; /** isotropy assumed **/	//Same as epsilon_eff_x, from lateral symmetry
		epsilon_eff_lowerz=eps_x(Layering_e_x(pml_thickness+1)); /** isotropy assumed **/  //Effective permittivity in the lower xy interface plane
		epsilon_eff_upperz=eps_x(Layering_e_x(NCELLS_Z+pml_thickness+1)); /** isotropy assumed **/ //Effective permittivity in the upper xy interface plane
		sigma_coeff=1;		//sigma_max/sigma_optimum (sigma_opt as in Taflove 7.61)

		//Maximum sigma values in (x,y,z)-normal PML permittivity tensor (max. in back plane)
		sigma_x_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_x,0.5)*dx)*sigma_coeff;
		sigma_y_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_y,0.5)*dx)*sigma_coeff;
		sigma_lowerz_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_lowerz,0.5)*dx)*sigma_coeff;
		sigma_upperz_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_upperz,0.5)*dx)*sigma_coeff;

		//Maximum kappa values in (x,y,z)-normal PML permittivity tensor (max. in back plane)
		kappa_x_max=1;
		kappa_y_max=1;
		kappa_z_max=1;

		//Maximum alpha in (x,y,z)-normal PML permittivity tensor (max. in front plane)
		//alpha_opt = c*eps/w , where w is the size of the scatterer
		// (see Berenger, "Numerical reflection from FDTD-PMLs - a comparison of the split PML with the unsplit and CFS PMLs", IEEE Ant.Prop. vol 50, Mar. 2002)
		alpha_x_max=c*epsilon_0/(dx*SizeOfScatterer);
		alpha_y_max=c*epsilon_0/(dx*SizeOfScatterer);
		alpha_z_max=c*epsilon_0/(dx*SizeOfScatterer);
	// 	alpha_x_max=1e-7;
	// 	alpha_y_max=1e-7;
	// 	alpha_z_max=1e-7;

		//Minimum alpha in (x,y,z)-normal PML permittivity tensor (min. in back plane)
		alpha_x_min=0;
		alpha_y_min=0;
		alpha_z_min=0;

		//Allocate and initialize PML current arrays
		if (rank_x==0)
		{
			Psi_Eyx_BackX.resize(Range(1,pml_thickness+1),Range(jleft,jright),Range(klower,kupper+1));
			Psi_Ezx_BackX.resize(Range(1,pml_thickness+1),Range(jleft,jright+1),Range(klower,kupper));
			Psi_Hyx_BackX.resize(Range(1,pml_thickness),Range(jleft,jright+1),Range(klower,kupper));
			Psi_Hzx_BackX.resize(Range(1,pml_thickness),Range(jleft,jright),Range(klower,kupper+1));
			Psi_Eyx_BackX=0;
			Psi_Ezx_BackX=0;
			Psi_Hyx_BackX=0;
			Psi_Hzx_BackX=0;
		}

		if (rank_x==(nodes_x-1))
		{
			Psi_Eyx_FrontX.resize(Range(1,pml_thickness+1),Range(jleft,jright),Range(klower,kupper+1));
			Psi_Ezx_FrontX.resize(Range(1,pml_thickness+1),Range(jleft,jright+1),Range(klower,kupper));
			Psi_Hyx_FrontX.resize(Range(1,pml_thickness),Range(jleft,jright+1),Range(klower,kupper));
			Psi_Hzx_FrontX.resize(Range(1,pml_thickness),Range(jleft,jright),Range(klower,kupper+1));
			Psi_Eyx_FrontX=0;
			Psi_Ezx_FrontX=0;
			Psi_Hyx_FrontX=0;
			Psi_Hzx_FrontX=0;
		}

		if (rank_y==0)
		{
			Psi_Exy_LeftY.resize(Range(iback,ifront),Range(1,pml_thickness+1),Range(klower,kupper+1));
			Psi_Ezy_LeftY.resize(Range(iback,ifront+1),Range(1,pml_thickness+1),Range(klower,kupper));
			Psi_Hxy_LeftY.resize(Range(iback,ifront+1),Range(1,pml_thickness),Range(klower,kupper));
			Psi_Hzy_LeftY.resize(Range(iback,ifront),Range(1,pml_thickness),Range(klower,kupper+1));
			Psi_Exy_LeftY=0;
			Psi_Ezy_LeftY=0;
			Psi_Hxy_LeftY=0;
			Psi_Hzy_LeftY=0;
		}

		if (rank_y==(nodes_y-1))
		{
			Psi_Exy_RightY.resize(Range(iback,ifront),Range(1,pml_thickness+1),Range(klower,kupper+1));
			Psi_Ezy_RightY.resize(Range(iback,ifront+1),Range(1,pml_thickness+1),Range(klower,kupper));
			Psi_Hxy_RightY.resize(Range(iback,ifront+1),Range(1,pml_thickness),Range(klower,kupper));
			Psi_Hzy_RightY.resize(Range(iback,ifront),Range(1,pml_thickness),Range(klower,kupper+1));
			Psi_Exy_RightY=0;
			Psi_Ezy_RightY=0;
			Psi_Hxy_RightY=0;
			Psi_Hzy_RightY=0;
		}

		if (rank_z==0)
		{
			Psi_Exz_LowerZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(1,pml_thickness+1));
			Psi_Eyz_LowerZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(1,pml_thickness+1));
			Psi_Hxz_LowerZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(1,pml_thickness));
			Psi_Hyz_LowerZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(1,pml_thickness));
			Psi_Exz_LowerZ=0;
			Psi_Eyz_LowerZ=0;
			Psi_Hxz_LowerZ=0;
			Psi_Hyz_LowerZ=0;
		}

		if (rank_z==(nodes_z-1))
		{
			Psi_Exz_UpperZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(1,pml_thickness+1));
			Psi_Eyz_UpperZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(1,pml_thickness+1));
			Psi_Hxz_UpperZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(1,pml_thickness));
			Psi_Hyz_UpperZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(1,pml_thickness));
			Psi_Exz_UpperZ=0;
			Psi_Eyz_UpperZ=0;
			Psi_Hxz_UpperZ=0;
			Psi_Hyz_UpperZ=0;
		}

		//Modify the coordinate-stretching kappa parameters in the PML region
		//BackX and FrontX
		for (int i=1;i<=pml_thickness+1;i++)
		{
			//BackX
			kappa_e_x(i)=1+(kappa_x_max-1)*pow((pml_thickness-i+1.0)/pml_thickness,m);
			//FrontX
			kappa_e_x(NCELLS_X+2*pml_thickness+2-i)=kappa_e_x(i);
		}

		for (int i=1;i<=pml_thickness;i++)
		{
			//BackX
			kappa_h_x(i)=1+(kappa_x_max-1)*pow((pml_thickness-i+0.5)/pml_thickness,m);
			//FrontX
			kappa_h_x(NCELLS_X+2*pml_thickness+1-i)=kappa_h_x(i);
		}

		//LeftY and RightY
		for (int j=1;j<=pml_thickness+1;j++)
		{
			//LeftY
			kappa_e_y(j)=1+(kappa_y_max-1)*pow((pml_thickness-j+1.0)/pml_thickness,m);
			//RightY
			kappa_e_y(NCELLS_Y+2*pml_thickness+2-j)=kappa_e_y(j);
		}
		for (int j=1;j<=pml_thickness;j++)
		{
			//LeftY
			kappa_h_y(j)=1+(kappa_y_max-1)*pow((pml_thickness-j+0.5)/pml_thickness,m);
			//RightY
			kappa_h_y(NCELLS_Y+2*pml_thickness+1-j)=kappa_h_y(j);
		}

		//UpperZ and LowerZ
		for (int k=1;k<=pml_thickness+1;k++)
		{
			//LowerZ
			kappa_e_z(k)=1+(kappa_z_max-1)*pow((pml_thickness-k+1.0)/pml_thickness,m);
			//UpperZ
			kappa_e_z(NCELLS_Z+2*pml_thickness+2-k)=kappa_e_z(k);
		}
		for (int k=1;k<=pml_thickness;k++)
		{
			//LowerZ
			kappa_h_z(k)=1+(kappa_z_max-1)*pow((pml_thickness-k+0.5)/pml_thickness,m);
			//UpperZ
			kappa_h_z(NCELLS_Z+2*pml_thickness+1-k)=kappa_h_z(k);
		}

		//Allocate and initialize PML current update parameters
		//	for updating the E-fields perpendicular to stretching axis:
		//Evaluated at the edge centers
		Apml_e_backx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the back PML slab
		Apml_e_frontx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the front PML slab
		Apml_e_lefty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the left PML slab
		Apml_e_righty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the right PML slab
		Apml_e_lowerz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the lower PML slab
		Apml_e_upperz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the upper PML slab
		Bpml_e_backx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the back PML slab
		Bpml_e_frontx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the front PML slab
		Bpml_e_lefty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the left PML slab
		Bpml_e_righty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the right PML slab
		Bpml_e_lowerz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the lower PML slab
		Bpml_e_upperz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the upper PML slab
		//	for updating the H-fields perpendicular to stretching axis:
		//Evaluated at the face centers
		Apml_h_backx.resize(Range(1,pml_thickness));		//for Hy and Hz in the back PML slab
		Apml_h_frontx.resize(Range(1,pml_thickness));	//for Hy and Hz in the front PML slab
		Apml_h_lefty.resize(Range(1,pml_thickness));		//for Hx and Hz in the left PML slab
		Apml_h_righty.resize(Range(1,pml_thickness));	//for Hx and Hz in the right PML slab
		Apml_h_lowerz.resize(Range(1,pml_thickness));	//for Hx and Hy in the lower PML slab
		Apml_h_upperz.resize(Range(1,pml_thickness));	//for Hx and Hy in the upper PML slab
		Bpml_h_backx.resize(Range(1,pml_thickness));		//for Hy and Hz in the back PML slab
		Bpml_h_frontx.resize(Range(1,pml_thickness));	//for Hy and Hz in the front PML slab
		Bpml_h_lefty.resize(Range(1,pml_thickness));		//for Hx and Hz in the left PML slab
		Bpml_h_righty.resize(Range(1,pml_thickness));	//for Hx and Hz in the right PML slab
		Bpml_h_lowerz.resize(Range(1,pml_thickness));	//for Hx and Hy in the lower PML slab
		Bpml_h_upperz.resize(Range(1,pml_thickness));	//for Hx and Hy in the upper PML slab
		//Initialize PML current update parameters (1-D)
		//	for the E-field:
		//Back and Front PML slab
		for (int i=1;i<=pml_thickness+1;i++)
		{
			double sigma_x,kappa_x,alpha_x;
			double r = pow((i-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_x=r*sigma_x_max;
			kappa_x=1+(kappa_x_max-1)*r;
			alpha_x=alpha_x_max+(alpha_x_min-alpha_x_max)*r;
			//BackX
			Apml_e_backx(pml_thickness+2-i)=sigma_x/(sigma_x*kappa_x+kappa_x*kappa_x*alpha_x)*
				(exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0)-1.0);
			Bpml_e_backx(pml_thickness+2-i)=exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0);
			//FrontX
			Apml_e_frontx(i)=Apml_e_backx(pml_thickness+2-i);
			Bpml_e_frontx(i)=Bpml_e_backx(pml_thickness+2-i);
		}
		//Left and Right PML slab
		for (int j=1;j<=pml_thickness+1;j++)
		{
			double sigma_y,kappa_y,alpha_y;
			double r = pow((j-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_y=r*sigma_y_max;
			kappa_y=1+(kappa_y_max-1)*r;
			alpha_y=alpha_y_max+(alpha_y_min-alpha_y_max)*r;
			//Left
			Apml_e_lefty(pml_thickness+2-j)=sigma_y/(sigma_y*kappa_y+kappa_y*kappa_y*alpha_y)*
				(exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0)-1.0);
			Bpml_e_lefty(pml_thickness+2-j)=exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0);
			//Right
			Apml_e_righty(j)=Apml_e_lefty(pml_thickness+2-j);
			Bpml_e_righty(j)=Bpml_e_lefty(pml_thickness+2-j);
		}
		//Lower and Upper PML slab
		for (int k=1;k<=pml_thickness+1;k++)
		{
			double sigma_z,kappa_z,alpha_z;
			double r = pow((k-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_z=r*sigma_lowerz_max;
			kappa_z=1+(kappa_z_max-1)*r;
			alpha_z=alpha_z_max+(alpha_z_min-alpha_z_max)*r;
			//Lower
			Apml_e_lowerz(pml_thickness+2-k)=sigma_z/(sigma_z*kappa_z+kappa_z*kappa_z*alpha_z)*
				(exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0)-1.0);
			Bpml_e_lowerz(pml_thickness+2-k)=exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0);
			//Upper
			Apml_e_upperz(k)=Apml_e_lowerz(pml_thickness+2-k);
			Bpml_e_upperz(k)=Bpml_e_lowerz(pml_thickness+2-k);
		}

		//	for the H-field:
		//Back and Front PML slab
		for (int i=1;i<=pml_thickness;i++)
		{
			double sigma_x,kappa_x,alpha_x;
			double r = pow((i-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_x=r*sigma_x_max;
			kappa_x=1+(kappa_x_max-1)*r;
			alpha_x=alpha_x_max*(1-r);
			//BackX
			Apml_h_backx(pml_thickness+1-i)=sigma_x/(sigma_x*kappa_x+kappa_x*kappa_x*alpha_x)*
				(exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0)-1.0);
			Bpml_h_backx(pml_thickness+1-i)=exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0);
			//FrontX
			Apml_h_frontx(i)=Apml_h_backx(pml_thickness+1-i);
			Bpml_h_frontx(i)=Bpml_h_backx(pml_thickness+1-i);
		}
		//Left and Right PML slab
		for (int j=1;j<=pml_thickness;j++)
		{
			double sigma_y,kappa_y,alpha_y;
			double r = pow((j-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_y=r*sigma_y_max;
			kappa_y=1+(kappa_y_max-1)*r;
			alpha_y=alpha_y_max*(1-r);
			//Left
			Apml_h_lefty(pml_thickness+1-j)=sigma_y/(sigma_y*kappa_y+kappa_y*kappa_y*alpha_y)*
				(exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0)-1.0);
			Bpml_h_lefty(pml_thickness+1-j)=exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0);
			//Right
			Apml_h_righty(j)=Apml_h_lefty(pml_thickness+1-j);
			Bpml_h_righty(j)=Bpml_h_lefty(pml_thickness+1-j);
		}
		//Lower and Upper PML slab
		for (int k=1;k<=pml_thickness;k++)
		{
			double sigma_z,kappa_z,alpha_z;
			double r = pow((k-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_z=r*sigma_lowerz_max;
			kappa_z=1+(kappa_z_max-1)*r;
			alpha_z=alpha_z_max*(1-r);
			//Lower
			Apml_h_lowerz(pml_thickness+1-k)=sigma_z/(sigma_z*kappa_z+kappa_z*kappa_z*alpha_z)*
				(exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0)-1.0);
			Bpml_h_lowerz(pml_thickness+1-k)=exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0);
			//Upper
			Apml_h_upperz(k)=Apml_h_lowerz(pml_thickness+1-k);
			Bpml_h_upperz(k)=Bpml_h_lowerz(pml_thickness+1-k);
		}
	}
}

void Cpml::UpdateE()
{
	if (pml_thickness>0)
	{
		UpdateE_X();
		UpdateE_Y();
		UpdateE_Z();
	}
}

void Cpml::UpdateH()
{
	if (pml_thickness>0)
	{
		UpdateH_X();
		UpdateH_Y();
		UpdateH_Z();
	}
}

void Cpml::UpdateE_X()
{
	if (rank_x==0)
	{
		//Psi_Eyx_BackX
		for (i=1; i<=pml_thickness; i++){
			for (j=jleft; j<=jright; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//BackX
					Psi_Eyx_BackX(i+1,j,k)=Bpml_e_backx(i+1)*Psi_Eyx_BackX(i+1,j,k)
						+ Apml_e_backx(i+1)*
						(Hz(i+1,j,k)-Hz(i,j,k))/dx;
					Ey(i+1,j,k) +=
						- Cb_Y(Media_Ey(i+1,j,k))*dx*Psi_Eyx_BackX(i+1,j,k);
				}
			}
		}
		//Psi_Ezx_BackX
		for (i=1; i<=pml_thickness; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=klower; k<=kupper; k++)
				{
					//BackX
					Psi_Ezx_BackX(i+1,j,k)=Bpml_e_backx(i+1)*Psi_Ezx_BackX(i+1,j,k)
						+ Apml_e_backx(i+1)*
						(Hy(i+1,j,k)-Hy(i,j,k))/dx;
					Ez(i+1,j,k) +=
						+ Cb_Z(Media_Ez(i+1,j,k))*dx*Psi_Ezx_BackX(i+1,j,k);
				}
			}
		}
	}
	if (rank_x==(nodes_x-1))
	{
		//Psi_Eyx_FrontX
		for (i=1; i<=pml_thickness; i++){
			for (j=jleft; j<=jright; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//FrontX
					Psi_Eyx_FrontX(i,j,k)=Bpml_e_frontx(i)*Psi_Eyx_FrontX(i,j,k)
						+ Apml_e_frontx(i)*
						(Hz(i+NCELLS_X+pml_thickness,j,k)-Hz(i+NCELLS_X+pml_thickness-1,j,k))/dx;
					Ey(i+NCELLS_X+pml_thickness,j,k) +=
						- Cb_Y(Media_Ey(i+NCELLS_X+pml_thickness,j,k))*dx*Psi_Eyx_FrontX(i,j,k);
				}
			}
		}
		//Psi_Ezx_FrontX
		for (i=1; i<=pml_thickness; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=klower; k<=kupper; k++)
				{
					//FrontX
					Psi_Ezx_FrontX(i,j,k)=Bpml_e_frontx(i)*Psi_Ezx_FrontX(i,j,k)
						+ Apml_e_frontx(i)*
						(Hy(i+NCELLS_X+pml_thickness,j,k)-Hy(i+NCELLS_X+pml_thickness-1,j,k))/dx;
					Ez(i+NCELLS_X+pml_thickness,j,k) +=
						+ Cb_Z(Media_Ez(i+NCELLS_X+pml_thickness,j,k))*dx*Psi_Ezx_FrontX(i,j,k);
				}
			}
		}
	}
}

void Cpml::UpdateH_X()
{
	if (rank_x==0)
	{
		//Psi_Hyx_BackX
		for (i=1; i<=pml_thickness; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=klower; k<=kupper; k++)
				{
					//BackX
					Psi_Hyx_BackX(i,j,k)=Bpml_h_backx(i)*Psi_Hyx_BackX(i,j,k)
						+ Apml_h_backx(i)*
						(Ez(i+1,j,k)-Ez(i,j,k))/dx;
					Hy(i,j,k) +=
						+ Db_Y(Media_Hy(i,j,k))*dx*Psi_Hyx_BackX(i,j,k);
				}
			}
		}
		//Psi_Hzx_BackX
		for (i=1; i<=pml_thickness; i++){
			for (j=jleft; j<=jright; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//BackX
					Psi_Hzx_BackX(i,j,k)=Bpml_h_backx(i)*Psi_Hzx_BackX(i,j,k)
						+ Apml_h_backx(i)*
						(Ey(i+1,j,k)-Ey(i,j,k))/dx;
					Hz(i,j,k) +=
						- Db_Z(Media_Hz(i,j,k))*dx*Psi_Hzx_BackX(i,j,k);
				}
			}
		}
	}
	if (rank_x==(nodes_x-1))
	{
		//Psi_Hyx_FrontX
		for (i=1; i<=pml_thickness; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=klower; k<=kupper; k++)
				{
					//FrontX
					Psi_Hyx_FrontX(i,j,k)=Bpml_h_frontx(i)*Psi_Hyx_FrontX(i,j,k)
						+ Apml_h_frontx(i)*
						(Ez(i+NCELLS_X+pml_thickness+1,j,k)-Ez(i+NCELLS_X+pml_thickness,j,k))/dx;
					Hy(i+NCELLS_X+pml_thickness,j,k) +=
						+ Db_Y(Media_Hy(i+NCELLS_X+pml_thickness,j,k))*dx*Psi_Hyx_FrontX(i,j,k);
				}
			}
		}
		//Psi_Hzx_FrontX
		for (i=1; i<=pml_thickness; i++){
			for (j=jleft; j<=jright; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//FrontX
					Psi_Hzx_FrontX(i,j,k)=Bpml_h_frontx(i)*Psi_Hzx_FrontX(i,j,k)
						+ Apml_h_frontx(i)*
						(Ey(i+NCELLS_X+pml_thickness+1,j,k)-Ey(i+NCELLS_X+pml_thickness,j,k))/dx;
					Hz(i+NCELLS_X+pml_thickness,j,k) +=
						- Db_Z(Media_Hz(i+NCELLS_X+pml_thickness,j,k))*dx*Psi_Hzx_FrontX(i,j,k);
				}
			}
		}
	}
}

void Cpml::UpdateE_Y()
{
	if (rank_y==0)
	{
		//Psi_Exy_LeftY
		for (i=iback; i<=ifront; i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//LeftY
					Psi_Exy_LeftY(i,j+1,k)=Bpml_e_lefty(j+1)*Psi_Exy_LeftY(i,j+1,k)
						+ Apml_e_lefty(j+1)*
						(Hz(i,j+1,k)-Hz(i,j,k))/dx;
					Ex(i,j+1,k) +=
						+ Cb_X(Media_Ex(i,j+1,k))*dx*Psi_Exy_LeftY(i,j+1,k);
				}
			}
		}
		//Psi_Ezy_LeftY
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=klower; k<=kupper; k++)
				{
					//LeftY
					Psi_Ezy_LeftY(i,j+1,k)=Bpml_e_lefty(j+1)*Psi_Ezy_LeftY(i,j+1,k)
						+ Apml_e_lefty(j+1)*
						(Hx(i,j+1,k)-Hx(i,j,k))/dx;
					Ez(i,j+1,k) +=
						- Cb_Z(Media_Ez(i,j+1,k))*dx*Psi_Ezy_LeftY(i,j+1,k);
				}
			}
		}
	}
	if (rank_y==(nodes_y-1))
	{
		//Psi_Exy_RightY
		for (i=iback; i<=ifront; i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//RightY
					Psi_Exy_RightY(i,j,k)=Bpml_e_righty(j)*Psi_Exy_RightY(i,j,k)
						+ Apml_e_righty(j)*
						(Hz(i,j+NCELLS_Y+pml_thickness,k)-Hz(i,j+NCELLS_Y+pml_thickness-1,k))/dx;
					Ex(i,j+NCELLS_Y+pml_thickness,k) +=
						+ Cb_X(Media_Ex(i,j+NCELLS_Y+pml_thickness,k))*dx*Psi_Exy_RightY(i,j,k);
				}
			}
		}
		//Psi_Ezy_RightY
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=klower; k<=kupper; k++)
				{
					//RightY
					Psi_Ezy_RightY(i,j,k)=Bpml_e_righty(j)*Psi_Ezy_RightY(i,j,k)
						+ Apml_e_righty(j)*
						(Hx(i,j+NCELLS_Y+pml_thickness,k)-Hx(i,j+NCELLS_Y+pml_thickness-1,k))/dx;
					Ez(i,j+NCELLS_Y+pml_thickness,k) +=
						- Cb_Z(Media_Ez(i,j+NCELLS_Y+pml_thickness,k))*dx*Psi_Ezy_RightY(i,j,k);
				}
			}
		}
	}
}

void Cpml::UpdateH_Y()
{
	if (rank_y==0)
	{
		//Psi_Hxy_LeftY
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=klower; k<=kupper; k++)
				{
					//LeftY
					Psi_Hxy_LeftY(i,j,k)=Bpml_h_lefty(j)*Psi_Hxy_LeftY(i,j,k)
						+ Apml_h_lefty(j)*
						(Ez(i,j+1,k)-Ez(i,j,k))/dx;
					Hx(i,j,k) +=
						- Db_X(Media_Hx(i,j,k))*dx*Psi_Hxy_LeftY(i,j,k);
				}
			}
		}
		//Psi_Hzy_LeftY
		for (i=iback; i<=ifront; i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//LeftY
					Psi_Hzy_LeftY(i,j,k)=Bpml_h_lefty(j)*Psi_Hzy_LeftY(i,j,k)
						+ Apml_h_lefty(j)*
						(Ex(i,j+1,k)-Ex(i,j,k))/dx;
					Hz(i,j,k) +=
						+ Db_Z(Media_Hz(i,j,k))*dx*Psi_Hzy_LeftY(i,j,k);
				}
			}
		}
	}
	if (rank_y==(nodes_y-1))
	{
		//Psi_Hxy_RightY
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=klower; k<=kupper; k++)
				{
					//RightY
					Psi_Hxy_RightY(i,j,k)=Bpml_h_righty(j)*Psi_Hxy_RightY(i,j,k)
						+ Apml_h_righty(j)*
						(Ez(i,j+NCELLS_Y+pml_thickness+1,k)-Ez(i,j+NCELLS_Y+pml_thickness,k))/dx;
					Hx(i,j+NCELLS_Y+pml_thickness,k) +=
						- Db_X(Media_Hx(i,j+NCELLS_Y+pml_thickness,k))*dx*Psi_Hxy_RightY(i,j,k);
				}
			}
		}
		//Psi_Hzy_RightY
		for (i=iback; i<=ifront; i++){
			for (j=1; j<=pml_thickness; j++){
				for (k=max(2,klower); k<=min(kupper+1,NCELLS_Z+2*pml_thickness); k++)
				{
					//RightY
					Psi_Hzy_RightY(i,j,k)=Bpml_h_righty(j)*Psi_Hzy_RightY(i,j,k)
						+ Apml_h_righty(j)*
						(Ex(i,j+NCELLS_Y+pml_thickness+1,k)-Ex(i,j+NCELLS_Y+pml_thickness,k))/dx;
					Hz(i,j+NCELLS_Y+pml_thickness,k) +=
						+ Db_Z(Media_Hz(i,j+NCELLS_Y+pml_thickness,k))*dx*Psi_Hzy_RightY(i,j,k);
				}
			}
		}
	}
}

void Cpml::UpdateE_Z()
{
	if (rank_z==0)
	{
		//Psi_Exz_LowerZ
		for (i=iback; i<=ifront; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//LowerZ
					Psi_Exz_LowerZ(i,j,k+1)=Bpml_e_lowerz(k+1)*Psi_Exz_LowerZ(i,j,k+1)
						+ Apml_e_lowerz(k+1)*
						(Hy(i,j,k+1)-Hy(i,j,k))/dx;
					Ex(i,j,k+1) = Ex(i,j,k+1)
						- Cb_X(Media_Ex(i,j,k+1))*dx*Psi_Exz_LowerZ(i,j,k+1);
				}
			}
		}
		//Psi_Eyz_LowerZ
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=jleft; j<=jright; j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//LowerZ
					Psi_Eyz_LowerZ(i,j,k+1)=Bpml_e_lowerz(k+1)*Psi_Eyz_LowerZ(i,j,k+1)
						+ Apml_e_lowerz(k+1)*
						(Hx(i,j,k+1)-Hx(i,j,k))/dx;
					Ey(i,j,k+1) = Ey(i,j,k+1)
						+ Cb_Y(Media_Ey(i,j,k+1))*dx*Psi_Eyz_LowerZ(i,j,k+1);
				}
			}
		}
	}

	if (rank_z==(nodes_z-1))
	{
		//Psi_Exz_UpperZ
		for (i=iback; i<=ifront; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//UpperZ
					Psi_Exz_UpperZ(i,j,k)=Bpml_e_upperz(k)*Psi_Exz_UpperZ(i,j,k)
						+ Apml_e_upperz(k)*
						(Hy(i,j,k+NCELLS_Z+pml_thickness)-Hy(i,j,k+NCELLS_Z+pml_thickness-1))/dx;
					Ex(i,j,k+NCELLS_Z+pml_thickness) +=
						- Cb_X(Media_Ex(i,j,k+NCELLS_Z+pml_thickness))*dx*Psi_Exz_UpperZ(i,j,k);
				}
			}
		}
		//Psi_Eyz_UpperZ
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=jleft; j<=jright; j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//UpperZ
					Psi_Eyz_UpperZ(i,j,k)=Bpml_e_upperz(k)*Psi_Eyz_UpperZ(i,j,k)
						+ Apml_e_upperz(k)*
						(Hx(i,j,k+NCELLS_Z+pml_thickness)-Hx(i,j,k+NCELLS_Z+pml_thickness-1))/dx;
					Ey(i,j,k+NCELLS_Z+pml_thickness) +=
						+ Cb_Y(Media_Ey(i,j,k+NCELLS_Z+pml_thickness))*dx*Psi_Eyz_UpperZ(i,j,k);
				}
			}
		}
	}
}

void Cpml::UpdateH_Z()
{
	if (rank_z==0)
	{
		//Psi_Hxz_LowerZ
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=jleft; j<=jright; j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//LowerZ
					Psi_Hxz_LowerZ(i,j,k)=Bpml_h_lowerz(k)*Psi_Hxz_LowerZ(i,j,k)
						+ Apml_h_lowerz(k)*
						(Ey(i,j,k+1)-Ey(i,j,k))/dx;
					Hx(i,j,k) = Hx(i,j,k)
						+ Db_X(Media_Hx(i,j,k))*dx*Psi_Hxz_LowerZ(i,j,k);
				}
			}
		}
		//Psi_Hyz_LowerZ
		for (i=iback; i<=ifront; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//LowerZ
					Psi_Hyz_LowerZ(i,j,k)=Bpml_h_lowerz(k)*Psi_Hyz_LowerZ(i,j,k)
						+ Apml_h_lowerz(k)*
						(Ex(i,j,k+1)-Ex(i,j,k))/dx;
					Hy(i,j,k) = Hy(i,j,k)
						- Db_Y(Media_Hy(i,j,k))*dx*Psi_Hyz_LowerZ(i,j,k);
				}
			}
		}
	}
	if (rank_z==(nodes_z-1))
	{
		//Psi_Hxz_UpperZ
		for (i=max(2,iback); i<=min(ifront+1,NCELLS_X+2*pml_thickness); i++){
			for (j=jleft; j<=jright; j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//UpperZ
					Psi_Hxz_UpperZ(i,j,k)=Bpml_h_upperz(k)*Psi_Hxz_UpperZ(i,j,k)
						+ Apml_h_upperz(k)*
						(Ey(i,j,k+NCELLS_Z+pml_thickness+1)-Ey(i,j,k+NCELLS_Z+pml_thickness))/dx;
					Hx(i,j,k+NCELLS_Z+pml_thickness) +=
						+ Db_X(Media_Hx(i,j,k+NCELLS_Z+pml_thickness))*dx*Psi_Hxz_UpperZ(i,j,k);
				}
			}
		}
		//Psi_Hyz_UpperZ
		for (i=iback; i<=ifront; i++){
			for (j=max(2,jleft); j<=min(jright+1,NCELLS_Y+2*pml_thickness); j++){
				for (k=1; k<=pml_thickness; k++)
				{
					//UpperZ
					Psi_Hyz_UpperZ(i,j,k)=Bpml_h_upperz(k)*Psi_Hyz_UpperZ(i,j,k)
						+ Apml_h_upperz(k)*
						(Ex(i,j,k+NCELLS_Z+pml_thickness+1)-Ex(i,j,k+NCELLS_Z+pml_thickness))/dx;
					Hy(i,j,k+NCELLS_Z+pml_thickness) +=
						- Db_Y(Media_Hy(i,j,k+NCELLS_Z+pml_thickness))*dx*Psi_Hyz_UpperZ(i,j,k);
				}
			}
		}
	}
}

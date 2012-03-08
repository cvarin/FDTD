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

//Defines functions that initialize the geometry in the grid.

#include "headers.h"

#include "initgeom.h"

#include "time_axis.h"

//for definition of MaterialId
#include "material_id.h"

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern Array<ElectricMaterialIndexType_X,3> Media_Ex;
extern Array<ElectricMaterialIndexType_Y,3> Media_Ey;
extern Array<ElectricMaterialIndexType_Z,3> Media_Ez;
extern Array<MagneticMaterialIndexType_X,3> Media_Hx;
extern Array<MagneticMaterialIndexType_Y,3> Media_Hy;
extern Array<MagneticMaterialIndexType_Z,3> Media_Hz;

extern int NumOfElectricMaterials_X;
extern int NumOfElectricMaterials_Y;
extern int NumOfElectricMaterials_Z;
extern int NumOfMagneticMaterials_X;
extern int NumOfMagneticMaterials_Y;
extern int NumOfMagneticMaterials_Z;
extern const int PEC;
extern const int vacuum;

extern Array<ElectricMaterialIndexType_X,1> Layering_e_x;
extern Array<ElectricMaterialIndexType_Y,1> Layering_e_y;
extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;
extern Array<MagneticMaterialIndexType_X,1> Layering_h_x;
extern Array<MagneticMaterialIndexType_Y,1> Layering_h_y;
extern Array<MagneticMaterialIndexType_Z,1> Layering_h_z;

extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<double,1> eps_x,eps_y,eps_z;
extern Array<double,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<double,1> mu_x,mu_y,mu_z;
extern Array<double,1> cond_h_x,cond_h_y,cond_h_z;
extern Array<double,1> kappa_e_x,kappa_e_y,kappa_e_z,kappa_h_x,kappa_h_y,kappa_h_z;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;
extern double c_upper,c_lower;

extern double max_field_value;
extern bool max_field_value_set_in_configfile;
extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;


void init_geom()
{//resets the material information and assigns vacuum to the entire grid
	//the following is not automatic, since named materials are handled in the object Cmats
//	/** FIXME: This is not the most ideal solution **/
//	//initialize the array of named (tagged) materials to zero length
//	InitializeNamedMaterialArray();
//	/** FIXME: This is not the most ideal solution **/

	//set the initial time value to 0
//	allow_setting_initial_time_value();
	set_initial_time_value(0);

	//Initialize some global variables
	//initially, there are only 2 materials: vacuum and PEC
	NumOfElectricMaterials_X = 2;
	NumOfElectricMaterials_Y = 2;
	NumOfElectricMaterials_Z = 2;
	NumOfMagneticMaterials_X = 2;
	NumOfMagneticMaterials_Y = 2;
	NumOfMagneticMaterials_Z = 2;

	epsilon_r_upper=1; 	//relative permittivity of the uppermost layer (changed only through PlaceSlab in geometry.cpp)
	epsilon_r_lower=1; 	//relative permittivity of the lowermost layer (changed only through PlaceSlab in geometry.cpp)
	c_upper=299792458.;	//velocity of propagation in the uppermost layer (changed only through PlaceSlab in geometry.cpp)
	c_lower=299792458.;	//velocity of propagation in the lowermost layer (changed only through PlaceSlab in geometry.cpp)

	//set the maximum electric field value encountered in grid to 1, if not set in the config file already
	if (!max_field_value_set_in_configfile) max_field_value = 1;

	if (!dB_accuracy_set_in_configfile) dB_accuracy = -60;		//useful accuracy range of the grid (in -dB), default to -60 dB

	//Initialize field arrays
	Ex=0;
	Ey=0;
	Ez=0;
	Hx=0;
	Hy=0;
	Hz=0;
	//Place materials
	Media_Ex=vacuum;
	Media_Ey=vacuum;
	Media_Ez=vacuum;
	Media_Hx=vacuum;
	Media_Hy=vacuum;
	Media_Hz=vacuum;

	//Allocate and initialize constitutive parameter arrays
	eps_x.resize(NumOfElectricMaterials_X);
	eps_y.resize(NumOfElectricMaterials_Y);
	eps_z.resize(NumOfElectricMaterials_Z);
	mu_x.resize(NumOfMagneticMaterials_X);
	mu_y.resize(NumOfMagneticMaterials_Y);
	mu_z.resize(NumOfMagneticMaterials_Z);
	cond_e_x.resize(NumOfElectricMaterials_X);
	cond_e_y.resize(NumOfElectricMaterials_Y);
	cond_e_z.resize(NumOfElectricMaterials_Z);
	cond_h_x.resize(NumOfMagneticMaterials_X);
	cond_h_y.resize(NumOfMagneticMaterials_Y);
	cond_h_z.resize(NumOfMagneticMaterials_Z);

	//Initialize constitutive parameter vectors:
	//relative permittivities
	eps_x = 1;
	eps_y = 1;
	eps_z = 1;
	//relative permeabilities
	mu_x = 1;
	mu_y = 1;
	mu_z = 1;
	//electric conductivities
	cond_e_x = 0;
	cond_e_y = 0;
	cond_e_z = 0;
	//magnetic conductivities
	cond_h_x = 0;
	cond_h_y = 0;
	cond_h_z = 0;

	//electric permittivity for PEC material (does not matter unless permittivity profile is visualized)
	//This is supposed to represent "infinity". Better and more portable value for the "maximum double value" can be found later.
	eps_x(PEC) = LIBSTD_DBL_MAX;
	eps_y(PEC) = LIBSTD_DBL_MAX;
	eps_z(PEC) = LIBSTD_DBL_MAX;
	mu_x(PEC) = 1;
	mu_y(PEC) = 1;
	mu_z(PEC) = 1;
	//electric conductivity for PEC material (does not matter unless conductivity profile is visualized)
	//This is supposed to represent "infinity". Better and more portable value for the "maximum double value" can be found later.
	cond_e_x(PEC) = LIBSTD_DBL_MAX;
	cond_e_y(PEC) = LIBSTD_DBL_MAX;
	cond_e_z(PEC) = LIBSTD_DBL_MAX;
	cond_h_x(PEC) = 1;
	cond_h_y(PEC) = 1;
	cond_h_z(PEC) = 1;

	//Allocate and initialize PML kappa parameters
	kappa_e_x.resize(Range(1,NCELLS_X+2*NPML+1));
	kappa_e_y.resize(Range(1,NCELLS_Y+2*NPML+1));
	kappa_e_z.resize(Range(1,NCELLS_Z+2*NPML+1));
	kappa_h_x.resize(Range(1,NCELLS_X+2*NPML));
	kappa_h_y.resize(Range(1,NCELLS_Y+2*NPML));
	kappa_h_z.resize(Range(1,NCELLS_Z+2*NPML));
	kappa_e_x=1;
	kappa_h_x=1;
	kappa_e_y=1;
	kappa_h_y=1;
	kappa_e_z=1;
	kappa_h_z=1;

	//Initialize the layering arrays (layering defined in the z-direction)
	//x component of the E-field
	Layering_e_x.resize(Range(1,NCELLS_Z+2*NPML+1));
	//y component of the E-field
	Layering_e_y.resize(Range(1,NCELLS_Z+2*NPML+1));
	//z component of the E-field
	Layering_e_z.resize(Range(1,NCELLS_Z+2*NPML));
	//x component of the H-field
	Layering_h_x.resize(Range(1,NCELLS_Z+2*NPML));
	//y component of the E-field
	Layering_h_y.resize(Range(1,NCELLS_Z+2*NPML));
	//z component of the E-field
	Layering_h_z.resize(Range(1,NCELLS_Z+2*NPML+1));
	Layering_e_x = vacuum;
	Layering_e_y = vacuum;
	Layering_e_z = vacuum;
	Layering_h_x = vacuum;
	Layering_h_y = vacuum;
	Layering_h_z = vacuum;

	Ca_X.resize(NumOfElectricMaterials_X);
	Cb_X.resize(NumOfElectricMaterials_X);
	Ca_Y.resize(NumOfElectricMaterials_Y);
	Cb_Y.resize(NumOfElectricMaterials_Y);
	Ca_Z.resize(NumOfElectricMaterials_Z);
	Cb_Z.resize(NumOfElectricMaterials_Z);
	Da_X.resize(NumOfMagneticMaterials_X);
	Db_X.resize(NumOfMagneticMaterials_X);
	Da_Y.resize(NumOfMagneticMaterials_Y);
	Db_Y.resize(NumOfMagneticMaterials_Y);
	Da_Z.resize(NumOfMagneticMaterials_Z);
	Db_Z.resize(NumOfMagneticMaterials_Z);

	// Update coefficients for vacuum material
	double eaf_x = dt*cond_e_x(vacuum)/(2.0*eps_x(vacuum)*epsilon_0);
	Ca_X(vacuum)=(1-eaf_x)/(1+eaf_x);
	Cb_X(vacuum)=dt/eps_x(vacuum)/epsilon_0/dx/(1+eaf_x);
	double eaf_y = dt*cond_e_y(vacuum)/(2.0*eps_y(vacuum)*epsilon_0);
	Ca_Y(vacuum)=(1-eaf_y)/(1+eaf_y);
	Cb_Y(vacuum)=dt/eps_y(vacuum)/epsilon_0/dx/(1+eaf_y);
	double eaf_z = dt*cond_e_z(vacuum)/(2.0*eps_z(vacuum)*epsilon_0);
	Ca_Z(vacuum)=(1-eaf_z)/(1+eaf_z);
	Cb_Z(vacuum)=dt/eps_z(vacuum)/epsilon_0/dx/(1+eaf_z);
	double haf_x = dt*cond_h_x(vacuum)/(2.0*mu_x(vacuum)*mu_0);
	Da_X(vacuum)=(1-haf_x)/(1+haf_x);
	Db_X(vacuum)=dt/mu_x(vacuum)/mu_0/dx/(1+haf_x);
	double haf_y = dt*cond_h_y(vacuum)/(2.0*mu_y(vacuum)*mu_0);
	Da_Y(vacuum)=(1-haf_y)/(1+haf_y);
	Db_Y(vacuum)=dt/mu_y(vacuum)/mu_0/dx/(1+haf_y);
	double haf_z = dt*cond_h_z(vacuum)/(2.0*mu_z(vacuum)*mu_0);
	Da_Z(vacuum)=(1-haf_z)/(1+haf_z);
	Db_Z(vacuum)=dt/mu_z(vacuum)/mu_0/dx/(1+haf_z);
	// Update coefficients for PEC material
	Ca_X(PEC)=-1;
	Cb_X(PEC)=0;
	Ca_Y(PEC)=-1;
	Cb_Y(PEC)=0;
	Ca_Z(PEC)=-1;
	Cb_Z(PEC)=0;
	Da_X(PEC)=Da_X(vacuum);
	Db_X(PEC)=Db_X(vacuum);
	Da_Y(PEC)=Da_Y(vacuum);
	Db_Y(PEC)=Db_Y(vacuum);
	Da_Z(PEC)=Da_Z(vacuum);
	Db_Z(PEC)=Db_Z(vacuum);
} //init_geom

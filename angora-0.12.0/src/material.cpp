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

//Includes the definitions of the data structures associated with material indices

#include "headers.h"

#include "material.h"

extern double dx,dt;

extern int NumOfElectricMaterials_X;
extern int NumOfElectricMaterials_Y;
extern int NumOfElectricMaterials_Z;
extern int NumOfMagneticMaterials_X;
extern int NumOfMagneticMaterials_Y;
extern int NumOfMagneticMaterials_Z;

extern const int vacuum;

extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<double,1> eps_x,eps_y,eps_z;
extern Array<double,1> mu_x,mu_y,mu_z;
extern Array<double,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<double,1> cond_h_x,cond_h_y,cond_h_z;


void AddIsotropicElectricMaterial(MaterialId& NewMaterialId, const double& epsilon_r, const double& sigma_e)
{//convenience function for adding an isotropic electric material
	AddDiagonalElectricMaterial(NewMaterialId,epsilon_r,epsilon_r,epsilon_r,sigma_e,sigma_e,sigma_e);
}

void AddDiagonalElectricMaterial(MaterialId& NewMaterialId,
								const double& epsilon_r_x, const double& epsilon_r_y, const double& epsilon_r_z,
								const double& sigma_e_x, const double& sigma_e_y, const double& sigma_e_z)
//Adds a new diagonally-anisotropic electric material into the material list, returns the material index of the created material
{
	//increment the number of materials by 1 in x,y,z directions
	NumOfElectricMaterials_X+=1;
	NumOfElectricMaterials_Y+=1;
	NumOfElectricMaterials_Z+=1;
	//Resize the constitutive parameter arrays to include the new material type
	eps_x.resizeAndPreserve(NumOfElectricMaterials_X);
	eps_y.resizeAndPreserve(NumOfElectricMaterials_Y);
	eps_z.resizeAndPreserve(NumOfElectricMaterials_Z);
	cond_e_x.resizeAndPreserve(NumOfElectricMaterials_X);
	cond_e_y.resizeAndPreserve(NumOfElectricMaterials_Y);
	cond_e_z.resizeAndPreserve(NumOfElectricMaterials_Z);
	//Resize the update coefficient arrays to include the new material type
	Ca_X.resizeAndPreserve(NumOfElectricMaterials_X);
	Cb_X.resizeAndPreserve(NumOfElectricMaterials_X);
	Ca_Y.resizeAndPreserve(NumOfElectricMaterials_Y);
	Cb_Y.resizeAndPreserve(NumOfElectricMaterials_Y);
	Ca_Z.resizeAndPreserve(NumOfElectricMaterials_Z);
	Cb_Z.resizeAndPreserve(NumOfElectricMaterials_Z);

	//indices of the newest material
	int newmaterial_e_x = NumOfElectricMaterials_X-1;
	int newmaterial_e_y = NumOfElectricMaterials_Y-1;
	int newmaterial_e_z = NumOfElectricMaterials_Z-1;

	//New constitutive parameters for the new material
	eps_x(newmaterial_e_x) = epsilon_r_x;
	eps_y(newmaterial_e_y) = epsilon_r_y;
	eps_z(newmaterial_e_z) = epsilon_r_z;
	cond_e_x(newmaterial_e_x) = sigma_e_x;
	cond_e_y(newmaterial_e_y) = sigma_e_y;
	cond_e_z(newmaterial_e_z) = sigma_e_z;
	//New update coefficients for the new material
	double eaf_x = dt*cond_e_x(newmaterial_e_x)/(2.0*eps_x(newmaterial_e_x)*epsilon_0);
	Ca_X(newmaterial_e_x)=(1-eaf_x)/(1+eaf_x);
	Cb_X(newmaterial_e_x)=dt/eps_x(newmaterial_e_x)/epsilon_0/dx/(1+eaf_x);
	double eaf_y = dt*cond_e_y(newmaterial_e_y)/(2.0*eps_y(newmaterial_e_y)*epsilon_0);
	Ca_Y(newmaterial_e_y)=(1-eaf_y)/(1+eaf_y);
	Cb_Y(newmaterial_e_y)=dt/eps_y(newmaterial_e_y)/epsilon_0/dx/(1+eaf_y);
	double eaf_z = dt*cond_e_z(newmaterial_e_z)/(2.0*eps_z(newmaterial_e_z)*epsilon_0);
	Ca_Z(newmaterial_e_z)=(1-eaf_z)/(1+eaf_z);
	Cb_Z(newmaterial_e_z)=dt/eps_z(newmaterial_e_z)/epsilon_0/dx/(1+eaf_z);

	//write the indices of the new electric material into the material identifier object
	NewMaterialId.ElectricIndex_X = newmaterial_e_x;
	NewMaterialId.ElectricIndex_Y = newmaterial_e_y;
	NewMaterialId.ElectricIndex_Z = newmaterial_e_z;
	//magnetic properties are assumed the same as vacuum
	NewMaterialId.MagneticIndex_X = vacuum;
	NewMaterialId.MagneticIndex_Y = vacuum;
	NewMaterialId.MagneticIndex_Z = vacuum;
}

void AddIsotropicMagneticMaterial(MaterialId& NewMaterialId, const double& mu_r, const double& sigma_h)
{//convenience function for adding an isotropic magnetic material
	AddDiagonalMagneticMaterial(NewMaterialId,mu_r,mu_r,mu_r,sigma_h,sigma_h,sigma_h);
}

void AddDiagonalMagneticMaterial(MaterialId& NewMaterialId,
								const double& mu_r_x, const double& mu_r_y, const double& mu_r_z,
								const double& sigma_h_x, const double& sigma_h_y, const double& sigma_h_z)
//Adds a new diagonally-anisotropic magnetic material into the material list, returns the material index of the created material
{
	//increment the number of materials by 1 in x,y,z directions
	NumOfMagneticMaterials_X+=1;
	NumOfMagneticMaterials_Y+=1;
	NumOfMagneticMaterials_Z+=1;
	//Resize the constitutive parameter arrays to include the new material type
	mu_x.resizeAndPreserve(NumOfMagneticMaterials_X);
	mu_y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	mu_z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	cond_h_x.resizeAndPreserve(NumOfMagneticMaterials_X);
	cond_h_y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	cond_h_z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	//Resize the update coefficient arrays to include the new material type
	Da_X.resizeAndPreserve(NumOfMagneticMaterials_X);
	Db_X.resizeAndPreserve(NumOfMagneticMaterials_X);
	Da_Y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	Db_Y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	Da_Z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	Db_Z.resizeAndPreserve(NumOfMagneticMaterials_Z);

	//indices of the newest material
	int newmaterial_h_x = NumOfMagneticMaterials_X-1;
	int newmaterial_h_y = NumOfMagneticMaterials_Y-1;
	int newmaterial_h_z = NumOfMagneticMaterials_Z-1;

	//New constitutive parameters for the new material
	mu_x(newmaterial_h_x) = mu_r_x;
	mu_y(newmaterial_h_y) = mu_r_y;
	mu_z(newmaterial_h_z) = mu_r_z;
	cond_h_x(newmaterial_h_x) = sigma_h_x;
	cond_h_y(newmaterial_h_y) = sigma_h_y;
	cond_h_z(newmaterial_h_z) = sigma_h_z;
	//New update coefficients for the new material
	double haf_x = dt*cond_h_x(newmaterial_h_x)/(2.0*mu_x(newmaterial_h_x)*mu_0);
	Da_X(newmaterial_h_x)=(1-haf_x)/(1+haf_x);
	Db_X(newmaterial_h_x)=dt/mu_x(newmaterial_h_x)/mu_0/dx/(1+haf_x);
	double haf_y = dt*cond_h_y(newmaterial_h_y)/(2.0*mu_y(newmaterial_h_y)*mu_0);
	Da_Y(newmaterial_h_y)=(1-haf_y)/(1+haf_y);
	Db_Y(newmaterial_h_y)=dt/mu_y(newmaterial_h_y)/mu_0/dx/(1+haf_y);
	double haf_z = dt*cond_h_z(newmaterial_h_z)/(2.0*mu_z(newmaterial_h_z)*mu_0);
	Da_Z(newmaterial_h_z)=(1-haf_z)/(1+haf_z);
	Db_Z(newmaterial_h_z)=dt/mu_z(newmaterial_h_z)/mu_0/dx/(1+haf_z);
	//write the indices of the new Magnetic material into the material identifier object
	NewMaterialId.MagneticIndex_X = newmaterial_h_x;
	NewMaterialId.MagneticIndex_Y = newmaterial_h_y;
	NewMaterialId.MagneticIndex_Z = newmaterial_h_z;
	//electric properties are assumed the same as vacuum
	NewMaterialId.ElectricIndex_X = vacuum;
	NewMaterialId.ElectricIndex_Y = vacuum;
	NewMaterialId.ElectricIndex_Z = vacuum;
}

void AddIsotropicMaterial(MaterialId& NewMaterialId, const double& epsilon_r, const double& mu_r, const double& sigma_e, const double& sigma_h)
//convenience function for creating an isotropic material
{
	AddDiagonalMaterial(NewMaterialId,epsilon_r,epsilon_r,epsilon_r,mu_r,mu_r,mu_r,sigma_e,sigma_e,sigma_e,sigma_h,sigma_h,sigma_h);
}

void AddDiagonalMaterial(MaterialId& NewMaterialId,
							const double& epsilon_r_x, const double& epsilon_r_y, const double& epsilon_r_z,
							const double& mu_r_x, const double& mu_r_y, const double& mu_r_z,
							const double& sigma_e_x, const double& sigma_e_y, const double& sigma_e_z,
							const double& sigma_h_x, const double& sigma_h_y, const double& sigma_h_z)
//Adds a new diagonally-anisotropic material into the material list, returns the material index of the created material
{
	//increment the number of materials by 1 in x,y,z directions
	NumOfElectricMaterials_X+=1;
	NumOfElectricMaterials_Y+=1;
	NumOfElectricMaterials_Z+=1;
	NumOfMagneticMaterials_X+=1;
	NumOfMagneticMaterials_Y+=1;
	NumOfMagneticMaterials_Z+=1;
	//Resize the constitutive parameter arrays to include the new material type
	eps_x.resizeAndPreserve(NumOfElectricMaterials_X);
	eps_y.resizeAndPreserve(NumOfElectricMaterials_Y);
	eps_z.resizeAndPreserve(NumOfElectricMaterials_Z);
	mu_x.resizeAndPreserve(NumOfMagneticMaterials_X);
	mu_y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	mu_z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	cond_e_x.resizeAndPreserve(NumOfElectricMaterials_X);
	cond_e_y.resizeAndPreserve(NumOfElectricMaterials_Y);
	cond_e_z.resizeAndPreserve(NumOfElectricMaterials_Z);
	cond_h_x.resizeAndPreserve(NumOfMagneticMaterials_X);
	cond_h_y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	cond_h_z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	//Resize the update coefficient arrays to include the new material type
	Ca_X.resizeAndPreserve(NumOfElectricMaterials_X);
	Cb_X.resizeAndPreserve(NumOfElectricMaterials_X);
	Ca_Y.resizeAndPreserve(NumOfElectricMaterials_Y);
	Cb_Y.resizeAndPreserve(NumOfElectricMaterials_Y);
	Ca_Z.resizeAndPreserve(NumOfElectricMaterials_Z);
	Cb_Z.resizeAndPreserve(NumOfElectricMaterials_Z);
	Da_X.resizeAndPreserve(NumOfMagneticMaterials_X);
	Db_X.resizeAndPreserve(NumOfMagneticMaterials_X);
	Da_Y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	Db_Y.resizeAndPreserve(NumOfMagneticMaterials_Y);
	Da_Z.resizeAndPreserve(NumOfMagneticMaterials_Z);
	Db_Z.resizeAndPreserve(NumOfMagneticMaterials_Z);

	//indices of the newest material
	int newmaterial_e_x = NumOfElectricMaterials_X-1;
	int newmaterial_e_y = NumOfElectricMaterials_Y-1;
	int newmaterial_e_z = NumOfElectricMaterials_Z-1;
	int newmaterial_h_x = NumOfMagneticMaterials_X-1;
	int newmaterial_h_y = NumOfMagneticMaterials_Y-1;
	int newmaterial_h_z = NumOfMagneticMaterials_Z-1;

	//New constitutive parameters for the new material
	eps_x(newmaterial_e_x) = epsilon_r_x;
	eps_y(newmaterial_e_y) = epsilon_r_y;
	eps_z(newmaterial_e_z) = epsilon_r_z;
	mu_x(newmaterial_h_x) = mu_r_x;
	mu_y(newmaterial_h_y) = mu_r_y;
	mu_z(newmaterial_h_z) = mu_r_z;
	cond_e_x(newmaterial_e_x) = sigma_e_x;
	cond_e_y(newmaterial_e_y) = sigma_e_y;
	cond_e_z(newmaterial_e_z) = sigma_e_z;
	cond_h_x(newmaterial_h_x) = sigma_h_x;
	cond_h_y(newmaterial_h_y) = sigma_h_y;
	cond_h_z(newmaterial_h_z) = sigma_h_z;
	//New update coefficients for the new material
	double eaf_x = dt*cond_e_x(newmaterial_e_x)/(2.0*eps_x(newmaterial_e_x)*epsilon_0);
	Ca_X(newmaterial_e_x)=(1-eaf_x)/(1+eaf_x);
	Cb_X(newmaterial_e_x)=dt/eps_x(newmaterial_e_x)/epsilon_0/dx/(1+eaf_x);
	double eaf_y = dt*cond_e_y(newmaterial_e_y)/(2.0*eps_y(newmaterial_e_y)*epsilon_0);
	Ca_Y(newmaterial_e_y)=(1-eaf_y)/(1+eaf_y);
	Cb_Y(newmaterial_e_y)=dt/eps_y(newmaterial_e_y)/epsilon_0/dx/(1+eaf_y);
	double eaf_z = dt*cond_e_z(newmaterial_e_z)/(2.0*eps_z(newmaterial_e_z)*epsilon_0);
	Ca_Z(newmaterial_e_z)=(1-eaf_z)/(1+eaf_z);
	Cb_Z(newmaterial_e_z)=dt/eps_z(newmaterial_e_z)/epsilon_0/dx/(1+eaf_z);
	//New update coefficients for the new material
	double haf_x = dt*cond_h_x(newmaterial_h_x)/(2.0*mu_x(newmaterial_h_x)*mu_0);
	Da_X(newmaterial_h_x)=(1-haf_x)/(1+haf_x);
	Db_X(newmaterial_h_x)=dt/mu_x(newmaterial_h_x)/mu_0/dx/(1+haf_x);
	double haf_y = dt*cond_h_y(newmaterial_h_y)/(2.0*mu_y(newmaterial_h_y)*mu_0);
	Da_Y(newmaterial_h_y)=(1-haf_y)/(1+haf_y);
	Db_Y(newmaterial_h_y)=dt/mu_y(newmaterial_h_y)/mu_0/dx/(1+haf_y);
	double haf_z = dt*cond_h_z(newmaterial_h_z)/(2.0*mu_z(newmaterial_h_z)*mu_0);
	Da_Z(newmaterial_h_z)=(1-haf_z)/(1+haf_z);
	Db_Z(newmaterial_h_z)=dt/mu_z(newmaterial_h_z)/mu_0/dx/(1+haf_z);

	//write the indices of the new electric material into the material identifier object
	NewMaterialId.ElectricIndex_X = newmaterial_e_x;
	NewMaterialId.ElectricIndex_Y = newmaterial_e_y;
	NewMaterialId.ElectricIndex_Z = newmaterial_e_z;
	NewMaterialId.MagneticIndex_X = newmaterial_h_x;
	NewMaterialId.MagneticIndex_Y = newmaterial_h_y;
	NewMaterialId.MagneticIndex_Z = newmaterial_h_z;
}

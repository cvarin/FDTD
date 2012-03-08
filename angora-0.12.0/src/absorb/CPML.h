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

#ifndef CPML_H
#define CPML_H

//Declaration of the PML object "CPML"


class Cpml
{
 public:
	 Cpml(const int& my_pml_thickness, const double& SizeOfScatterer);
	 void UpdateE();		//Do the E-field updates in the PML region
	 void UpdateH();		//Do the H-field updates in the PML region
 private:
	 //PML Current Arrays
	 //FrontX : (NPML)x(NCELLS_Y+NPML)x(NCELLS_Z+NPML)	cells
	 //RightY : (NCELLS_X+NPML)x(NPML)x(NCELLS_Z+NPML)	cells
	 //UpperZ : (NCELLS_X+NPML)x(NCELLS_Y+NPML)xNPML	cells
	 Array<double,3>
		 Psi_Eyx_BackX,Psi_Ezx_BackX,Psi_Eyx_FrontX,Psi_Ezx_FrontX,Psi_Exy_LeftY,Psi_Ezy_LeftY,Psi_Exy_RightY,Psi_Ezy_RightY,Psi_Exz_UpperZ,Psi_Eyz_UpperZ,Psi_Exz_LowerZ,Psi_Eyz_LowerZ,
		 Psi_Hyx_BackX,Psi_Hzx_BackX,Psi_Hyx_FrontX,Psi_Hzx_FrontX,Psi_Hxy_LeftY,Psi_Hzy_LeftY,Psi_Hxy_RightY,Psi_Hzy_RightY,Psi_Hxz_UpperZ,Psi_Hyz_UpperZ,Psi_Hxz_LowerZ,Psi_Hyz_LowerZ;
	 //PML Parameters
	 int m;	//Order for polynomial grading
	 double epsilon_eff_x,epsilon_eff_y,epsilon_eff_lowerz,epsilon_eff_upperz;
	 double sigma_coeff,sigma_x_max,sigma_y_max,sigma_lowerz_max,sigma_upperz_max;
	 double kappa_x_max,kappa_y_max,kappa_z_max;
	 double alpha_x_max,alpha_y_max,alpha_z_max;
	 double alpha_x_min,alpha_y_min,alpha_z_min;
	 //PML Update Arrays (1-D)
	 Array<double,1> Apml_e_backx,Apml_e_frontx,Apml_e_lefty,Apml_e_righty,Apml_e_lowerz,Apml_e_upperz,
		 Bpml_e_backx,Bpml_e_frontx,Bpml_e_lefty,Bpml_e_righty,Bpml_e_lowerz,Bpml_e_upperz;
	 Array<double,1> Apml_h_backx,Apml_h_frontx,Apml_h_lefty,Apml_h_righty,Apml_h_lowerz,Apml_h_upperz,
		 Bpml_h_backx,Bpml_h_frontx,Bpml_h_lefty,Bpml_h_righty,Bpml_h_lowerz,Bpml_h_upperz;

	 const int pml_thickness;

	 void UpdateE_X();
	 void UpdateH_X();
	 void UpdateE_Y();
	 void UpdateH_Y();
	 void UpdateE_Z();
	 void UpdateH_Z();

};
#endif

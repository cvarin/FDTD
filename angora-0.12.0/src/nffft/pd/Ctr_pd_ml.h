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

#ifndef CTR_PD_ML_H
#define CTR_PD_ML_H

//Declaration of the class "Ctr_pd_ml" for a PHASOR-DOMAIN near-field-to-far-field transformer in a general N-layered medium

#include "headers.h"

#include "Ctr_pd.h"

#define UPPERMOST_LAYER_INDEX number_of_layers-1
#define LOWERMOST_LAYER_INDEX 0

#define UPPERMOST_INTERFACE_INDEX max(0,number_of_layers-2)
#define LOWERMOST_INTERFACE_INDEX 0


class Ctr_pd_ml : public Ctr_pd
{
 public:
	 Ctr_pd_ml(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index);		//constructor
	 virtual ~Ctr_pd_ml() {};  //virtual destructor for virtual base class

	 virtual void ConstructFarField();	//(virtual from Ctr_pd) Construct far-field from the stored equivalent-current phasors

 protected:
 	 virtual void ConstructPotential_A(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the magnetic potential A
 	 virtual void ConstructPotential_F(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the electric potential F

 	 virtual void Calculate_TL_GreenFunctions(const int& l, const int& d1, const int& d2);	//calculates the TL Green functions for the given l,d1,d2 triplet

 	 int TemporaryFarFieldOriginZ; //z-index of the cell used temporarily as the origin of the NFFFT, until the post-processing step

	 double threshold;	//what is the threshold for sin(theta) that determines upper or lower half space?

	 //the following are for a given triplet (l,d1,d2) representing an observation direction
	 double epsilon_r_o;	//relative permittivity in the observation half space
	 double mu_r_o;	//relative permeability in the observation half space
	 //wavenumber at the observation half space, corrected for grid-dispersion
	 //(first dim: frequency, second&third dims: first&second direction parameters)
	 double kk_o_g;		//in m^-1

     //temporary variables for fast post-processing
     //sines and cosines
     double cosP,sinP,cosT,sinT;
     double sinTcosP,sinTsinP;
     //wavenumber in the observation half space, multiplied by the grid spacing
     double ko_dx;

	 //Theoretical far-field calculations
	 virtual complex<double> TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical theta-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
	 virtual complex<double> TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical phi-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)

 private:
	int number_of_interfaces;

	Array<double,1> epsilon_r;
	Array<double,1> mu_r;
	Array<double,1> cond_e;
	Array<double,1> grounded;

	Array<complex<double>,1> eps_r,m_r;

	complex<double> obliquityfactor;
	Array<complex<double>,1> kt_g;

	Array<complex<double>,1> eta_e,eta_h; //transmission-line impedances (TM and TE)

	Array<complex<double>,1> Z_e,Z_h; //(total) impedances at each planar interface (for TM and TE waves, respectively)

	Array<complex<double>,1> V_v_e_minus,V_v_e_plus,V_v_h_minus,V_v_h_plus,V_i_e_minus,V_i_e_plus,V_i_h_minus,V_i_h_plus;

	Array<complex<double>,1> kt_dx;

	Array<complex<double>,1> Green_A_theta_J_x_minus,Green_A_theta_J_x_plus,Green_A_theta_J_y_minus,Green_A_theta_J_y_plus,Green_A_theta_J_z_minus,Green_A_theta_J_z_plus,Green_A_phi_J_x_minus,Green_A_phi_J_x_plus,Green_A_phi_J_y_minus,Green_A_phi_J_y_plus,Green_F_theta_M_x_minus,Green_F_theta_M_x_plus,Green_F_theta_M_y_minus,Green_F_theta_M_y_plus,Green_F_theta_M_z_minus,Green_F_theta_M_z_plus,Green_F_phi_M_x_minus,Green_F_phi_M_x_plus,Green_F_phi_M_y_minus,Green_F_phi_M_y_plus;

	double z_pr_minus_z_n;	//z_pr minus the coordinate of the next interface below (always positive)
	double z_pr_minus_z_n_1;	//z_pr minus the coordinate of the next interface above (always negative)

	Array<double,1> interface_coord; //z coordinates of each layer interface with respect the temporary NFFFT origin
	Array<double,1> d_n; //thicknesses of each layer (all nonnegative), defined for n=0...N-1, arbitrary values at n=0 and n=N-1 (they do not affect the recursion)

	double recursion_max_error;	//the maximum allowable error that determines the stopping condition in the recursion algorithms

	// Hertzian-dipole variables:
	complex<double> theta_component,phi_component,total_theta_component,total_phi_component;	//field components are complex for freq.domain NFFFT
	string Orientation;
	int x,y,z;
	double SourceX,SourceY; //x and y coordinates of the Hertzian source (with respect to the grid origin)
	double SourceZ_minus_z_n; //z-coordinate of the Hertzian source minus the coordinate of the next interface below (always positive)
	double SourceZ_minus_z_n_1; //z-coordinate of the Hertzian source minus the coordinate of the next interface above (always negative)
};
#endif // CTR_PD_ML_H

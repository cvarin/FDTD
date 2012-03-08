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

#ifndef CTR_PD_2L_H
#define CTR_PD_2L_H

//Declaration of the class "Ctr_pd_2l" for a PHASOR-DOMAIN near-field-to-far-field transformer in a 2-layered medium
//Derived from the abstract base class "Ctr_pd".

#include "Ctr_pd_ml.h"


class Ctr_pd_2l : public Ctr_pd_ml
{
 public:
	 Ctr_pd_2l(const TrDataType_pd& MyData, const string& FarFieldFileName, const int& Index);		//constructor

 private:
 	 virtual void ConstructPotential_A(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the magnetic potential A
 	 virtual void ConstructPotential_F(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) computes the theta and phi components of the electric potential F

 	 virtual void Calculate_TL_GreenFunctions(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd_ml) calculates the TL Green functions for the given l,d1,d2 triplet

     //relative permittivities of the upper and lower half spaces
	 double epsilon_r_0,epsilon_r_1;
	 //relative permeabilities of the upper and lower half spaces
	 double mu_r_0,mu_r_1;
	 //electrical conductivities of the upper and lower half spaces
	 double cond_e_0,cond_e_1;

	 int LowerHalfSpaceUpper;		//index of the uppermost cell in the lower half space layer in the 2-layer medium

	 //epsilon may be complex for frequency-selective materials (lossy, etc)
	 complex<double> eps_r_0,eps_r_1;	//relative permittivities of the upper and lower half spaces, normalized by epsilon_r_o
	 double m_r_0,m_r_1;	//relative permeabilities of the upper and lower half spaces, normalized by mu_r_o

	 complex<double> obliquityfactor_0,obliquityfactor_1; //obliquity factors that are multiplied by the (dispersion-corrected) wave propagation constant kk_o_g to yield the TL propagation constants kt_g_0,kt_g_1.
 	 complex<double> kt_g_0,kt_g_1;		//transmission-line propagation constants in the upper and lower half spaces (grid-dispersion-corrected)
 	 //NOTE: The type is "complex", because the propagation constant can become complex in the following cases:
	 // (a) The layers are lossy
	 // (b) (for some observation angles in the observation space) the layer in which the source resides has lower permittivity than the observation space -- e.g. For a source in the lower half space, the propagation constant can be complex in the lower half space in total-internal-reflection (TIR) range, if the upper half space has higher permittivity than the lower)
	 // The latter [(b)] only happens when the observation angle is in the "opposite" of the source space. For example, if the source is in the upper half space, the TL impedances are complex in the TIR region in the lower half space. Creation of propagating plane waves in this region can occur, provided that the evanescent waves in the upper half space can reach the interface. Because the evanescent waves are bounded in space, this coupling only happens for sources very close to the interface.

	 double z_pr;	//vertical distance (in cells) from the NFFFT origin

	 //in the following, TE and TM refer to "with respect to the axis of invariance (x axis)"
	 //this is the opposite of the definition adopted in the 3D transformer code
	 complex<double> Z_e_0,Z_e_1;	//TE transmission-line impedances in the upper and lower half spaces
	 complex<double> Z_h_0,Z_h_1;	//TM transmission-line impedances in the upper and lower half spaces

	 complex<double> Refl_e_10,Trans_e_10;	//TE refl./trans. coefficients passing from the lower to upper half space
	 complex<double> Refl_e_01,Trans_e_01;	//TE refl./trans. coefficients passing from the upper to lower half space

	 complex<double> Refl_h_10,Trans_h_10;	//TM refl./trans. coefficients passing from the lower to upper half space
	 complex<double> Refl_h_01,Trans_h_01;	//TM refl./trans. coefficients passing from the upper to lower half space

	 //transmission-line Green functions
	 complex<double> V_v_e_upward_0,V_v_e_downward_0,V_v_h_upward_0,V_v_h_downward_0;
	 complex<double> V_i_e_upward_0,V_i_e_downward_0,V_i_h_upward_0,V_i_h_downward_0;
	 complex<double> V_v_e_upward_1,V_v_e_downward_1,V_v_h_upward_1,V_v_h_downward_1;
	 complex<double> V_i_e_upward_1,V_i_e_downward_1,V_i_h_upward_1,V_i_h_downward_1;

     //transmission-line wavenumbers in different layers, multiplied by the grid spacing
     complex<double> kt0_dx,kt1_dx;

	 //contributions to the potentials A_theta,A_phi,F_theta,F_phi by the equivalent currents J,M
	 //these are temporary variables used in the post-processing loop
	 //upward contributions have exp(+jk^{p}z') dependency, while downward contributions have exp(-jk^{p}z') dependency
	 //upper half space
	 complex<double> Green_A_theta_J_x_upward_0,Green_A_theta_J_x_downward_0,
			 		Green_A_theta_J_y_upward_0,Green_A_theta_J_y_downward_0,
  					Green_A_theta_J_z_upward_0,Green_A_theta_J_z_downward_0;
	 complex<double> Green_A_phi_J_x_upward_0,Green_A_phi_J_x_downward_0,
  					Green_A_phi_J_y_upward_0,Green_A_phi_J_y_downward_0;
	 complex<double> Green_F_theta_M_x_upward_0,Green_F_theta_M_x_downward_0,
			 		Green_F_theta_M_y_upward_0,Green_F_theta_M_y_downward_0,
  					Green_F_theta_M_z_upward_0,Green_F_theta_M_z_downward_0;
	 complex<double> Green_F_phi_M_x_upward_0,Green_F_phi_M_x_downward_0,
  					Green_F_phi_M_y_upward_0,Green_F_phi_M_y_downward_0;
	 //lower half space
	 complex<double> Green_A_theta_J_x_upward_1,Green_A_theta_J_x_downward_1,
			 		Green_A_theta_J_y_upward_1,Green_A_theta_J_y_downward_1,
  					Green_A_theta_J_z_upward_1,Green_A_theta_J_z_downward_1;
	 complex<double> Green_A_phi_J_x_upward_1,Green_A_phi_J_x_downward_1,
  					Green_A_phi_J_y_upward_1,Green_A_phi_J_y_downward_1;
	 complex<double> Green_F_theta_M_x_upward_1,Green_F_theta_M_x_downward_1,
			 		Green_F_theta_M_y_upward_1,Green_F_theta_M_y_downward_1,
  					Green_F_theta_M_z_upward_1,Green_F_theta_M_z_downward_1;
	 complex<double> Green_F_phi_M_x_upward_1,Green_F_phi_M_x_downward_1,
  					Green_F_phi_M_y_upward_1,Green_F_phi_M_y_downward_1;

	 //Theoretical far-field calculations
	 virtual complex<double> TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical theta-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
	 virtual complex<double> TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2);	//(virtual from Ctr_pd) the theoretical phi-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)

	// Hertzian-dipole variables:
	 complex<double> theta_component,phi_component,total_theta_component,total_phi_component;	//field components are complex for freq.domain NFFFT
//	 double temp;
	 string Orientation;
	 int x,y,z;
	 double SourceX,SourceY; //x and y coordinates of the Hertzian source (with respect to the grid origin)
	 double SourceZ; //z-component of the Hertzian source (with respect to the uppermost interface)
};
#endif

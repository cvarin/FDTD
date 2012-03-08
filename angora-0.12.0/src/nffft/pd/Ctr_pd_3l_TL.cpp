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

//Defines the transmission-line parameters in the PHASOR-DOMAIN near-field-to-far-field transformer object "Ctr_pd_3l" in a 3-layered medium

//(see Ctr_pd_ml_TL.cpp for details on the Green's function definitions)

#include "headers.h"

#include "Ctr_pd_3l.h"

extern double dx;

extern int number_of_layers;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;


void Ctr_pd_3l::Calculate_TL_GreenFunctions(const int& l, const int& d1, const int& d2)	//calculates the TL Green functions (functional form, does not require storage)
{
    //temporary variables for fast post-processing
    //sines and cosines
    cosP = CosP(l,d1,d2);
    sinP = SinP(l,d1,d2);
    cosT = CosT(l,d1,d2);
    sinT = SinT(l,d1,d2);
    sinTcosP = sinT*cosP;
    sinTsinP = sinT*sinP;

	if (cosT>abs(threshold))		//upper half space
	{
		epsilon_r_o = epsilon_r_upper;	//rel. permittivity of the upper half space
		mu_r_o = mu_r_upper;	//rel. permeability of the upper half space
	}
	else if (cosT<-abs(threshold))		//lower half space
	{
		epsilon_r_o = epsilon_r_lower;	//rel. permittivity of the lower half space
		mu_r_o = mu_r_lower;	//rel. permeability of the lower half space
	}
	//wavenumber at the observation half space, corrected for grid-dispersion  (in m^-1)
	/** FIXME: I don't know why this works (or seems to work) for lossy observation space **/
	/** It seems as though kk_o_g (as well as epsilon_r_o) should be complex in that case **/
	kk_o_g = (ww(l)/c/GridVelocity(ww(l),THETA(l,d1,d2),PHI(l,d1,d2)))*sqrt(epsilon_r_o*mu_r_o);

	//(complex) rel. permittivities of the layers,  normalized by epsilon_r_o
	eps_r_0 = (epsilon_r_0-ii*cond_e_0/epsilon_0/ww(l))/epsilon_r_o;
	eps_r_1 = (epsilon_r_1-ii*cond_e_1/epsilon_0/ww(l))/epsilon_r_o;
	eps_r_2 = (epsilon_r_2-ii*cond_e_2/epsilon_0/ww(l))/epsilon_r_o;
	m_r_0 = mu_r_0/mu_r_o;		//rel. permeability of the upper half space, normalized by mu_r_o
	m_r_1 = mu_r_1/mu_r_o;		//rel. permeability of the middle layer, normalized by mu_r_o
	m_r_2 = mu_r_2/mu_r_o;		//rel. permeability of the lower half space, normalized by mu_r_o

	obliquityfactor_0 = sqrt(complex<double>(eps_r_0*m_r_0-pow2(sinT)));
	obliquityfactor_1 = sqrt(complex<double>(eps_r_1*m_r_1-pow2(sinT)));
	obliquityfactor_2 = sqrt(complex<double>(eps_r_2*m_r_2-pow2(sinT)));
	// kt_g_0,kt_g_1,kt_g_2 should have negative imaginary parts; therefore so should these obliquity factors.
	// Note that imaginary kt_g_% occurs when the source is radiating to the "opposite" half space. In this case, for example, the "downward" Green function contributes to the far field if the source is in the upper half space, and the obs. half space is the lower half space. In the reverse case, the "upward" Green function contributes. In both cases, it turns out that the longitudinal (z-dependent) phase term is in the form of exp(-ii*kt_g_%*dx*|z_pr|), where |z_pr| is the absolute value of the longitudinal distance from the interface. Therefore, kt_g_% should have negative imaginary part in order to make the term decay with |z_pr|.
	if (imag(obliquityfactor_0)>0)
	{
		obliquityfactor_0 *= -1.0;
	}
	if (imag(obliquityfactor_1)>0)
	{
		obliquityfactor_1 *= -1.0;
	}
	if (imag(obliquityfactor_2)>0)
	{
		obliquityfactor_2 *= -1.0;
	}

	//transmission-line propagation constants
	kt_g_0 = obliquityfactor_0*kk_o_g;
	kt_g_1 = obliquityfactor_1*kk_o_g;
	kt_g_2 = obliquityfactor_2*kk_o_g;

	//Transmission-line impedances (normalized by Z, since only V_i_{e,h}/Z is used in the formulas)
	//These are defined in an indirect manner in terms of kt_g_0,kt_g_1,kt_g_2; the advantage being that the sign of the imaginary part of Z_{e,h}_{0,1,2} is automatically determined by that of kt_g_0,kt_g_1,kt_g_2. The choice of sign of the latter is explained above.
	//TM impedances
	Z_e_0 = obliquityfactor_0/eps_r_0;
	Z_e_1 = obliquityfactor_1/eps_r_1;
	Z_e_2 = obliquityfactor_2/eps_r_2;
	//TE impedances
	Z_h_0 = pow(obliquityfactor_0,-1)*m_r_0;
	Z_h_1 = pow(obliquityfactor_1,-1)*m_r_1;
	Z_h_2 = pow(obliquityfactor_2,-1)*m_r_2;

	//Transmission and reflection coefficients
	//TM transmission and reflection coefficients
	Refl_e_10 = ((Z_e_0-Z_e_1)/(Z_e_0+Z_e_1));
	Trans_e_10 = 1.0 + Refl_e_10;
	Refl_e_01 = -Refl_e_10;
	Trans_e_01 = 1.0 + Refl_e_01;

	if (grounded)
	{
		Refl_e_21 = -1;
		Trans_e_21 = 0;
		Refl_e_12 = -1;
		Trans_e_12 = 0;
	}
	else
	{
		Refl_e_21 = ((Z_e_1-Z_e_2)/(Z_e_1+Z_e_2));
		Trans_e_21 = 1.0 + Refl_e_21;
		Refl_e_12 = -Refl_e_21;
		Trans_e_12 = 1.0 + Refl_e_12;
	}

	//TE transmission and reflection coefficients
	Refl_h_10 = ((Z_h_0-Z_h_1)/(Z_h_0+Z_h_1));
	Trans_h_10 = 1.0 + Refl_h_10;
	Refl_h_01 = -Refl_h_10;
	Trans_h_01 = 1.0 + Refl_h_01;

	if (grounded)
	{
		Refl_h_21 = -1;
		Trans_h_21 = 0;
		Refl_h_12 = -1;
		Trans_h_12 = 0;
	}
	else
	{
		Refl_h_21 = ((Z_h_1-Z_h_2)/(Z_h_1+Z_h_2));
		Trans_h_21 = 1.0 + Refl_h_21;
		Refl_h_12 = -Refl_h_21;
		Trans_h_12 = 1.0 + Refl_h_12;
	}

	//Fill the auxiliary arrays that hold the contributions of the infinite reverberations inside the middle layer
	//if the observation angle is in the upper half space:
	V_e_upper = Trans_e_10/(1.0-Refl_e_10*Refl_e_12*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));
	V_h_upper = Trans_h_10/(1.0-Refl_h_10*Refl_h_12*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));
	//if the observation angle is in the lower half space:
	V_e_lower = Trans_e_12/(1.0-Refl_e_10*Refl_e_12*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));
	V_h_lower = Trans_h_12/(1.0-Refl_h_10*Refl_h_12*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));
	//(only difference is in the prefactor (Trans_e_10 or Trans_e_12))

	//TL Green functions V_v_e,V_v_h
	if (cosT>abs(threshold))		//observation angle is in the upper half space AND the imaginary permittivity in that half space is negligible
	{
		V_v_e_upward_0 = 1;
		V_v_e_downward_0 = -(Refl_e_01+Trans_e_01*Refl_e_12*V_e_upper*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));
		V_v_h_upward_0 = 1;
		V_v_h_downward_0 = -(Refl_h_01+Trans_h_01*Refl_h_12*V_h_upper*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness)));

		V_v_e_upward_1 = V_e_upper;
		V_v_e_downward_1 = -Refl_e_12*V_e_upper*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness));
		V_v_h_upward_1 = V_h_upper;
		V_v_h_downward_1 = -Refl_h_12*V_h_upper*exp(-ii*kt_g_1*(2*dx*MiddleLayerThickness));

		V_v_e_upward_2 = Trans_e_21*V_e_upper*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_e_downward_2 = 0;
		V_v_h_upward_2 = Trans_h_21*V_h_upper*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_h_downward_2 = 0;
	}
	else if (cosT<-abs(threshold))		//observation angle is in the lower half space
	{
		V_v_e_upward_0 = 0;
		V_v_e_downward_0 = Trans_e_01*V_e_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_h_upward_0 = 0;
		V_v_h_downward_0 = Trans_h_01*V_h_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));

		V_v_e_upward_1 = -Refl_e_10*V_e_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_e_downward_1 = V_e_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_h_upward_1 = -Refl_h_10*V_h_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));
		V_v_h_downward_1 = V_h_lower*exp(ii*(kt_g_2-kt_g_1)*(dx*MiddleLayerThickness));

		V_v_e_upward_2 = -(Refl_e_21*exp(ii*kt_g_2*(2*dx*MiddleLayerThickness))
				+Trans_e_21*Refl_e_10*V_e_lower*exp(ii*(kt_g_2-kt_g_1)*(2*dx*MiddleLayerThickness)));
		V_v_e_downward_2 = 1;
		V_v_h_upward_2 = -(Refl_h_21*exp(ii*kt_g_2*(2*dx*MiddleLayerThickness))
				+Trans_h_21*Refl_h_10*V_h_lower*exp(ii*(kt_g_2-kt_g_1)*(2*dx*MiddleLayerThickness)));
		V_v_h_downward_2 = 1;
	}
	else
	{
		V_v_e_upward_0 = 0;
		V_v_e_downward_0 = 0;
		V_v_h_upward_0 = 0;
		V_v_h_downward_0 = 0;

		V_v_e_upward_1 = 0;
		V_v_e_downward_1 = 0;
		V_v_h_upward_1 = 0;
		V_v_h_downward_1 = 0;

		V_v_e_upward_2 = 0;
		V_v_e_downward_2 = 0;
		V_v_h_upward_2 = 0;
		V_v_h_downward_2 = 0;
	}

	//TL Green functions V_i_e,V_i_h
	V_i_e_upward_0 = V_v_e_upward_0*Z_e_0;
	V_i_e_downward_0 = -V_v_e_downward_0*Z_e_0;
	V_i_h_upward_0 = V_v_h_upward_0*Z_h_0;
	V_i_h_downward_0 = -V_v_h_downward_0*Z_h_0;

	V_i_e_upward_1 = V_v_e_upward_1*Z_e_1;
	V_i_e_downward_1 = -V_v_e_downward_1*Z_e_1;
	V_i_h_upward_1 = V_v_h_upward_1*Z_h_1;
	V_i_h_downward_1 = -V_v_h_downward_1*Z_h_1;

	V_i_e_upward_2 = V_v_e_upward_2*Z_e_2;
	V_i_e_downward_2 = -V_v_e_downward_2*Z_e_2;
	V_i_h_upward_2 = V_v_h_upward_2*Z_h_2;
	V_i_h_downward_2 = -V_v_h_downward_2*Z_h_2;

    //contributions to the potentials A_theta,A_phi,F_theta,F_phi by the equivalent currents J,M
    // upper half space
    // Green_A_theta_J = (-1/Z0)*G_{theta}^{J}
    Green_A_theta_J_x_upward_0 =  V_i_e_upward_0*cosP;
    Green_A_theta_J_x_downward_0 =  V_i_e_downward_0*cosP;
    Green_A_theta_J_y_upward_0 =  V_i_e_upward_0*sinP;
    Green_A_theta_J_y_downward_0 =  V_i_e_downward_0*sinP;
    Green_A_theta_J_z_upward_0 =  (-V_v_e_upward_0)*sinT/eps_r_0;
    Green_A_theta_J_z_downward_0 =  (-V_v_e_downward_0)*sinT/eps_r_0;
    // Green_A_phi_J = (-1/Z0)*G_{phi}^{J}
    Green_A_phi_J_x_upward_0 =  (-V_i_h_upward_0)*sinP*cosT;
    Green_A_phi_J_x_downward_0 =  (-V_i_h_downward_0)*sinP*cosT;
    Green_A_phi_J_y_upward_0 =  V_i_h_upward_0*cosP*cosT;
    Green_A_phi_J_y_downward_0 =  V_i_h_downward_0*cosP*cosT;
    // Green_F_theta_M = G_{phi}^{M}
    Green_F_theta_M_x_upward_0 =  V_v_h_upward_0*cosP*cosT;
    Green_F_theta_M_x_downward_0 =  V_v_h_downward_0*cosP*cosT;
    Green_F_theta_M_y_upward_0 =  V_v_h_upward_0*sinP*cosT;
    Green_F_theta_M_y_downward_0 =  V_v_h_downward_0*sinP*cosT;
    Green_F_theta_M_z_upward_0 =  (-V_i_h_upward_0)*sinT*cosT/m_r_0;
    Green_F_theta_M_z_downward_0 =  (-V_i_h_downward_0)*sinT*cosT/m_r_0;
    // Green_F_phi_M = -G_{theta}^{M}
    Green_F_phi_M_x_upward_0 =  (-V_v_e_upward_0)*sinP;
    Green_F_phi_M_x_downward_0 =  (-V_v_e_downward_0)*sinP;
    Green_F_phi_M_y_upward_0 =  V_v_e_upward_0*cosP;
    Green_F_phi_M_y_downward_0 =  V_v_e_downward_0*cosP;

    // middle layer
    // Green_A_theta_J = (-1/Z0)*G_{theta}^{J}
    Green_A_theta_J_x_upward_1 =  V_i_e_upward_1*cosP;
    Green_A_theta_J_x_downward_1 =  V_i_e_downward_1*cosP;
    Green_A_theta_J_y_upward_1 =  V_i_e_upward_1*sinP;
    Green_A_theta_J_y_downward_1 =  V_i_e_downward_1*sinP;
    Green_A_theta_J_z_upward_1 =  (-V_v_e_upward_1)*sinT/eps_r_1;
    Green_A_theta_J_z_downward_1 =  (-V_v_e_downward_1)*sinT/eps_r_1;
    // Green_A_phi_J = (-1/Z0)*G_{phi}^{J}
    Green_A_phi_J_x_upward_1 =  (-V_i_h_upward_1)*sinP*cosT;
    Green_A_phi_J_x_downward_1 =  (-V_i_h_downward_1)*sinP*cosT;
    Green_A_phi_J_y_upward_1 =  V_i_h_upward_1*cosP*cosT;
    Green_A_phi_J_y_downward_1 =  V_i_h_downward_1*cosP*cosT;
    // Green_F_theta_M = G_{phi}^{M}
    Green_F_theta_M_x_upward_1 =  V_v_h_upward_1*cosP*cosT;
    Green_F_theta_M_x_downward_1 =  V_v_h_downward_1*cosP*cosT;
    Green_F_theta_M_y_upward_1 =  V_v_h_upward_1*sinP*cosT;
    Green_F_theta_M_y_downward_1 =  V_v_h_downward_1*sinP*cosT;
    Green_F_theta_M_z_upward_1 =  (-V_i_h_upward_1)*sinT*cosT/m_r_1;
    Green_F_theta_M_z_downward_1 =  (-V_i_h_downward_1)*sinT*cosT/m_r_1;
    // Green_F_phi_M = -G_{theta}^{M}
    Green_F_phi_M_x_upward_1 =  (-V_v_e_upward_1)*sinP;
    Green_F_phi_M_x_downward_1 =  (-V_v_e_downward_1)*sinP;
    Green_F_phi_M_y_upward_1 =  V_v_e_upward_1*cosP;
    Green_F_phi_M_y_downward_1 =  V_v_e_downward_1*cosP;

    // lower half space
    // Green_A_theta_J = (-1/Z0)*G_{theta}^{J}
    Green_A_theta_J_x_upward_2 =  V_i_e_upward_2*cosP;
    Green_A_theta_J_x_downward_2 =  V_i_e_downward_2*cosP;
    Green_A_theta_J_y_upward_2 =  V_i_e_upward_2*sinP;
    Green_A_theta_J_y_downward_2 =  V_i_e_downward_2*sinP;
    Green_A_theta_J_z_upward_2 =  (-V_v_e_upward_2)*sinT/eps_r_2;
    Green_A_theta_J_z_downward_2 =  (-V_v_e_downward_2)*sinT/eps_r_2;
    // Green_A_phi_J = (-1/Z0)*G_{phi}^{J}
    Green_A_phi_J_x_upward_2 =  (-V_i_h_upward_2)*sinP*cosT;
    Green_A_phi_J_x_downward_2 =  (-V_i_h_downward_2)*sinP*cosT;
    Green_A_phi_J_y_upward_2 =  V_i_h_upward_2*cosP*cosT;
    Green_A_phi_J_y_downward_2 =  V_i_h_downward_2*cosP*cosT;
    // Green_F_theta_M = G_{phi}^{M}
    Green_F_theta_M_x_upward_2 =  V_v_h_upward_2*cosP*cosT;
    Green_F_theta_M_x_downward_2 =  V_v_h_downward_2*cosP*cosT;
    Green_F_theta_M_y_upward_2 =  V_v_h_upward_2*sinP*cosT;
    Green_F_theta_M_y_downward_2 =  V_v_h_downward_2*sinP*cosT;
    Green_F_theta_M_z_upward_2 =  (-V_i_h_upward_2)*sinT*cosT/m_r_2;
    Green_F_theta_M_z_downward_2 =  (-V_i_h_downward_2)*sinT*cosT/m_r_2;
    // Green_F_phi_M = -G_{theta}^{M}
    Green_F_phi_M_x_upward_2 =  (-V_v_e_upward_2)*sinP;
    Green_F_phi_M_x_downward_2 =  (-V_v_e_downward_2)*sinP;
    Green_F_phi_M_y_upward_2 =  V_v_e_upward_2*cosP;
    Green_F_phi_M_y_downward_2 =  V_v_e_downward_2*cosP;

    //some other temporary variables for fast post-processing
    //wavenumber in the observation half space, multiplied by the grid spacing
    ko_dx = kk_o_g*dx;
    //transmission-line wavenumbers in different layers, multiplied by the grid spacing
    kt0_dx = kt_g_0*dx;
    kt1_dx = kt_g_1*dx;
    kt2_dx = kt_g_2*dx;
}

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

// For the transmission-line Green functions, we use the phasor-domain equivalents of formulas (54)-(55) in Capoglu thesis pg. 28.
// From these, we construct the potential "dyadic Green functions", the "dyadic" components of which are represented by the arrays "Green_A_theta_J_x_minus_0" etc.
// (NOTE: They might not be "dyadic" functions in the exact sense of the word, but they operate on a field vector and yield a vector.)
// These "dyadic" functions are related to the functions G_{theta,phi}^{E,M} in formulas (21)-(24) in Capoglu thesis pg. 17.
//.In order to find these relations, we remember the expression for the radiated electric field in terms of the electric and magnetic potentials:
// (in three dimensions)
// According to Smith, "Classical EM Radiation" (pg. 215) :
// E_{theta or phi} = (jk)/(2*pi) * [(sin/cos of some angle) * (the component of the PW spectrum that points in the observation direction)]
// According to Harrington "Time-Harmonic EM Fields" (pg. 132) :
// E_{theta or phi} = (jk)/(4*pi) * (-mu^o*A -+ F), where A and F are the magnetic and electric vector potentials at the observation direction, mu^o is the permeability in the observation direction
//
// The expression in the square brackets in the Smith formulation is given by formulas (20)-(24) in Capoglu thesis pg. 17.
// In this code, we use the Harrington notation consistently; so it is useful to document the relations between these two formulations.
// By direct comparison between the Smith and Harrington formulations, one obtains
//   Green_A_theta_J = (-2/Z^o)*G_{theta}^{J}
//   Green_A_phi_J = (-2/Z^o)*G_{phi}^{J}
//   Green_F_theta_M = 2G_{phi}^{M}
//   Green_F_phi_M = -2G_{theta}^{M}
// In the above, Z^o is the wave impedance in the observation direction.
// The functions G_{theta,phi}^{E,M} above are given by formulas (21)-(24) in Capoglu thesis pg. 17.

#include "headers.h"

#include "Ctr_pd_ml.h"

extern double dx;

extern int number_of_layers;

extern Array<int,1> LayerLowerZIndices;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;

namespace{
	int pln,interf;
};

complex<double> stable_tan(const complex<double>& x);


void Ctr_pd_ml::Calculate_TL_GreenFunctions(const int& l, const int& d1, const int& d2)
{
    //temporary variables for fast post-processing
    //sines and cosines
    cosP = CosP(l,d1,d2);
    sinP = SinP(l,d1,d2);
    cosT = CosT(l,d1,d2);
    sinT = SinT(l,d1,d2);
    sinTcosP = sinT*cosP;
    sinTsinP = sinT*sinP;

	//the constitutive parameters in the observation half space (uppermost or lowermost)
	if (cosT>abs(threshold))
	{//upper half space
		epsilon_r_o = epsilon_r_upper;
		mu_r_o = mu_r_upper;
	}
	else if (cosT<-abs(threshold))
	{//lower half space
		epsilon_r_o = epsilon_r_lower;
		mu_r_o = mu_r_lower;
	}
	//wavenumber at the observation half space, corrected for grid-dispersion  (in m^-1)
	/** FIXME: I don't know why this works (or seems to work) for lossy observation space **/
	/** It seems as though kk_o_g (as well as epsilon_r_o) should be complex in that case **/
	kk_o_g = (ww(l)/c/GridVelocity(ww(l),THETA(l,d1,d2),PHI(l,d1,d2)))*sqrt(epsilon_r_o*mu_r_o);

	for (pln=0; pln<number_of_layers; pln++)
	{
		//(complex) rel. permittivities of the layers,  normalized by epsilon_r_o
		eps_r(pln) = (epsilon_r(pln)-ii*cond_e(pln)/epsilon_0/ww(l))/epsilon_r_o;
		m_r(pln) = mu_r(pln)/mu_r_o;		//rel. permeability of the upper half space, normalized by mu_r_o

		obliquityfactor = sqrt(complex<double>(eps_r(pln)*m_r(pln)-pow2(sinT)));
		// kt_g should have negative imaginary parts; therefore so should these obliquity factors.
		// Note that imaginary kt_g occurs when the source is radiating to the "opposite" half space. In this case, for example, the "downward" Green function contributes to the far field if the source is in the upper half space, and the obs. half space is the lower half space. In the reverse case, the "upward" Green function contributes. In both cases, it turns out that the longitudinal (z-dependent) phase term is in the form of exp(-ii*kt_g*dx*|z_pr|), where |z_pr| is the absolute value of the longitudinal distance from the interface. Therefore, kt_g should have negative imaginary part in order to make the term decay with |z_pr|.

		if (imag(obliquityfactor)>0)
		{
			obliquityfactor *= -1.0;
		}

		kt_g(pln) = obliquityfactor*kk_o_g;

		//Transmission-line impedances (normalized by Z, since only V_i_{e,h}/Z is used in the formulas)
		//These are defined in an indirect manner in terms of kt_g; the advantage being that the sign of the imaginary part of Z_{e,h} is automatically determined by that of kt_g. The choice of sign of the latter is explained above.
		//TM impedances
		eta_e(pln) = obliquityfactor/eps_r(pln);
		//TE impedances
		eta_h(pln) = pow(obliquityfactor,-1)*m_r(pln);
	}

	/** these are defined at interfaces between layers, therefore they have size (#-of-layers-1) **/
	if (cosT>abs(threshold))
	{//observation angle is in the upper half space
		//upward recursion for the impedances at the interfaces
		//initial values at the lowermost interface
		Z_e(LOWERMOST_INTERFACE_INDEX) = -eta_e(LOWERMOST_LAYER_INDEX);
		Z_h(LOWERMOST_INTERFACE_INDEX) = -eta_h(LOWERMOST_LAYER_INDEX);
		for (interf=1; interf<number_of_interfaces; interf++)
		{
			Z_e(interf) = eta_e(interf)*((Z_e(interf-1)-ii*eta_e(interf)*stable_tan(kt_g(interf)*d_n(interf)))
										/(eta_e(interf)-ii*Z_e(interf-1)*stable_tan(kt_g(interf)*d_n(interf))));
			Z_h(interf) = eta_h(interf)*((Z_h(interf-1)-ii*eta_h(interf)*stable_tan(kt_g(interf)*d_n(interf)))
										/(eta_h(interf)-ii*Z_h(interf-1)*stable_tan(kt_g(interf)*d_n(interf))));
		}
	}
	else if (cosT<-abs(threshold))
	{//observation angle is in the lower half space
		//downward recursion for the impedances at the interfaces
		//initial values at the uppermost interface
		Z_e(UPPERMOST_INTERFACE_INDEX) = eta_e(UPPERMOST_LAYER_INDEX);
		Z_h(UPPERMOST_INTERFACE_INDEX) = eta_h(UPPERMOST_LAYER_INDEX);
		for (interf=UPPERMOST_INTERFACE_INDEX-1; interf>=0; interf--)
		{
			Z_e(interf) = eta_e(interf+1)*((Z_e(interf+1)+ii*eta_e(interf+1)*stable_tan(kt_g(interf+1)*d_n(interf+1)))
										/(eta_e(interf+1)+ii*Z_e(interf+1)*stable_tan(kt_g(interf+1)*d_n(interf+1))));
			Z_h(interf) = eta_h(interf+1)*((Z_h(interf+1)+ii*eta_h(interf+1)*stable_tan(kt_g(interf+1)*d_n(interf+1)))
										/(eta_h(interf+1)+ii*Z_h(interf+1)*stable_tan(kt_g(interf+1)*d_n(interf+1))));
		}
	}
	/** these are defined at interfaces between layers, therefore they have size (#-of-layers-1) **/

	//TL Green functions V_i_e,V_i_h
	if (cosT>abs(threshold))
	{//observation angle is in the upper half space
		//downward recursion for V_i_*_minus
		V_i_e_minus(UPPERMOST_LAYER_INDEX) = eta_e(UPPERMOST_LAYER_INDEX)*exp(ii*kt_g(UPPERMOST_LAYER_INDEX)*dx*interface_coord(UPPERMOST_INTERFACE_INDEX+1));
		V_i_h_minus(UPPERMOST_LAYER_INDEX) = eta_h(UPPERMOST_LAYER_INDEX)*exp(ii*kt_g(UPPERMOST_LAYER_INDEX)*dx*interface_coord(UPPERMOST_INTERFACE_INDEX+1));
		for (pln=UPPERMOST_LAYER_INDEX-1; pln>=0; pln--)
		{
			V_i_e_minus(pln) = V_i_e_minus(pln+1)*(Z_e(pln)-eta_e(pln))/(Z_e(pln)-eta_e(pln+1))*exp(-ii*kt_g(pln+1)*d_n(pln+1));
			V_i_h_minus(pln) = V_i_h_minus(pln+1)*(Z_h(pln)-eta_h(pln))/(Z_h(pln)-eta_h(pln+1))*exp(-ii*kt_g(pln+1)*d_n(pln+1));
		}

		//V_i_*_plus is obtained from V_i_*_minus
		V_i_e_plus(LOWERMOST_LAYER_INDEX) = 0;
		V_i_h_plus(LOWERMOST_LAYER_INDEX) = 0;
		for (pln=LOWERMOST_LAYER_INDEX+1; pln<number_of_layers; pln++)
		{
			V_i_e_plus(pln) = V_i_e_minus(pln)*(Z_e(pln-1)+eta_e(pln))/(Z_e(pln-1)-eta_e(pln))*exp(-ii*kt_g(pln)*d_n(pln));
			V_i_h_plus(pln) = V_i_h_minus(pln)*(Z_h(pln-1)+eta_h(pln))/(Z_h(pln-1)-eta_h(pln))*exp(-ii*kt_g(pln)*d_n(pln));
		}
	}
	else if (cosT<-abs(threshold))
	{//observation angle is in the lower half space
		//upward recursion for V_i_*_plus
		V_i_e_plus(LOWERMOST_LAYER_INDEX) = eta_e(LOWERMOST_LAYER_INDEX)*exp(-ii*kt_g(LOWERMOST_LAYER_INDEX)*dx*interface_coord(LOWERMOST_INTERFACE_INDEX-1));
		V_i_h_plus(LOWERMOST_LAYER_INDEX) = eta_h(LOWERMOST_LAYER_INDEX)*exp(-ii*kt_g(LOWERMOST_LAYER_INDEX)*dx*interface_coord(LOWERMOST_INTERFACE_INDEX-1));
		for (pln=LOWERMOST_LAYER_INDEX+1; pln<number_of_layers; pln++)
		{
			V_i_e_plus(pln) = V_i_e_plus(pln-1)*(Z_e(pln-1)+eta_e(pln))/(Z_e(pln-1)+eta_e(pln-1))*exp(-ii*kt_g(pln-1)*d_n(pln-1));
			V_i_h_plus(pln) = V_i_h_plus(pln-1)*(Z_h(pln-1)+eta_h(pln))/(Z_h(pln-1)+eta_h(pln-1))*exp(-ii*kt_g(pln-1)*d_n(pln-1));
		}

		//V_i_*_minus is obtained from V_i_*_plus
		V_i_e_minus(UPPERMOST_LAYER_INDEX) = 0;
		V_i_h_minus(UPPERMOST_LAYER_INDEX) = 0;
		for (pln=UPPERMOST_LAYER_INDEX-1; pln>=0; pln--)
		{
			V_i_e_minus(pln) = V_i_e_plus(pln)*(Z_e(pln)-eta_e(pln))/(Z_e(pln)+eta_e(pln))*exp(-ii*kt_g(pln)*d_n(pln));
			V_i_h_minus(pln) = V_i_h_plus(pln)*(Z_h(pln)-eta_h(pln))/(Z_h(pln)+eta_h(pln))*exp(-ii*kt_g(pln)*d_n(pln));
		}
	}
	else
	{
		V_i_e_minus = 0;
		V_i_e_plus = 0;
		V_i_h_minus = 0;
		V_i_h_plus = 0;
	}

	//the TL Green's functions change sign if the obs. direction is in the lower half space
	if (cosT<-abs(threshold))
	{
		V_i_e_minus *= -1;
		V_i_e_plus *= -1;
		V_i_h_minus *= -1;
		V_i_h_plus *= -1;
	}

	//TL Green functions V_v_e,V_v_h
	V_v_e_minus = V_i_e_minus/eta_e;
	V_v_e_plus = -V_i_e_plus/eta_e;
	V_v_h_minus = V_i_h_minus/eta_h;
	V_v_h_plus = -V_i_h_plus/eta_h;

    //contributions to the potentials A_theta,A_phi,F_theta,F_phi by the equivalent currents J,M
    // upper half space
    // Green_A_theta_J = (-1/Z0)*G_{theta}^{J}
    Green_A_theta_J_x_minus =  V_i_e_minus*cosP;
    Green_A_theta_J_x_plus =  V_i_e_plus*cosP;
    Green_A_theta_J_y_minus =  V_i_e_minus*sinP;
    Green_A_theta_J_y_plus =  V_i_e_plus*sinP;
    Green_A_theta_J_z_minus =  (-V_v_e_minus)*sinT/eps_r;
    Green_A_theta_J_z_plus =  (-V_v_e_plus)*sinT/eps_r;
    // Green_A_phi_J = (-1/Z0)*G_{phi}^{J}
    Green_A_phi_J_x_minus =  (-V_i_h_minus)*sinP*cosT;
    Green_A_phi_J_x_plus =  (-V_i_h_plus)*sinP*cosT;
    Green_A_phi_J_y_minus =  V_i_h_minus*cosP*cosT;
    Green_A_phi_J_y_plus =  V_i_h_plus*cosP*cosT;
    // Green_F_theta_M = G_{phi}^{M}
    Green_F_theta_M_x_minus =  V_v_h_minus*cosP*cosT;
    Green_F_theta_M_x_plus =  V_v_h_plus*cosP*cosT;
    Green_F_theta_M_y_minus =  V_v_h_minus*sinP*cosT;
    Green_F_theta_M_y_plus =  V_v_h_plus*sinP*cosT;
    Green_F_theta_M_z_minus =  (-V_i_h_minus)*sinT*cosT/m_r;
    Green_F_theta_M_z_plus =  (-V_i_h_plus)*sinT*cosT/m_r;
    // Green_F_phi_M = -G_{theta}^{M}
    Green_F_phi_M_x_minus =  (-V_v_e_minus)*sinP;
    Green_F_phi_M_x_plus =  (-V_v_e_plus)*sinP;
    Green_F_phi_M_y_minus =  V_v_e_minus*cosP;
    Green_F_phi_M_y_plus =  V_v_e_plus*cosP;

    //some other temporary variables for fast post-processing
    //wavenumber in the observation half space, multiplied by the grid spacing
    ko_dx = kk_o_g*dx;
    //transmission-line wavenumbers in different layers, multiplied by the grid spacing
    kt_dx = kt_g*dx;

/** Remove later **/
int layer_index = UPPERMOST_LAYER_INDEX;
//////cout << layer_index << endl;
//cout << "downward:" << V_i_e_minus(layer_index) << "upward:" << V_i_e_plus(layer_index) << endl;
//cout << "downward:" << V_i_h_minus(layer_index) << "upward:" << V_i_h_plus(layer_index) << endl;
//cout << "downward:" << V_i_e_minus(layer_index-1) << "upward:" << V_i_e_plus(layer_index-1) << endl;
//cout << "downward:" << V_i_h_minus(layer_index-1) << "upward:" << V_i_h_plus(layer_index-1) << endl;
//cout << Green_A_theta_J_x_minus(layer_index) << endl; //V_i_e_minus*cosP;
//cout << Green_A_theta_J_x_plus(layer_index) << endl; //V_i_e_plus*cosP;
//cout << Green_A_theta_J_y_minus(layer_index) << endl; //V_i_e_minus*sinP;
//cout << Green_A_theta_J_y_plus(layer_index) << endl; //V_i_e_plus*sinP;
//cout << Green_A_theta_J_z_minus(layer_index) << endl; //(-V_v_e_minus)*sinT/eps_r;
//cout << Green_A_theta_J_z_plus(layer_index) << endl; //(-V_v_e_plus)*sinT/eps_r;
//cout << Green_A_phi_J_x_minus(layer_index) << endl; //(-V_i_h_minus)*sinP*cosT;
//cout << Green_A_phi_J_x_plus(layer_index) << endl; //(-V_i_h_plus)*sinP*cosT;
//cout << Green_A_phi_J_y_minus(layer_index) << endl; //V_i_h_minus*cosP*cosT;
//cout << Green_A_phi_J_y_plus(layer_index) << endl; //V_i_h_plus*cosP*cosT;
//cout << Green_F_theta_M_x_minus(layer_index) << endl; //V_v_h_minus*cosP*cosT;
//cout << Green_F_theta_M_x_plus(layer_index) << endl; //V_v_h_plus*cosP*cosT;
//cout << Green_F_theta_M_y_minus(layer_index) << endl; //V_v_h_minus*sinP*cosT;
//cout << Green_F_theta_M_y_plus(layer_index) << endl; //V_v_h_plus*sinP*cosT;
//cout << Green_F_theta_M_z_minus(layer_index) << endl; //(-V_i_h_minus)*sinT*cosT/m_r;
//cout << Green_F_theta_M_z_plus(layer_index) << endl; //(-V_i_h_plus)*sinT*cosT/m_r;
//cout << Green_F_phi_M_x_minus(layer_index) << endl; //(-V_v_e_minus)*sinP;
//cout << Green_F_phi_M_x_plus(layer_index) << endl; //(-V_v_e_plus)*sinP;
//cout << Green_F_phi_M_y_minus(layer_index) << endl; //V_v_e_minus*cosP;
//cout << Green_F_phi_M_y_plus(layer_index) << endl; //V_v_e_plus*cosP;
////cout << kt_dx << endl;
//exit(-1);
/** Remove later **/
}

inline complex<double> stable_tan(const complex<double>& x)
{
	return (sin(real(x))/(cos(real(x))-ii*sin(real(x))*tanh(imag(x)))
			+ii*cos(real(x))/(cos(real(x))/tanh(imag(x))-ii*sin(real(x))));
}

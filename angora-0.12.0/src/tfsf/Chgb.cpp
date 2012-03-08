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

//Defines the class "Chgb" for a TF/SF Hermite-Gaussian source in free space

#include "headers.h"

#include "Chgb.h"

#include "Cpw_fs.h"
#include "Cpw_2l.h"
#include "Cpw_ml.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

extern double dx;
extern int NCELLS_X,NCELLS_Y;

extern double epsilon_r_upper,mu_r_upper;

extern int rank;

extern int number_of_layers;

//gauss-legendre quadrature rule generator
extern void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w);

extern double Hermite(const int& N, const double& x);

extern void MPI_exit(const int& exitcode);

//factorial function (from Cgsmb)
extern int factorial (int num);


Chgb::Chgb(const HGBDataType& MyData, const int& Index)
		:Data(MyData), HGPIndex(Index)
{
	if ((Data.direction!="+z")&&(Data.direction!="-z"))
	/** FIXME: This check needs to be done by the HGBDataType class **/
	{
		if (rank==0)
		{
			cout << "Error: Invalid propagation direction (" << Data.direction << ") for the Hermite-Gaussian beam " << HGPIndex << endl << endl;
		}
		MPI_exit(-1);
	}

	//half width of the angular spectrum is assumed independent of frequency
	//consequently, the spatial beam width scales with the wavelength
	// "beam_halfwidth" is given for the center wavelength of the excitation
	double k_center = Data.waveform->W_0()/c;  //this is for free space
	//the beam is assumed to be incident from the upper half space, so use epsilon_r_upper,mu_r_upper to find the wavenumber in that material
	double k_times_beam_halfwidth = k_center*sqrt(epsilon_r_upper*mu_r_upper)*Data.beam_halfwidth;  //THIS IS ASSUMED TO BE CONSTANT FOR ALL WAVELENGTHS!
	double s_halfwidth = 1/k_times_beam_halfwidth; //half-width in the angle space

	//approximate maxima of the Hermite polynomials
	double hermite_max_x = sqrt(pow(2.0,Data.x_order)*factorial(Data.x_order));
	double hermite_max_y = sqrt(pow(2.0,Data.y_order)*factorial(Data.y_order));

	double angular_threshold = 1e-3;

	//the maximum direction cosine that needs to be considered (it could be at most 1)
	//do rough, brute-force calculation
	double sx_loop=(4+Data.x_order)*s_halfwidth;  //As the order increases, the decay is slower [4+Data.x_order is still a bit too far but..]
	double sx_step = 0.01*sx_loop; //step backward by %1 of the starting point
	double sy_loop=(4+Data.y_order)*s_halfwidth;  //As the order increases, the decay is slower [4+Data.y_order is still a bit too far but..]
	double sy_step = 0.01*sy_loop; //step backward by %1 of the starting point
	double ratio_to_maximum_x,ratio_to_maximum_y; //ratio of the spectrum at a certain point to the (roughly) maximum value
	do
	{//decrease sx until the threshold is found
		sx_loop-=sx_step;
		ratio_to_maximum_x = abs(Hermite(Data.x_order,sx_loop/s_halfwidth)*exp(-pow2(sx_loop/s_halfwidth)/2)/hermite_max_x);
	} while (ratio_to_maximum_x<angular_threshold);
	do
	{//decrease sy until the threshold is found
		sy_loop-=sy_step;
		ratio_to_maximum_y = abs(Hermite(Data.y_order,sy_loop/s_halfwidth)*exp(-pow2(sy_loop/s_halfwidth)/2)/hermite_max_y);
	} while (ratio_to_maximum_y<angular_threshold);

	double s_max = min(max(sx_loop,sy_loop),1.0);

	if (sx_loop>1.0)
	{
		if (rank==0)
		{
			cout << "Warning: The spatial spectrum of Hermite-Gaussian beam " << HGPIndex << " has non-propagating components (spectrum falls to " << abs(Hermite(Data.x_order,s_max/s_halfwidth)*exp(-pow2(s_max/s_halfwidth)/2)/hermite_max_x)*100 << "% of its maximum at the edge of the propagation circle in the x direction)" << endl << endl;
		}
	}
	if (sy_loop>1.0)
	{
		if (rank==0)
		{
			cout << "Warning: The spatial spectrum of Hermite-Gaussian beam " << HGPIndex << " has non-propagating components (spectrum falls to " << abs(Hermite(Data.y_order,s_max/s_halfwidth)*exp(-pow2(s_max/s_halfwidth)/2)/hermite_max_y)*100 << "% of its maximum at the edge of the propagation circle in the y direction)" << endl << endl;
		}
	}

	//now, for the maximum allowable direction cosine spacing:
	//the beam is assumed to be incident from the upper half space, so use epsilon_r_upper,mu_r_upper to find the real wavelengths
	//k_center and lambda_center are for free space:
//	double lambda_center = (2*M_PI/k_center)/sqrt(epsilon_r_upper*mu_r_upper);
	//minimum and maximum wavelengths in the excitation waveform
	// for lambda_min and lambda_max, -40 dB might be too tight, so we use -20dB
//	double lambda_min = c/sqrt(epsilon_r_upper*mu_r_upper)/(Data.waveform->w_max_20()/2/M_PI);
//	double lambda_max = c/sqrt(epsilon_r_upper*mu_r_upper)/(Data.waveform->w_min_20()/2/M_PI);
	//maximum k in the excitation waveform
	double k_max = (Data.waveform->w_max_20())/(c/sqrt(epsilon_r_upper*mu_r_upper));

	//maximum beam half-width among all the wavelengths (direction-cosine spacing is determined by this width)

	double minimum_necessary_shift_in_beam_halfwidth = 10;
	double k_times_beam_waist = minimum_necessary_shift_in_beam_halfwidth*k_times_beam_halfwidth;

	if (minimum_necessary_shift_in_beam_halfwidth<8)
		if (rank==0) cout << "Warning: There may be aliasing!!!!" << endl;

	//x and y extents are determined either by the beam waist or the size of the TF/SF box
	double TFSF_x_extent = (NCELLS_X-(Data.HGBMarginFrontX+Data.HGBMarginBackX))*dx;
	double TFSF_y_extent = (NCELLS_Y-(Data.HGBMarginRightY+Data.HGBMarginLeftY))*dx;
	double k_times_W_x = max(k_times_beam_waist,k_max*TFSF_x_extent);
	double k_times_W_y = max(k_times_beam_waist,k_max*TFSF_y_extent);

	//uninitialized in config file, use sampling theorem for default spacing
	dsx = 2*M_PI/k_times_W_x;
	//The range 0->(dsx*N_X1) is divided into N_X1 regions and sx is placed at the midpoint of each region.
	//This provides consistency with previous publications.
	//The point (sx,sy) on the 2D plane may fall out of the sin(th)<1 circle for N_X1=N_X2=2, but neither N_X1 nor N_X2 are supposed to be that small anyway. Note that dsx<=2*pi/k_times_beam_waist=2*pi/(minimum_necessary_shift_in_beam_halfwidth*k_times_beam_halfwidth)
	//						   <(lambda_center/beam_halfwidth)/(minimum_necessary_shift_in_beam_halfwidth)
	//						   <1/minimum_necessary_shift_in_beam_halfwidth
	// Therefore N_X1 is larger than (minimum_necessary_shift_in_beam_halfwidth), which is larger than 2.
	N_X1 = int(2*s_max/dsx)+1;

	//uninitialized in config file, use sampling theorem for default spacing
	dsy = 2*M_PI/k_times_W_y;
	//(see note above for sx)
	N_X2 = int(2*s_max/dsy)+1;

	Array<double,1> sx_array,sy_array;
	sx_array.resize(N_X1);
	sy_array.resize(N_X2);

	for (int i=0;i<N_X1;i++)
	{
		sx_array(i) = -s_max+(i+0.5)*dsx;	//evaluated at the center of each piece
	}

	for (int j=0;j<N_X2;j++)
	{
		sy_array(j) = -s_max+(j+0.5)*dsy;	//evaluated at the center of each piece
	}

	double sx,sy;
	//start adding the plane waves
	for (int i=0;i<N_X1;i++)
	{
		sx = sx_array(i);
		for (int j=0;j<N_X2;j++)
		{
			sy = sy_array(j);

			if ((pow2(sx)+pow2(sy))<=1)  //is there a direction corresponding to sx,sy?
			{
				/** determine the theta, phi and psi angles of the individual plane wave **/
				/** NOTE: sx, sy are just angular Fourier parameters. The actual direction cosines of the *incident* plane wave might be different. **/
				//Let's find the direction cosines of the PWs:
				//the plane wave is INCIDENT from the opposite direction in which it PROPAGATES:
				double sx_PW = -sx;
				double sy_PW = -sy;
				//abs_sz is always positive: The actual sign of sz depends of the "direction" parameter
				double abs_sz_PW = real(sqrt((complex<double>)(1-pow2(sx)-pow2(sy)))); //cast into complex, in case there are roundoff errors for sx^2+sy^2=1
				double sz_PW;
				if (Data.direction=="+z")
				{
					sz_PW = -abs_sz_PW;	//the plane wave should be incident from the lower half space
				}
				else if (Data.direction=="-z")
				{
					sz_PW = abs_sz_PW;	//the plane wave should be incident from the upper half space
				}

				theta = acos(sz_PW);	//sz = cos(theta), result of acos is always in the [0,pi] range

				if ((abs(sx_PW)<LIBSTD_DBL_EPSILON*100)&&(abs(sy_PW)<LIBSTD_DBL_EPSILON*100)) //is theta=0?
				{
					phi = 0; //just assign an arbitrary value, since phi is undefined here
				}
				else
				{
					phi = atan2(sy_PW,sx_PW);
				}
				//atan2 is between [-pi,pi]
				//if less than 0, make phi positive
				if (phi<0)
				{
					phi += 2*M_PI;
				}

				//psi is adjusted to give zero cross-pol component
				double relative_phi = phi - Data.polarization; //phi angle relative to the polarization
				psi = M_PI + atan2(cos(relative_phi),cos(theta)*sin(relative_phi)); //output of atan2 is between [-pi,pi]
				/** theta, phi and psi angles of the individual plane wave are determined. **/

				PWDataType PWData;
				PWData.THETA = theta;
				PWData.PHI = phi;
				PWData.PSI = psi;

				PWData.PWMarginBackX = Data.HGBMarginBackX;
				PWData.PWMarginFrontX = Data.HGBMarginFrontX;
				PWData.PWMarginLeftY = Data.HGBMarginLeftY;
				PWData.PWMarginRightY = Data.HGBMarginRightY;
				PWData.PWMarginLowerZ = Data.HGBMarginLowerZ;
				PWData.PWMarginUpperZ = Data.HGBMarginUpperZ;
				PWData.PWOriginX = Data.HGBOriginX;
				PWData.PWOriginY = Data.HGBOriginY;
				PWData.PWOriginZ = Data.HGBOriginZ;
				PWData.L_req = Data.L_req;

				PWData.DisplayWarnings = Data.DisplayWarnings;
				PWData.DisplayCellsPerLambda = Data.DisplayCellsPerLambda;

				PWData.E0 = Data.E0*dsx*dsy*
						sqrt(pow2(sin(relative_phi))+pow2(cos(relative_phi)/cos(theta)))* //amplitude of pw component such that no cross-pol component exists
						pow2(k_center/(2*M_PI))* //k^2 is to convert the angular spectrum into inverse (spatial) Fourier transform, (2pi)^2 is the 2D inverse-Fourier transform coefficient
						(2*M_PI*pow2(Data.beam_halfwidth))*  //(pi/c) factor in the angular spectrum [the i^{m+n} factor is taken care of below]
						Hermite(Data.x_order,sx/s_halfwidth)*exp(-pow2(sx/s_halfwidth)/2)*
						Hermite(Data.y_order,sy/s_halfwidth)*exp(-pow2(sy/s_halfwidth)/2);
//cout << Data.E0 << endl;
//cout << dsx*dsy << endl;
//cout << sqrt(pow2(cos(theta))*pow2(sin(relative_phi))+pow2(cos(relative_phi))) << endl;
//cout << (k_center*k_center) << endl;
//cout << (-1/pow2(2*M_PI)/cos(theta)) << endl;
//cout << (2*M_PI*pow2(Data.beam_halfwidth)) << endl;
//cout << Hermite(Data.x_order,sx/s_halfwidth) << endl;
//cout << exp(-pow2(sx/s_halfwidth)/2) << endl;
//cout << Hermite(Data.y_order,sy/s_halfwidth) << endl;
//cout << exp(-pow2(sy/s_halfwidth)/2) << endl;
//cout << PWData.E0 << endl;
//exit(-1);

				//waveform of the constituent plane waves is the integral of the image-plane waveform
				// every (i) multiplier in the i^{m+n} factor corresponds to a 90deg phase shift (-1 times the Hilbert transform)
				switch((Data.x_order+Data.y_order) % 4) {
				case  0  : PWData.waveform = Data.waveform;
						break;
				case  1  : PWData.waveform = Data.waveform->HilbertTransform();
						PWData.E0 *= -1;
						break;
				case  2  : PWData.waveform = Data.waveform;
						PWData.E0 *= -1;
						break;
				case  3  : PWData.waveform = Data.waveform->HilbertTransform();
						break;
				}

				//Add the plane wave source to the Hermite-Gaussian-beam object
				if (number_of_layers==1)
				{//attach a free-space plane-wave source
					boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_fs(PWData));
					PlaneWaves.push_back(new_pw_ptr);
				}
				else if (number_of_layers==2)
				{//attach a 2-layered-medium plane-wave source
					boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_2l(PWData));
					PlaneWaves.push_back(new_pw_ptr);
				}
				else
				{//attach a multilayered-medium plane-wave source
					boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_ml(PWData));
					PlaneWaves.push_back(new_pw_ptr);
				}
			}
		}
	}
}

void Chgb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the Hermite-Gaussian beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Chgb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the Hermite-Gaussian beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Chgb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
}

void Chgb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Chgb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
}

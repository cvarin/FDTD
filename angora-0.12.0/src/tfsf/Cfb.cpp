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

//Defines the class "Cfb" for a TF/SF focused-beam source in free space

#include "headers.h"

#include "Cfb.h"

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


Cfb::Cfb(const FBDataType& MyData, const int& Index)
		:Data(MyData), FPIndex(Index)
{
	if ((cos(Data.theta_max)<0)&&(sin(Data.theta_max)<0))
	{
		if (rank==0)
		{
			cout << "Error: Theta half-range (theta_max) should be between 0 and 90deg." << endl;
		}
		exit(-1);
	}

	//maximum x and y extents of the TF/SF box
	// for lambda_min and lambda_max, -40 dB might be too tight, so we use -20dB
	double lambda_min = c/sqrt(epsilon_r_upper*mu_r_upper)/(Data.waveform->w_max_20()/2/M_PI);	//the beam is assumed to be incident from the upper half space, so use epsilon_r_upper,mu_r_upper
	double lambda_max = c/sqrt(epsilon_r_upper*mu_r_upper)/(Data.waveform->w_min_20()/2/M_PI);	//the beam is assumed to be incident from the upper half space, so use epsilon_r_upper,mu_r_upper

	double minimum_aliasing_distance_in_beam_width = 5; //must be larger than 2
	double beam_waist = minimum_aliasing_distance_in_beam_width*lambda_max/abs(sin(Data.theta_max));  //FIXME:should find a better way for this later

	if (minimum_aliasing_distance_in_beam_width<5)
		if (rank==0) cout << "Warning: There may be aliasing!!!!" << endl;

	//x and y extents are determined either by the beam waist or the size of the TF/SF box
	double W_x = max(beam_waist,(NCELLS_X-(Data.FBMarginFrontX+Data.FBMarginBackX))*dx);
	double W_y = max(beam_waist,(NCELLS_Y-(Data.FBMarginRightY+Data.FBMarginLeftY))*dx);

	if (Data.angular_discretization=="cartesian")
	{
		if (Data.n_1==-1)
		{
			//uninitialized in config file, use sampling theorem for default spacing
			dsx = lambda_min/W_x;
			//The range 0->(dsx*N_X1) is divided into N_X1 regions and sx is placed at the midpoint of each region.
			//This provides consistency with previous publications.
			//The point (sx,sy) on the 2D plane may fall out of the sin(th)<1 circle for N_X1=N_X2=2, but neither N_X1 nor N_X2 are supposed to be that small anyway. Note that dsx<=lambda_min/beam_waist=sin(th_max)/(minimum_aliasing_distance_in_beam_width)*(lambda_min/lambda_max)
			//						   <sin(th_max)/(minimum_aliasing_distance_in_beam_width)
			// Therefore N_X1 is larger than (minimum_aliasing_distance_in_beam_width), which is larger than 2.
			N_X1 = int(2*abs(sin(Data.theta_max))/dsx)+1;
		}
		else
		{
			N_X1 = Data.n_1;
			if (N_X1<=0)
			{
				if (rank==0)
				{
					cout << "Error: Invalid number of x-direction cosines (" << N_X1 << ")" << endl;
				}
				exit(-1);
			}
			dsx = 2*abs(sin(Data.theta_max))/N_X1;
		}
		if (Data.n_2==-1)
		{
			//uninitialized in config file, use sampling theorem for default spacing
			dsy = lambda_min/W_y;
			//(see note above for sx)
			N_X2 = int(2*abs(sin(Data.theta_max))/dsy)+1;
		}
		else
		{
			N_X2 = Data.n_2;
			if (N_X2<=0)
			{
				if (rank==0)
				{
					cout << "Error: Invalid number of y-direction cosines (" << N_X2 << ")" << endl;
				}
				exit(-1);
			}
			dsy = 2*abs(sin(Data.theta_max))/N_X2;
		}

		//form the 2D theta and phi arrays for the incidence directions
		angle_within_ill_cone.resize(N_X1,N_X2);
		angle_within_ill_cone = true;
		theta_array.resize(N_X1,N_X2);
		phi_array.resize(N_X1,N_X2);
		pwfactor.resize(N_X1,N_X2);

		for (int sx_index=0; sx_index<N_X1; sx_index++)
		{
			for (int sy_index=0; sy_index<N_X2; sy_index++)
			{
				//(see note above for the placement of sx and sy)
				sx = dsx*(sx_index-(N_X1-1.0)/2.0);
				sy = dsy*(sy_index-(N_X2-1.0)/2.0);
				if (((pow2(sx)+pow2(sy))>1) //dir1^2+dir2^2 = sin^2(theta) should be less than or equal to 1 for a valid sin(theta)
				   ||((pow2(sx)+pow2(sy))>pow2(abs(sin(Data.theta_max)))))	// it should also be less than [maximum allowable direction cosine]^2
				{
					angle_within_ill_cone(sx_index,sy_index) = false;
				}
				else
				{
					theta = asin(sqrt(pow2(sx)+pow2(sy)));	//sx^2+sy^2 = sin^2(theta)
						 //^ select the [0,pi/2] range
					if ((abs(sx)<LIBSTD_DBL_EPSILON*100)&&(abs(sy)<LIBSTD_DBL_EPSILON*100)) //is THETA=0?
					{
						phi = 0; //just assign an arbitrary value, since phi is undefined here
					}
					else
					{
						phi = atan2(sy,sx);
					}
					//atan2 is between [-pi,pi]
					//if less than 0, make phi positive
					if (phi<0)
					{
						phi += 2*M_PI;
					}
					theta_array(sx_index,sy_index) = theta;
					phi_array(sx_index,sy_index) = phi;
				}
			}
		}

		//amplitude factor for every plane wave
		pwfactor = dsx*dsy;
	}
	else if (Data.angular_discretization=="radial")
	{
		// Data.n_1 and Data.n_2 *must* be initialized: This should be ensured in read_tfsf
		N_X1 = Data.n_1;
		N_X2 = Data.n_2;
		//spacing of phi (between [0,pi])
		d_phi = M_PI/N_X2;

		//Gauss-Legendre parameters for rho quadrature [ rho=sin(theta), -sin(theta_max)<rho<sin(theta_max) ]
		GLx.resize(N_X1); //Gauss-Legendre coordinates
		GLw.resize(N_X1); //Gauss-Legendre weights

		//calculate GL quadrature rule for the sin(theta) values (between [-sin(theta_max),sin(theta_max)])
		gaussquadrule(N_X1,GLx,GLw);
		//rescale the positions and weights to the correct ranges of sin(theta) (between [-sin(theta_max),sin(theta_max)])
		//see http://en.wikipedia.org/wiki/Gaussian_quadrature for the following formulas
		//GLx and GLw are for sin(theta)
		GLx = abs(sin(Data.theta_max))*(-GLx);	//-GLx because it is calculated in reverse direction in gaussquadrule
		GLw = abs(sin(Data.theta_max))*GLw;

		//form the 2D theta and phi arrays for the incidence directions
		angle_within_ill_cone.resize(N_X1,N_X2); //always true for "radial"
		angle_within_ill_cone = true;
		theta_array.resize(N_X1,N_X2);
		phi_array.resize(N_X1,N_X2);
		pwfactor.resize(N_X1,N_X2);

		for(int rho_index=0; rho_index<N_X1; rho_index++)
		{
			for(int phi_index=0; phi_index<N_X2; phi_index++)
			{
				theta = asin(GLx(rho_index)); //GLx and GLw are for sin(theta), result is between [-pi/2,pi/2]
				phi = (phi_index+0.5)*d_phi;

				if (theta<0)
				{// convert to traditional spherical angles
					theta *= -1;
					phi += M_PI;
				}
				theta_array(rho_index,phi_index) = theta;
				phi_array(rho_index,phi_index) = phi;
				//amplitude factor for every plane wave
				pwfactor(rho_index,phi_index) = abs(sin(theta))*GLw(rho_index)*d_phi; //could also use GLx(rho_index) for sin(theta)
			}
		}
	}
	else
	{
		if (rank==0)
		{
			cout << "Developer error: Angular discretization type (angular_discretization) should be either \"cartesian\" or \"radial\"" << endl;
		}
		exit(-1);
	}

// // cout << theta_array << endl;
// // cout << phi_array << endl;
// ofstream testfile;
// testfile.open("test.bin",ios::binary);
// testfile.write((char*)&N_X1,sizeof(N_X1));
// testfile.write((char*)&N_X2,sizeof(N_X2));
// testfile.write((char*)&dsx,sizeof(dsx));
// testfile.write((char*)&dsy,sizeof(dsy));
// testfile.write((char*)angle_within_ill_cone.data(),angle_within_ill_cone.size()*sizeof(angle_within_ill_cone(0)));
// testfile.write((char*)theta_array.data(),theta_array.size()*sizeof(theta_array(0)));
// testfile.write((char*)phi_array.data(),phi_array.size()*sizeof(phi_array(0)));
// testfile.write((char*)pwfactor.data(),pwfactor.size()*sizeof(pwfactor(0)));
// // testfile << phi_array << endl;
// testfile.close();

	//start adding the plane waves
	for (int i=0;i<N_X1;i++)
	{
		for (int j=0;j<N_X2;j++)
		{
			if (angle_within_ill_cone(i,j))
			{//do not add PW if it is outside illumination cone
				theta = theta_array(i,j);
				phi = phi_array(i,j);

				psi = Data.pw_pol + 3*M_PI/2 - phi; //pw_pol is ccw wrt the +x axis, like the spherical phi angle

				PWDataType PWData;
				PWData.THETA = theta;
				PWData.PHI = phi;
				PWData.PSI = psi;

				PWData.PWMarginBackX = Data.FBMarginBackX;
				PWData.PWMarginFrontX = Data.FBMarginFrontX;
				PWData.PWMarginLeftY = Data.FBMarginLeftY;
				PWData.PWMarginRightY = Data.FBMarginRightY;
				PWData.PWMarginLowerZ = Data.FBMarginLowerZ;
				PWData.PWMarginUpperZ = Data.FBMarginUpperZ;
				PWData.PWOriginX = Data.FBOriginX;
				PWData.PWOriginY = Data.FBOriginY;
				PWData.PWOriginZ = Data.FBOriginZ;
				PWData.L_req = Data.L_req;

				PWData.DisplayWarnings = Data.DisplayWarnings;
				PWData.DisplayCellsPerLambda = Data.DisplayCellsPerLambda;

				PWData.E0 = Data.E0*pwfactor(i,j)*(Data.f/(2*M_PI*c))*sqrt(abs(cos(theta)));
																		//cos^1/2(theta) from Richards&Wolf
																		//sin(theta) is in the integrand (2.26 in Richards&Wolf)
																		//abs(sin(theta)) because sin(theta) can go negative in radial quadrature
				//waveform of the constituent plane waves is the derivative of the incident plane-wave waveform
				PWData.waveform = Data.waveform->Derivative();

				//Add the plane wave source to the focused-beam object
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

void Cfb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the focused beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Cfb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the focused beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Cfb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
}

void Cfb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Cfb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
}

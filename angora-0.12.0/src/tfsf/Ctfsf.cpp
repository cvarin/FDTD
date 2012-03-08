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

//Defines the TF/SF object "Ctfsf"
//An TF/SF object may comprise multiple incident wave sources.

#include "headers.h"

#include "Ctfsf.h"

#include "Cpw_fs.h"
#include "Cpw_2l.h"
#include "Cpw_ml.h"
#include "Cfb.h"
#include "Chgb.h"
#include "Cgsmb.h"
#include "Ckohler.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern int GridIndex;

extern int number_of_layers;

extern int rank;

extern void add_slash_to_path(string& path);

extern int create_path(const string& path);


int Ctfsf::AddPlaneWave(const PWDataType& MyData)
{//Adds a plane wave source to the TF/SF object

	if (number_of_layers==1)
	{//attach a free-space plane-wave source
		boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_fs(MyData));
		PlaneWaves.push_back(new_pw_ptr);
	}
	else if (number_of_layers==2)
	{//attach a 2-layered-medium plane-wave source
		try {
		boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_2l(MyData));
		PlaneWaves.push_back(new_pw_ptr);
		}
		catch (EvanescentPWException& evanpw)
		{// catch exception if there is an evanescent PW in the lower layer
			if (rank==0)
			{
				cout << evanpw.what() << endl;
			}
			//insert a Cpw_ml object instead, since they can handle evanescent PWs
			boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_ml(MyData));
			PlaneWaves.push_back(new_pw_ptr);
		}
	}
	else
	{//attach a multilayered-medium plane-wave source
		boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_ml(MyData));
		PlaneWaves.push_back(new_pw_ptr);
	}
	return NumberOfPlaneWaves()-1;
}

int Ctfsf::AddFocusedBeam(const FBDataType& MyData)
{//adds a focused beam to the TF/SF object
	boost::shared_ptr<Cfb> new_fp_ptr(new Cfb(MyData,NumberOfFocusedBeams()));
	FocusedBeams.push_back(new_fp_ptr);
	return NumberOfFocusedBeams()-1;
}

int Ctfsf::AddHermiteGaussianBeam(const HGBDataType& MyData)
{//adds a Hermite-Gaussian beam to the TF/SF object
	boost::shared_ptr<Chgb> new_hgp_ptr(new Chgb(MyData,NumberOfHermiteGaussianBeams()));
	HermiteGaussianBeams.push_back(new_hgp_ptr);
	return NumberOfHermiteGaussianBeams()-1;
}

int Ctfsf::AddGaussianSchellModelBeam(const GSMBDataType& MyData)
{//adds a Gaussian Schell-model beam to the TF/SF object
	boost::shared_ptr<Cgsmb> new_gsmp_ptr(new Cgsmb(MyData,NumberOfGaussianSchellModelBeams()));
	GaussianSchellModelBeams.push_back(new_gsmp_ptr);
	return NumberOfGaussianSchellModelBeams()-1;
}

int Ctfsf::AddKohlerBeam(const KBDataType& MyData)
{//adds Kohler beam to the TF/SF object
	boost::shared_ptr<Ckohler> new_ki_ptr(new Ckohler(MyData,NumberOfKohlerBeams()));
	KohlerBeams.push_back(new_ki_ptr);
	return NumberOfKohlerBeams()-1;
}

void Ctfsf::CorrectE(const int& n)
{
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
	for (int i=0; i<NumberOfFocusedBeams(); i++)
	{
		FocusedBeams[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
	for (int i=0; i<NumberOfGaussianSchellModelBeams(); i++)
	{
		GaussianSchellModelBeams[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
	for (int i=0; i<NumberOfKohlerBeams(); i++)
	{
		KohlerBeams[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Ctfsf::CorrectH(const int& n)
{
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
	for (int i=0; i<NumberOfFocusedBeams(); i++)
	{
		FocusedBeams[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
	for (int i=0; i<NumberOfGaussianSchellModelBeams(); i++)
	{
		GaussianSchellModelBeams[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
	for (int i=0; i<NumberOfKohlerBeams(); i++)
	{
		KohlerBeams[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Ctfsf::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
	for (int i=0; i<NumberOfFocusedBeams(); i++)
	{
		FocusedBeams[i]->WriteScatteredPWDirections(PW_THETA,PW_PHI);
	}
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWDirections(PW_THETA,PW_PHI);
	}
	for (int i=0; i<NumberOfGaussianSchellModelBeams(); i++)
	{
		GaussianSchellModelBeams[i]->WriteScatteredPWDirections(PW_THETA,PW_PHI);
	}
	for (int i=0; i<NumberOfKohlerBeams(); i++)
	{
		KohlerBeams[i]->WriteScatteredPWDirections(PW_THETA,PW_PHI);
	}
}

void Ctfsf::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
	for (int i=0; i<NumberOfFocusedBeams(); i++)
	{
		FocusedBeams[i]->WriteScatteredPWDelaysFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWDelaysFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
	for (int i=0; i<NumberOfGaussianSchellModelBeams(); i++)
	{
		GaussianSchellModelBeams[i]->WriteScatteredPWDelaysFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
	for (int i=0; i<NumberOfKohlerBeams(); i++)
	{
		KohlerBeams[i]->WriteScatteredPWDelaysFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Ctfsf::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
	for (int i=0; i<NumberOfFocusedBeams(); i++)
	{
		FocusedBeams[i]->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);
	}
	for (int i=0; i<NumberOfHermiteGaussianBeams(); i++)
	{
		HermiteGaussianBeams[i]->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);
	}
	for (int i=0; i<NumberOfGaussianSchellModelBeams(); i++)
	{
		GaussianSchellModelBeams[i]->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);
	}
	for (int i=0; i<NumberOfKohlerBeams(); i++)
	{
		KohlerBeams[i]->WriteScatteredPWFieldAmplitudes(E_x_array,E_y_array);
	}
}

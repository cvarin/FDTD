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

//Definitions for the class "Cimg" for synthesizing an optical image
//This file includes the implementation of the post-processing steps for the construction of the final image.

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

#include "Cimg.h"

//definition of Ctr_pd needed
#include "nffft/pd/Ctr_pd.h"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//include the HDF5 header, if not disabled
#ifndef HDF5_DISABLE
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
#endif

//fft routines
#include "fft/fft.h"

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern int rank;

extern void MPI_exit(const int& exitcode);

//// the minimum package version that the far-field file should have
//const int minimum_package_version = 400;


void Cimg::ConstructImage()
{
	//first, do the post-processing steps to construct the far field
	Transformer->ConstructFarField();

#ifndef MPI_DISABLE
	MPI_Barrier(MPI_CartSubComm);
#endif

	if (rank==0)
	{// only node 0 constructs the image
#ifndef HDF5_DISABLE
		//open the far-field HDF5 file with read-only access
		H5File far_field_file(FarFieldOutputFileName.c_str(), H5F_ACC_RDONLY);
		DataSet dataset;
		/** read whether the direction cosines are scaled with wavelength **/
		bool scale_with_wavelength;
		dataset = far_field_file.openDataSet("scale_with_wavelength");
		dataset.read(&scale_with_wavelength, PredType::NATIVE_CHAR);
		//if scale_with_wavelength is false, there is something wrong
		if (!scale_with_wavelength)
		{
		    throw AngoraDeveloperException("The direction cosines should be scaled by wavelength.");
		}
		/** read number of wavelengths **/
		int lambda_size;
		dataset = far_field_file.openDataSet("num_lambda");
		dataset.read(&lambda_size, PredType::NATIVE_INT);
		/** read number of x direction cosines **/
		int dir1_size;
		try{
		dataset = far_field_file.openDataSet("num_dircos_x");
		}
		catch (FileIException& exc)
		{
		    throw AngoraDeveloperException(exc.getDetailMsg());
		}
		dataset.read(&dir1_size, PredType::NATIVE_INT);
		/** read number of y direction cosines **/
		int dir2_size;
		try{
		dataset = far_field_file.openDataSet("num_dircos_y");
		}
		catch (FileIException& exc)
		{
		    throw AngoraDeveloperException(exc.getDetailMsg());
		}
		dataset.read(&dir2_size, PredType::NATIVE_INT);
		/** read wavelength array **/
		Array<double,1> lambda(lambda_size);
		dataset = far_field_file.openDataSet("lambda");
		dataset.read(lambda.data(), PredType::NATIVE_DOUBLE);
		/** read x direction cosines **/
		Array<double,1> dir1(dir1_size);
		try{
		dataset = far_field_file.openDataSet("dircos_x");
		}
		catch (FileIException& exc)
		{
		    throw AngoraDeveloperException(exc.getDetailMsg());
		}
		dataset.read(dir1.data(), PredType::NATIVE_DOUBLE);
		/** read y direction cosines **/
		Array<double,1> dir2(dir2_size);
		try{
		dataset = far_field_file.openDataSet("dircos_y");
		}
		catch (FileIException& exc)
		{
		    throw AngoraDeveloperException(exc.getDetailMsg());
		}
		dataset.read(dir2.data(), PredType::NATIVE_DOUBLE);

		/** read number of plane waves reaching the image space **/
		int numofPWs;
		dataset = far_field_file.openDataSet("num_pw");
		dataset.read(&numofPWs, PredType::NATIVE_INT);
		/** read the theta and pghi angles of the plane waves **/
		Array<double,1> PW_theta(numofPWs);
		Array<double,1> PW_phi(numofPWs);
		if (numofPWs>0)
		{//don't read if there are no plane waves
			dataset = far_field_file.openDataSet("pw_theta");
			dataset.read(PW_theta.data(), PredType::NATIVE_DOUBLE);
			dataset = far_field_file.openDataSet("pw_phi");
			dataset.read(PW_phi.data(), PredType::NATIVE_DOUBLE);
		}
		/** read E_theta and E_phi **/
		Array<complex<double>,3> E_theta(lambda_size,dir1_size,dir2_size);
		Array<complex<double>,3> E_phi(lambda_size,dir1_size,dir2_size);
		Array<double,3> temp(lambda_size, dir1_size, dir2_size);
		dataset = far_field_file.openDataSet("E_theta_r");
		dataset.read(temp.data(), PredType::NATIVE_DOUBLE);
		real(E_theta) = temp;
		dataset = far_field_file.openDataSet("E_theta_i");
		dataset.read(temp.data(), PredType::NATIVE_DOUBLE);
		imag(E_theta) = temp;
		dataset = far_field_file.openDataSet("E_phi_r");
		dataset.read(temp.data(), PredType::NATIVE_DOUBLE);
		real(E_phi) = temp;
		dataset = far_field_file.openDataSet("E_phi_i");
		dataset.read(temp.data(), PredType::NATIVE_DOUBLE);
		imag(E_phi) = temp;
		/** read the plane-wave phasors **/
		Array<complex<double>,2> PW_x_phasorarray(numofPWs,lambda_size);
		Array<complex<double>,2> PW_y_phasorarray(numofPWs,lambda_size);
		Array<double,2> temp2D(numofPWs, lambda_size);
		if (numofPWs>0)
		{//don't read if there are no plane waves
			dataset = far_field_file.openDataSet("E_x_pw_r");
			dataset.read(temp2D.data(), PredType::NATIVE_DOUBLE);
			real(PW_x_phasorarray) = temp2D;
			dataset = far_field_file.openDataSet("E_x_pw_i");
			dataset.read(temp2D.data(), PredType::NATIVE_DOUBLE);
			imag(PW_x_phasorarray) = temp2D;
			dataset = far_field_file.openDataSet("E_y_pw_r");
			dataset.read(temp2D.data(), PredType::NATIVE_DOUBLE);
			real(PW_y_phasorarray) = temp2D;
			dataset = far_field_file.openDataSet("E_y_pw_i");
			dataset.read(temp2D.data(), PredType::NATIVE_DOUBLE);
			imag(PW_y_phasorarray) = temp2D;
		}
		//close the far-field file
		far_field_file.close();

		/*************************/
		/**** COLLECTION PART ****/
		/*************************/

		//No need to check for NA, since it is assumed that the coll. NA is larger than the img. NA

		/*************************/
		/**** REFOCUSING PART ****/
		/*************************/

		Array<double,1> kk(nimg*2*M_PI/lambda);

		if (kk.size()>2)
		{
			if (abs(((kk(1)-kk(0))-(kk(2)-kk(1)))/(kk(1)-kk(0)))>1e-5)
			{
				if (rank==0) cout << "Warning: The wavenumber (k) vector is not linearly-spaced." << endl;
			}
		}

        double min_lambda = min(abs(lambda)); //minimum lambda
        double kmax = max(abs(kk)); //wavenumber (in the image space) corresponding to min_lambda

		Array<double,1> ones1(dir1.shape());
		Array<double,1> ones2(dir2.shape());
		ones1=1;ones2=1;
		Array<double,3> D1(lambda.size(),dir1.size(),dir2.size()),
                        D2(lambda.size(),dir1.size(),dir2.size()),
                        D3(lambda.size(),dir1.size(),dir2.size());
		for (int k=0; k<lambda_size; k++)
		{
            D1(k,Range::all(),Range::all()) = (lambda(k)/min_lambda)*dir1(tensor::i)*ones2(tensor::j);	//outer product, scaled by lambda
            D2(k,Range::all(),Range::all()) = (lambda(k)/min_lambda)*ones1(tensor::i)*dir2(tensor::j);	//outer product, scaled by lambda
		}
		if (Data.coll_half_space=="upper")
		{
			D3 = real(sqrt((Array<complex<double>,3>)(1-pow2(D1)-pow2(D2)))); //cast into complex, in case d1^2+d2^2>1
		}
		else if (Data.coll_half_space=="lower")
		{
			D3 = -real(sqrt((Array<complex<double>,3>)(1-pow2(D1)-pow2(D2)))); //cast into complex, in case d1^2+d2^2>1
		}

		Array<double,3> CosThetaObj(D3.shape()),invCosThetaObj(D3.shape());
		CosThetaObj = D3;
		invCosThetaObj = 0; //initialize to zero
		for (int k=0; k<lambda_size; k++)
		{
            for (int d1=0;d1<dir1_size;d1++)
            {
                for (int d2=0;d2<dir2_size;d2++)
                {
                    if (CosThetaObj(k,d1,d2)!=0)
                    {
                        invCosThetaObj(k,d1,d2) = 1/CosThetaObj(k,d1,d2);
                    }//for non-propagating angles, leave invCosThetaObj at zero
                }
            }
		}

		//direction cosines in the image space
		Array<double,1> dir1M(dir1.shape()),dir2M(dir2.shape());
		//direction cosines in the image space are demagnified by M'
		Mpr = (nimg/nobj)*M;
		dir1M = dir1/Mpr;
		dir2M = dir2/Mpr;
		Array<double,3> D1M(lambda.size(),dir1.size(),dir2.size()),
                        D2M(lambda.size(),dir1.size(),dir2.size()),
                        D3M(lambda.size(),dir1.size(),dir2.size());
		D1M = D1/Mpr;
		D2M = D2/Mpr;
		if (Data.coll_half_space=="upper")
		{
			D3M = real(sqrt((Array<complex<double>,3>)(1-pow2(D1M)-pow2(D2M)))); //cast into complex, in case d1M^2+d2M^2>1
		}
		else if (Data.coll_half_space=="lower")
		{
			D3M = -real(sqrt((Array<complex<double>,3>)(1-pow2(D1M)-pow2(D2M)))); //cast into complex, in case d1M^2+d2M^2>1
		}

		Array<double,3> CosThetaImg(D3M.shape()),invCosThetaImg(D3M.shape());
		CosThetaImg = D3M;
		invCosThetaImg = 0; //initialize to zero
		for (int k=0; k<lambda_size; k++)
		{
            for (int d1=0;d1<dir1_size;d1++)
            {
                for (int d2=0;d2<dir2_size;d2++)
                {
                    if (CosThetaImg(k,d1,d2)!=0)
                    {
                        invCosThetaImg(k,d1,d2) = 1/CosThetaImg(k,d1,d2);
                    }//for non-propagating angles, leave invCosThetaImg at zero
                }
            }
		}

		Array<double,3> SinThetaImg(D3M.shape()),CosPhiImg(D3M.shape()),SinPhiImg(D3M.shape());
		SinThetaImg = real(sqrt((Array<complex<double>,3>)(1-pow2(CosThetaImg)))); //cast into complex, in case there is a floating-point problem
		CosPhiImg = 1;
		SinPhiImg = 0;
		for (int k=0; k<lambda_size; k++)
		{
            for (int d1=0;d1<dir1_size;d1++)
            {
                for (int d2=0;d2<dir2_size;d2++)
                {
                    if (SinThetaImg(k,d1,d2)!=0)
                    {//if sin(theta)=0, define phi to be zero, leave cos(phi) at 1 and sin(phi) at 0
                        CosPhiImg(k,d1,d2) = D1M(k,d1,d2)/SinThetaImg(k,d1,d2);
                        SinPhiImg(k,d1,d2) = D2M(k,d1,d2)/SinThetaImg(k,d1,d2);
                    }
                }
            }
		}

		// definition of E_s'
		// (CosThetaImg and CosThetaObj should have the same sign, so the sqrt
		// should work out fine)
		// but check anyway:
		if (product((CosThetaImg*invCosThetaObj)<0)!=0) //are there negative elements in (CosThetaImg*invCosThetaObj)?
		{
			if (rank==0)
			{
				cout << "Developer error: CosThetaImg and invCosThetaObj should have the same sign!" << endl;
			}//maybe implement C++ exception?
			MPI_exit(-1);
		}
		Array<complex<double>,3> E_theta_img(lambda_size,dir1_size,dir2_size),E_phi_img(lambda_size,dir1_size,dir2_size);
//		Array<complex<double>,3> prefactor3D(prefactor.data(),shape(1,dir1_size,dir2_size),neverDeleteData);
//		for (int k=0; k<lambda_size; k++)
//		{
//			E_theta_img(k,Range::all(),Range::all()) = prefactor3D*E_theta(k,Range::all(),Range::all());
//			E_phi_img(k,Range::all(),Range::all()) = prefactor3D*E_phi(k,Range::all(),Range::all());
//		}
        E_theta_img = -M*sqrt(nimg/nobj)*sqrt(CosThetaImg*invCosThetaObj)*E_theta;
        E_phi_img = -M*sqrt(nimg/nobj)*sqrt(CosThetaImg*invCosThetaObj)*E_phi;

		// determine the FFT lengths
		int N_fft_x = (int)ceil(dir1_size*Data.img_oversamplingrate_x);
		int N_fft_y = (int)ceil(dir2_size*Data.img_oversamplingrate_y);

		if (Data.choose_smallest_power_of_2)
		{
			N_fft_x = (int)pow(2,ceil(log(N_fft_x)/log(2.0)));
			N_fft_y = (int)pow(2,ceil(log(N_fft_y)/log(2.0)));
		}
// cout << N_fft_x << "," << N_fft_y << endl;

		//the following are allocated to the FFT dimensions, because the KissFFT library does not support arbitrary FFT lengths
		Array<complex<double>,3> E_s_x(lambda_size,N_fft_x,N_fft_y),E_s_y(lambda_size,N_fft_x,N_fft_y),E_s_z(lambda_size,N_fft_x,N_fft_y);
		//initialize to zero (very important because this amounts to zero-padding the array to the FFT dimensions N_fft_x,N_fft_y)
		E_s_x = 0;
		E_s_y = 0;
		E_s_z = 0;
		//define E_s_x, E_s_y, E_s_z
		E_s_x(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) = CosThetaImg*CosPhiImg*E_theta_img - SinPhiImg*E_phi_img;
		E_s_y(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) = CosThetaImg*SinPhiImg*E_theta_img + CosPhiImg*E_phi_img;
		E_s_z(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) = -SinThetaImg*E_theta_img;

		//definition of G'
        E_s_x(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) *= invCosThetaImg;
        E_s_y(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) *= invCosThetaImg;
        E_s_z(Range(0,lambda_size-1),Range(0,dir1_size-1),Range(0,dir2_size-1)) *= invCosThetaImg;

		if ((dir1M.size()<2)||(dir2M.size()<2))
		{
			if (rank==0)
			{
				cout << "Developer error: dir1M and dir2M should have lengths > 1" << endl;
			}//maybe implement C++ exception?
			MPI_exit(-1);
		}
		double dksxpr = kmax*(dir1M(1)-dir1M(0));
		double dksypr = kmax*(dir2M(1)-dir2M(0));
		double x_max = 2*M_PI/dksxpr;
		double y_max = 2*M_PI/dksypr;
		double d_x = x_max/N_fft_x;
		double d_y = y_max/N_fft_y;
		// x and y ranges after fftshift'ing the arrays
		Array<double,1> x_range(N_fft_x);
		for (int i=0; i<N_fft_x; i++)  x_range(i) = (-ceil((N_fft_x-1)/2.0)+i)*d_x;
		Array<double,1> y_range(N_fft_y);
		for (int j=0; j<N_fft_y; j++)  y_range(j) = (-ceil((N_fft_y-1)/2.0)+j)*d_y;

		ones1.resize(N_fft_x);
		ones2.resize(N_fft_y);
		ones1=1;ones2=1;
		Array<double,2> XR(N_fft_x,N_fft_y),YR(N_fft_x,N_fft_y);
		XR = x_range(tensor::i)*ones2(tensor::j);	//outer product
		YR = ones1(tensor::i)*y_range(tensor::j);	//outer product

		// field components in the image space
		Array<complex<double>,3> Ex(lambda_size,N_fft_x,N_fft_y),Ey(lambda_size,N_fft_x,N_fft_y),Ez(lambda_size,N_fft_x,N_fft_y);
        Array<complex<double>,2> fft_output(N_fft_x,N_fft_y);
        //perform inverse 2D FFT
		for (int k=0; k<lambda_size; k++)
		{
		    //x component
		    perform_2D_FFT(E_s_x(k,Range::all(),Range::all()),fft_output,N_fft_x,N_fft_y,0);
		    fftshift_2D(fft_output);
		    Ex(k,Range::all(),Range::all()) = (ii/2.0/M_PI)*(dksxpr*dksypr)/kk(k)*fft_output;
		    //FFT assumes that the first element of E_s_* is at (0,0) but it is actually at (kmax*dir1M(0),kmax*dir2M(0))
		    //The resulting phase shift is recovered as follows:
            Ex(k,Range::all(),Range::all()) *= exp(-ii*(kmax*dir1M(0)*XR+kmax*dir2M(0)*YR));
		    //y component
		    perform_2D_FFT(E_s_y(k,Range::all(),Range::all()),fft_output,N_fft_x,N_fft_y,0);
		    fftshift_2D(fft_output);
		    Ey(k,Range::all(),Range::all()) = (ii/2.0/M_PI)*(dksxpr*dksypr)/kk(k)*fft_output;
            Ey(k,Range::all(),Range::all()) *= exp(-ii*(kmax*dir1M(0)*XR+kmax*dir2M(0)*YR)); //see explanation above
		    //z component
		    perform_2D_FFT(E_s_z(k,Range::all(),Range::all()),fft_output,N_fft_x,N_fft_y,0);
		    fftshift_2D(fft_output);
		    Ez(k,Range::all(),Range::all()) = (ii/2.0/M_PI)*(dksxpr*dksypr)/kk(k)*fft_output;
            Ez(k,Range::all(),Range::all()) *= exp(-ii*(kmax*dir1M(0)*XR+kmax*dir2M(0)*YR)); //see explanation above
		}

//cout << Ex << endl;
//cout << Ey << endl;
//cout << Ez << endl;

//// cout << array << endl;
//// cout << PW_x_phasorarray << endl;
//// cout << PW_y_phasorarray << endl;

        /** plane-wave contributions to the image **/

        double pw_theta,pw_phi;
        double dir1_pw,dir2_pw,dir3_pw,dir1M_pw,dir2M_pw,dir3M_pw;
        complex<double> PW_z_phasor;
        double radiometric_factor;
        bool pw_theta_within_NA;
        Array<complex<double>,3> E_pw_x(Ex.shape()),E_pw_y(Ey.shape()),E_pw_z(Ez.shape());

        for (int PW=0; PW<numofPWs; PW++)
        {
            pw_theta = PW_theta(PW);
            pw_phi = PW_phi(PW);

            /*************************/
            /**** COLLECTION PART ****/
            /*************************/

            //determine if the pw direction within the NA
            if (Data.coll_half_space=="upper")
            {
                pw_theta_within_NA = ((cos(pw_theta)>0)&&(abs(sin(pw_theta))<SinThetaMax));
            }
            else if (Data.coll_half_space=="lower")
            {
                pw_theta_within_NA = ((cos(pw_theta)<0)&&(abs(sin(pw_theta))<SinThetaMax));
            }

            /*************************/
            /**** REFOCUSING PART ****/
            /*************************/

            if (pw_theta_within_NA)
            {
                // original plane-wave direction cosines (without de-magnification)
                dir1_pw = sin(pw_theta)*cos(pw_phi);  // x-direction cosine
                dir2_pw = sin(pw_theta)*sin(pw_phi);  // y-direction cosine
                // z-direction cosine (discard the imaginary sq.roots)
                if (cos(pw_theta)>0)   // traveling into the +z half plane
                {
                    dir3_pw = real(sqrt((complex<double>)(1-pow2(dir1_pw)-pow2(dir2_pw))));
                }
                else if (cos(pw_theta)<0)   // traveling into the -z half plane
                {
                    dir3_pw = -real(sqrt((complex<double>)(1-pow2(dir1_pw)-pow2(dir2_pw))));
                }

                // define the de-magnified plane-wave direction cosines
                dir1M_pw = dir1_pw/Mpr;
                dir2M_pw = dir2_pw/Mpr;
                // z-direction cosine (discard the imaginary sq.roots)
                if (cos(pw_theta)>0)   // traveling into the +z half plane
                {
                    dir3M_pw = real(sqrt((complex<double>)(1-pow2(dir1M_pw)-pow2(dir2M_pw))));
                }
                else if (cos(pw_theta)<0)   // traveling into the -z half plane
                {
                    dir3M_pw = -real(sqrt((complex<double>)(1-pow2(dir1M_pw)-pow2(dir2M_pw))));
                }

                //calculate the contribution of the plane waves to the image field
                for (int k=0; k<lambda_size; k++)
                {
                    // To find the z component of the PW, again use the orthogonality of (PW_x,PW_y,PW_z) and (kx,ky,kz)
                    PW_z_phasor = -(1.0/dir3M_pw)*(dir1M_pw*PW_x_phasorarray(PW,k)+dir2M_pw*PW_y_phasorarray(PW,k));

					//No phase shift is needed, because the phasors are already referenced to the TRANSFORMER origin, which is also the origin of the image
//                    // apply the extra phases due to the shift in origin
//                    // Note that the original direction cosines (without de-magnification) are
//                    // used, because they determine the phases acquired before the light reaches
//                    // the lens.
//                    // NOTE: ImgOriginX etc. is not relative to the grid origin, so subtraction of OriginX etc. is necessary.
//                    origin_phase = exp(-ii*kk(k)*(dir1_pw*(Data.ImgOriginX-OriginX)*dx+dir2_pw*(Data.ImgOriginY-OriginY)*dx+dir3_pw*(Data.ImgOriginZ-OriginZ)*dx));

                    // divide the plane-wave phasors by M to keep the total energy in the plane wave (area*intensity) constant. Since
                    // the area gets magnified by M^2, the phasor amplitude should be reduced by M.
                    // Note that dir3_pw is cos(th) in object space, dir3M_pw is cos(th') in image space. They should have the same sign,
                    // so the sqrt shouldn't be a problem.
                    // But let's check anyway:
                    if (dir3_pw*dir3M_pw<0)
                    {
                        if (rank==0)
                        {
                            cout << "Developer error (Cimg class): There is a sign error." << endl;
                        }
                        MPI_exit(-1);
                    }
                    radiometric_factor = sqrt((nimg*dir3_pw)/(nobj*dir3M_pw))/M;

                    E_pw_x(k,Range::all(),Range::all()) += PW_x_phasorarray(PW,k)*exp(-ii*kk(k)*(dir1M_pw*XR+dir2M_pw*YR))*radiometric_factor;
                    E_pw_y(k,Range::all(),Range::all()) += PW_y_phasorarray(PW,k)*exp(-ii*kk(k)*(dir1M_pw*XR+dir2M_pw*YR))*radiometric_factor;
                    E_pw_z(k,Range::all(),Range::all()) += PW_z_phasor*exp(-ii*kk(k)*(dir1M_pw*XR+dir2M_pw*YR))*radiometric_factor;
                }
            }
        }

#else //#ifndef HDF5_DISABLE
		throw AngoraDeveloperException("Custom file format is now obsolete for far-field output files.");
#endif //#ifndef HDF5_DISABLE

		/****************************************************************/
		/****************************************************************/
		/****************************************************************/
		/*******       WRITE THE DATA TO IMAGE FILE       ***************/
		/****************************************************************/
		/****************************************************************/
		/****************************************************************/

#ifndef HDF5_DISABLE
		//create the image file
		H5File image_file(ImageFileName.c_str(), H5F_ACC_TRUNC);
		/** write the fdtd code version **/
		// Default property list
		DSetCreatPropList plist;
		//dataspace for single double value
		hsize_t dims=3;
		DataSpace three_val_dspace(1,&dims);
		dataset=image_file.createDataSet("angora_version", PredType::NATIVE_INT, three_val_dspace, plist);
		//write the package version (in int)
		int version_array[3] = {angora_version_major,angora_version_minor,angora_version_revision};
		dataset.write(&version_array,PredType::NATIVE_INT);

		//dataspace for the 3D image arrays
		hsize_t imagedims[] = {lambda_size, N_fft_x, N_fft_y};
		DataSpace dspace( 3, imagedims );

		//dataspace for single double value
		dims=1;
		DataSpace double_dspace(1,&dims);

		/** magnification **/
		dataset=image_file.createDataSet("magnification", PredType::NATIVE_DOUBLE, double_dspace, plist);
		dataset.write(&M,PredType::NATIVE_DOUBLE);
		/** aperture half angle (in rad) **/
		dataset=image_file.createDataSet("ap_half_angle", PredType::NATIVE_DOUBLE, double_dspace, plist);
		dataset.write(&theta_max,PredType::NATIVE_DOUBLE);
		/** image-space refractive index **/
		dataset=image_file.createDataSet("n_img", PredType::NATIVE_DOUBLE, double_dspace, plist);
		dataset.write(&nimg,PredType::NATIVE_DOUBLE);
		/** object-space refractive index **/
		dataset=image_file.createDataSet("n_obj", PredType::NATIVE_DOUBLE, double_dspace, plist);
		dataset.write(&nobj,PredType::NATIVE_DOUBLE);
		/** x and y ranges **/
		hsize_t x_range_length = x_range.size();
		DataSpace xspace(1,&x_range_length);
		dataset=image_file.createDataSet("x_range", PredType::NATIVE_DOUBLE, xspace, plist);
		dataset.write(x_range.data(),PredType::NATIVE_DOUBLE);
		hsize_t y_range_length = y_range.size();
		DataSpace yspace(1,&y_range_length);
		dataset=image_file.createDataSet("y_range", PredType::NATIVE_DOUBLE, yspace, plist);
		dataset.write(y_range.data(),PredType::NATIVE_DOUBLE);
		/** wavelength range **/
		hsize_t wv_range_length = lambda.size();
		DataSpace wvspace(1,&wv_range_length);
		dataset=image_file.createDataSet("wv_range", PredType::NATIVE_DOUBLE, wvspace, plist);
		dataset.write(lambda.data(),PredType::NATIVE_DOUBLE);
		//also write the wavenumber 2*pi/lambda
		hsize_t k_range_length = kk.size();
		DataSpace kspace(1,&k_range_length);
		dataset=image_file.createDataSet("k_range", PredType::NATIVE_DOUBLE, kspace, plist);
		dataset.write(kk.data(),PredType::NATIVE_DOUBLE);

		//create temporary arrays to hold the 3D and 2D data [less elegant, but more portable than creating strided dataspaces]
		temp.resize(lambda_size, N_fft_x, N_fft_y);
		temp2D.resize(N_fft_x, N_fft_y);

		/** x-component of the total E-field **/
		if (Data.output_data_item_exists("E_x_tot"))
		{
			dataset=image_file.createDataSet("E_x_tot_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ex+E_pw_x);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_x_tot_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ex+E_pw_x);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** y-component of the total E-field **/
		if (Data.output_data_item_exists("E_y_tot"))
		{
			dataset=image_file.createDataSet("E_y_tot_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ey+E_pw_y);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_y_tot_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ey+E_pw_y);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** z-component of the total E-field **/
		if (Data.output_data_item_exists("E_z_tot"))
		{
			dataset=image_file.createDataSet("E_z_tot_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ez+E_pw_z);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_z_tot_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ez+E_pw_z);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}

		/** x-component of the scattered E-field **/
		if (Data.output_data_item_exists("E_x_sca"))
		{
			dataset=image_file.createDataSet("E_x_sca_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ex);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_x_sca_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ex);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** y-component of the scattered E-field **/
		if (Data.output_data_item_exists("E_y_sca"))
		{
			dataset=image_file.createDataSet("E_y_sca_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ey);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_y_sca_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ey);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** z-component of the scattered E-field **/
		if (Data.output_data_item_exists("E_z_sca"))
		{
			dataset=image_file.createDataSet("E_z_sca_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(Ez);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_z_sca_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(Ez);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}

		/** x-component of the incident E-field **/
		if (Data.output_data_item_exists("E_x_inc"))
		{
			dataset=image_file.createDataSet("E_x_inc_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(E_pw_x);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_x_inc_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(E_pw_x);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** y-component of the incident E-field **/
		if (Data.output_data_item_exists("E_y_inc"))
		{
			dataset=image_file.createDataSet("E_y_inc_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(E_pw_y);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_y_inc_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(E_pw_y);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
		/** z-component of the incident E-field **/
		if (Data.output_data_item_exists("E_z_inc"))
		{
			dataset=image_file.createDataSet("E_z_inc_r", PredType::NATIVE_DOUBLE, dspace, plist);
			//write the real part
			temp = real(E_pw_z);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
			//write the imaginary part
			dataset=image_file.createDataSet("E_z_inc_i", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = imag(E_pw_z);
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}

		/** total intensity **/
		if (Data.output_data_item_exists("intensity_tot"))
		{
			dataset=image_file.createDataSet("intensity_tot", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = pow2(abs(Ex+E_pw_x))+pow2(abs(Ey+E_pw_y))+pow2(abs(Ez+E_pw_z));
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
//			/** total bright-field intensity **/
//			dataset=image_file.createDataSet("bfintensity_tot", PredType::NATIVE_DOUBLE, dspace2D, plist);
//			temp2D = sum(temp(tensor::k,tensor::i,tensor::j),tensor::k); //blitz can only reduce over the last dimension, so make lambda the last one
//			dataset.write(temp2D.data(),PredType::NATIVE_DOUBLE);

		/** scattered intensity **/
		if (Data.output_data_item_exists("intensity_sca"))
		{
			dataset=image_file.createDataSet("intensity_sca", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = pow2(abs(Ex))+pow2(abs(Ey))+pow2(abs(Ez));
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}

		/** incident intensity **/
		if (Data.output_data_item_exists("intensity_inc"))
		{
			dataset=image_file.createDataSet("intensity_inc", PredType::NATIVE_DOUBLE, dspace, plist);
			temp = pow2(abs(E_pw_x))+pow2(abs(E_pw_y))+pow2(abs(E_pw_z));
			dataset.write(temp.data(),PredType::NATIVE_DOUBLE);
		}
//			/** incident bright-field intensity **/
//			dataset=image_file.createDataSet("bfintensity_inc", PredType::NATIVE_DOUBLE, dspace2D, plist);
//			temp2D = sum(temp(tensor::k,tensor::i,tensor::j),tensor::k); //blitz can only reduce over the last dimension, so make lambda the last one
//			dataset.write(temp2D.data(),PredType::NATIVE_DOUBLE);

		temp.resize(0,0,0);  //deallocate memory for the temporary array
		temp2D.resize(0,0);  //deallocate memory for the temporary array
		image_file.close();

#else //#ifndef HDF5_DISABLE
		if (rank==0)
		{
			cout << "Error: Custom file format not implemented for optical images." << endl;
		}
		MPI_exit(-1);
#endif //#ifndef HDF5_DISABLE
	}
}

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

#ifndef RANDOM_H
#define RANDOM_H

//template declaration for the function that places a random material region

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

//for the gamma function in the cephes math library
#include "specfun/gamma/gamma.h"

//for reading a material region from a file
#include "matfile.h"

//fft routines
#include "fft/fft.h"

//for the blitz++ Gaussian-distributed random number generator
#include <random/normal.h>

#include <fstream>

//standard c time library for seeding the random-number generator
#include <ctime>

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

extern double dx;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif

extern void MPI_exit(const int& exitcode);



template <typename RandMatType>
void PlaceWMCorrRandomBlock
(//Places a rectangular block of inhomogeneous random material with Whittle-Matern correlation (permittivity, permeability or electric_conductivity)
			const string& constitutive_param_type,
			const double& mean, const double& std_dev, const double& corr_len, const double& m, //statistical parameters of the permittivity distr.
			const double& back_coord, const double& front_coord, //coordinates of the back and front surfaces of the block
			const double& left_coord, const double& right_coord, //coordinates of the left and right surfaces of the block
			const double& lower_coord, const double& upper_coord, //coordinates of the lower and upper surfaces of the block
			const int& random_seed) //random integer seed for the random-number generator
{
	//The following is copied from PlaceMaterialRegionFromFile in order to avoid calculating the material region unnecessarily
	if ((constitutive_param_type!="rel_permittivity")&&(constitutive_param_type!="rel_permeability")&&(constitutive_param_type!="electric_conductivity"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,constitutive_param_type,
			"(valid arguments are \"rel_permittivity\", \"rel_permeability\", or \"electric_conductivity\")");
	}

	if (random_seed<0)
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<int>(func_name,random_seed,
			"(random seed should be a nonnegative integer)");
	}

	int BlockBack = (int)ceil(back_coord)+1;
	int BlockFront = (int)floor(front_coord);
	int BlockLeft = (int)ceil(left_coord)+1;
	int BlockRight = (int)floor(right_coord);
	int BlockLower = (int)ceil(lower_coord)+1;
	int BlockUpper = (int)floor(upper_coord);

	int A = BlockFront-BlockBack+1;
	int B = BlockRight-BlockLeft+1;
	int C = BlockUpper-BlockLower+1;

///** FOR DEBUGGING ONLY **/
//A=300;
//B=300;
//C=300;
///** FOR DEBUGGING ONLY **/

	//this will only be modified by node 0, and broadcasted to the other nodes by the same node
//	char randfile_template[] = "/tmp/fdtd-XXXXXX"; /** GNU-specific **/ //tmp directory location
/** TODO: FInd a better solution for this!! **/
	char randfile_template[] = "./fdtd-XXXXXX"; /** GNU-specific **/ //UNTIL A BETTER SOLUTION IS FOUND!!!
/** TODO: FInd a better solution for the above!! **/

/*******************************************************************************************************************/
/** Write the random medium into temporary file. **/
	if (rank==0)
	{
		//first open the stream and write the preamble
		//C temporary file stream:
		FILE * randfilestream;
		int randfiledes = mkstemp(randfile_template); /** GNU-specific **/ //creates temporary file, returns low-level file descriptor (note that "file descriptors" are GNU-specific)
		if (randfiledes<0)
		{
			if (rank==0)
				cout << "Error: Cannot create temporary file for random medium" << endl;
			MPI_exit(-1);
		}
		randfilestream = fdopen(randfiledes,"w+"); /** GNU-specific **/ //convert file descriptor to higher-level stream, open stream for write/read (recreate if exists)
		//first, write the x,y,z extents (in cells)
		fwrite(&A, sizeof(A), 1, randfilestream);
		fwrite(&B, sizeof(B), 1, randfilestream);
		fwrite(&C, sizeof(C), 1, randfilestream);

/* Generate a three-dimensional (3-D) zero-mean real correlated array x[a,b,c] of dimensions AxBxC.
   std_dev is the variance of the random array : E{x^2}
   lc_cell is the correlation length in grid cells : nc = lc/dx, where dx is the spatial increment, and lc is the correlation length
                                                   in the continuous-domain correlation function.
   The correlation function (called the Whittle-Matern family of correlations) is given by
    Bn[a,b,c]=Bn[r]=2^(5/2-m)(std_dev^2)*(r/nc)^(m-3/2)*BesselK(m-3/2,r/nc)/
                 Gamma(m-3/2)
   Here, r = (a^2+b^2+c^2)^1/2. The correlation function is circularly symmetric in the (a,b,c) space.
   For m=2, this becomes the exponential correlation function: Bn[r,m=2]=std_dev^2*exp(-r/nc)
   To avoid aliasing in the spatial-frequency domain, nc>>(1/pi) must be satisfied, so that 1/(1+(K*nc).^2)^m
   decays sufficiently at K=pi.
   To avoid aliasing in the spatial domain, (A/nc)>>1, (B/nc)>>1, and (C/nc)>>1 must be satisfied separately, so that
   the correlation decays sufficiently at a=A and b=B and c=C.
*/
		double lc_cell = corr_len/dx;
		double std_dev_norm = std_dev/mean;

		//initialize the random number generator
		ranlib::Normal<RandMatType> NormalRV(0.0,1.0); //Gaussian r.v. with zero mean, unit standard deviation
		NormalRV.seed((unsigned int)random_seed); //seed using the given (positive) integer

		/**************************************************************************/
		/** Memory-inefficient but faster method that uses an auxiliary 3D array **/
//		Array<complex<RandMatType>,1> input(3);
//		Array<RandMatType,1> output(4);
//		input(0) = 10.0; input(1) = -2.0+2.0*ii; input(2) = -2.0;
//		output = 0;
//		int N_FFT = 4;
//		cout << input << endl;
//		perform_1D_FFT_sym_i(input,output,N_FFT,1);
//		output = output/double(N_FFT);
//		cout << (output) << endl;

//		Array<RandMatType,1> input(4);
//		Array<complex<RandMatType>,1> output(3);
//		input(0) = 1.0; input(1) = 2.0; input(2) = 3.0; input(3) = 4.0;
//		output = 0;
//		int N_FFT = 4;
//		cout << input << endl;
//		perform_1D_FFT_sym_r(input,output,N_FFT,0);
//		cout << output << endl;

//		int N1 = 4;
//		int N2 = 6;
//		Array<RandMatType,2> input(N1,N2);
//		Array<complex<RandMatType>,2> output(N1,N2/2+1);
//		input = 1,2,3,4,5,6,7,8,9,11,10,12,13,14,15,16,9,11,10,12,13,14,15,16;
//		output = 0;
//		cout << input << endl;
//		perform_2D_FFT_sym_r(input,output,N1,N2,0);
//		cout << output << endl;
//exit(-1);

//		Array<RandMatType,3> input(3,3,4);
//		Array<complex<RandMatType>,3> output(3,3,3);
//		input = 2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36;
//		output = 0;
//		int N1 = 3;
//		int N2 = 3;
//		int N3 = 4;
//		cout << input << endl;
//		perform_3D_FFT_sym_r(input,output,N1,N2,N3,0);
//		cout << output << endl;
//exit(-1);

		Array<RandMatType,3> pow_spec_array;
		Array<complex<RandMatType>,3> rand_output;

		try{//try to allocate memory for the 3D spectrum arrays

		//Create the array backward, since PlaceMaterialRegionFromFile reads the x-dimension first (the files used to be always created by matlab)
		int even_A = 2*((A+1)/2); //the last dim. has to be even for the KissFFT conjugate-symmetry functionality to work

		pow_spec_array.resize(C,B,even_A);  //real power-spectrum array
		rand_output.resize(C,B,even_A/2+1); //conj. symmetric output of half the size

		{double k_ijk;
		int i,j,k;
		for (k=0; k<C; k++)
		{
			for (j=0; j<B; j++)
			{
				for (i=0; i<even_A; i++)
				{
					k_ijk = sqrt(pow2(2*M_PI/C*(-floor(C/2.0)+k))+
								 pow2(2*M_PI/B*(-floor(B/2.0)+j))+
								 pow2(2*M_PI/even_A*(-floor(even_A/2.0)+i)));
					pow_spec_array(k,j,i) = NormalRV.random()*
								sqrt(C*B*even_A*  //each frequency is independent with power (C*B*even_A*powerspectrum(k_ijk))
								pow2(std_dev_norm)*(8*pow(M_PI,1.5)*pow((RandMatType)lc_cell,(RandMatType)3)*gamma(m)/(abs(gamma(m-1.5))*pow((RandMatType)(1+pow2(k_ijk)*pow2(lc_cell)),(RandMatType)m)))); //power spectral density of the Whittle-Matern-correlated medium
						//extra casts are needed for some versions of gcc
				}
			}
		}
		}
		//bring the zero-frequency component back to the beginning
		ifftshift_3D(pow_spec_array);
		//take 3D FFT of array, get conjugate-symmetric complex array of half the size
///** FOR DEBUGGING ONLY **/
//pow_spec_array = 0.807013,0.0854601,0.0282915,0.0109957,-0.00160347,0.599055,-0.158318,
//0.0307762,0.0238219,0.0179352,-0.0255191,0.00304019,-0.0103204,0.0608063,
//0.0608514,-0.0196143,0.00681328,0.078054,-0.205359,0.0716104,0.0592157,
//-0.0636565,-0.012153,0.0678146,-0.0703019,-0.0212802,-0.0307865,-0.0162769,
//-0.0417819,-0.0151659,-0.00488557,0.0521722,-0.0371249,-0.0240599,0.0351162,
//0.0191775,0.10874,0.0198407,-0.0487446,-0.0397812,-0.0457832,0.0263564,
//0.0858933,0.0360089,-0.0414271,-0.00413052,-0.0323868,0.119837,0.0874144,
//0.0205548,0.0304685,-0.00340674,0.0688936,-0.0691354;
///** FOR DEBUGGING ONLY **/
//cout << pow_spec_array << endl;
		perform_3D_FFT_sym_r(pow_spec_array, rand_output, C, B, even_A, 0); //does NOT work with inverse FFT, so take regular FFT
		rand_output = rand_output/(RandMatType)(C*B*even_A); //inverse FFT needs to be scaled
		//rand_output is of dimensions C x B x (even_A/2+1)
//cout << rand_output(0,1,2) << endl;
//cout << rand_output << endl;
//cout << "Finished computing" << endl;
//exit(-1);

		//scale appropriately to get the correct mean and std.dev.
		rand_output = (RandMatType)mean*((RandMatType)1.0+rand_output);

		//now, create real array from conjugate-symmetric complex array, write into temporary file
		{int i,j,k;
		Array<RandMatType,1> real_rand_column(A); //the actual dimension of the random array will be A
		Array<complex<RandMatType>,1> this_column(even_A/2+1),conjugate_column(even_A/2+1); //temporary complex columns from which the real data is extracted (these two are symmetrically positioned w.r.t. the center)
		for (k=0; k<C; k++)
		{
			for (j=0; j<B; j++)
			{
				this_column = rand_output(k,j,Range::all());
				conjugate_column = rand_output((C-k)%C,(B-j)%B,Range::all()); //symmetric position
//cout << this_column << endl;
				real_rand_column(Range(0,even_A/2)) = real(this_column(Range(0,even_A/2)))+imag(this_column(Range(0,even_A/2)));
				if (A>=3)
				real_rand_column(Range(even_A/2+1,A-1)) = real(conjugate_column(Range(even_A/2-1,even_A-A+1,-1)))-imag(conjugate_column(Range(even_A/2-1,even_A-A+1,-1)));
				fwrite(real_rand_column.data(),sizeof(RandMatType),A,randfilestream);
			}
		}
		}
//
//cout << randfile_template << endl;
//fclose(randfilestream);
//double aaa;cin >> aaa;
//cout << "Finished computing" << endl;
//exit(-1);
		}//end of try block
		/** Memory-inefficient but faster method that uses an auxiliary 3D array **/
		/**************************************************************************/

		catch (bad_alloc&)
		{/** Ooops, not enough memory. Why not use the disk then? Sure... **/
			/**********************************************************************************/
			/** MUCH slower method that uses a temporary auxiliary file and not much memory **/
			//secondary temporary file to hold the intermediate FFT information
			//C temporary file stream:
			FILE * tmpfilestream;
			char tmpfile_template[] = "/tmp/fdtd-XXXXXX"; /** GNU-specific **/ //tmp directory location
			int tmpfiledes = mkstemp(tmpfile_template); /** GNU-specific **/ //creates temporary file, returns low-level file descriptor (note that "file descriptors" are GNU-specific)
			if (tmpfiledes<0)
			{
				if (rank==0)
					cout << "Error: Cannot create temporary file for random medium" << endl;
				MPI_exit(-1);
			}
			tmpfilestream = fdopen(tmpfiledes,"w+"); /** GNU-specific **/ //convert file descriptor to higher-level stream, open stream for write/read (recreate if exists)

			{double k_ijk;
			 int i,j,k;
			 Array<complex<RandMatType>,1> last_dim(C);

			//PlaceMaterialRegionFromFile reads the x-dimension first (since the files have always been created by matlab)
			//So, create the z dimension first, place in strided fashion into file.
			for (i=0; i<A; i++)
			{
				for (j=0; j<B; j++)
				{
					for (k=0; k<C; k++)
					{
						k_ijk = sqrt(pow2(2*M_PI/A*(-floor(A/2.0)+i))+
									 pow2(2*M_PI/B*(-floor(B/2.0)+j))+
									 pow2(2*M_PI/C*(-floor(C/2.0)+k)));
						last_dim(k) = NormalRV.random()*
									sqrt(A*B*C*  //each frequency is independent with power (A*B*C*powerspectrum(k_ijk))
									pow2(std_dev_norm)*(8*pow(M_PI,1.5)*pow((RandMatType)lc_cell,(RandMatType)3)*gamma(m)/(abs(gamma(m-1.5))*pow((RandMatType)(1+pow2(k_ijk)*pow2(lc_cell)),(RandMatType)m)))); //power spectral density
							//extra casts are needed for some versions of gcc
					}
					//bring the zero-frequency component back to the beginning
					ifftshift_1D(last_dim);
					//take 1D FFT of last_dim, place in random_array
					perform_1D_FFT(last_dim, last_dim, C, 1); //take inverse FFT
					last_dim = last_dim/(complex<RandMatType>)C; //KissFFT inverse FFT does not normalize by the length
					//Place the Fourier transform of the 1D array in the temp file at the appropriate position
					/** (write into tmpfilestream...) **/
					if (!fseek(tmpfilestream, (i+A*j)*sizeof(complex<RandMatType>),SEEK_SET))
					{
						/** throw exception **/
					}
					for (k=0; k<C; k++)
					{
						if (fwrite(last_dim.data()+k,sizeof(complex<RandMatType>),1,tmpfilestream)!=1)
						{
							/** throw exception **/
						}
						if (!fseek(tmpfilestream, (A*B-1)*sizeof(complex<RandMatType>),SEEK_CUR))
						{
							/** throw exception **/
						}
					}
				}
			}
			//close temporary file stream
			fclose(tmpfilestream);
			//reopen temporary file stream for reading (use modified char* template directly)
			tmpfilestream = fopen(tmpfile_template,"r");

			//Transform the x and y dimensions next
			//We want column-order placement in the file, so higher rank indices are written/read first
			Array<complex<RandMatType>,2> first_two_dims(B,A);
			Array<RandMatType,2> first_two_dims_real(B,A);
			for (k=0; k<C; k++)
			{
				//read the first two dims from the secondary tmp file
				/** (read from tmpfilestream ...) **/
				if (fread((char*)first_two_dims.data(),sizeof(complex<RandMatType>),A*B,tmpfilestream) < A*B)
				{
					/** throw exception **/
				}
				//bring the zero-frequency component back to the corner
				ifftshift_2D(first_two_dims);
				perform_2D_FFT(first_two_dims, first_two_dims, B, A, 1); //take inverse FFT
				first_two_dims = first_two_dims/(complex<RandMatType>)(A*B); //KissFFT inverse FFT does not normalize by the length
				//add real and imaginary parts to get real array, then rescale to the correct mean and std. dev.
				first_two_dims_real = mean*(1.0+(real(first_two_dims)+imag(first_two_dims)));
				//write the 2D real array into randfilestream
				/** (write into randfilestream ...) **/
				if (fwrite(first_two_dims_real.data(),sizeof(RandMatType),A*B,randfilestream)!=A*B)
				{
					/** throw exception **/
				}
			}
			}

			//close temporary file stream (this time for good!)
			fclose(tmpfilestream);
			//delete the temporary file too (use modified char* template directly)
			if (unlink(tmpfile_template)<0)
			{
				if (rank==0)
					cout << "Warning: Cannot delete temporary file for random medium" << endl;
			}
		/** MUCH slower method that uses a temporary auxiliary file and not much memory **/
		/*********************************************************************************/
	} //end of catch block


		//finally, close temporary-file stream for the random-material file
		fclose(randfilestream);
	}
/** Random medium written into temporary file (by node 0).
/*******************************************************************************************************************/

#ifndef MPI_DISABLE
	//Get the temporary file name from node 0
	MPI_Bcast(randfile_template,strlen(randfile_template),MPI_CHAR,0,MPI_CartSubComm);
#endif
	//name of the temporary file (for later use by all nodes when calling PlaceMaterialRegionFromFile)
	string randfilename = randfile_template; //copy the filename into string (needed for PlaceMaterialRegionFromFile)

#ifndef MPI_DISABLE
	MPI_Barrier(MPI_CartSubComm);	//other nodes should wait until the master node finishes writing
#endif

	//all nodes should read random material region from file
	PlaceMaterialRegionFromFile<RandMatType>(randfilename,BlockBack-1,BlockLeft-1,BlockLower-1,"BLL",constitutive_param_type); //the coordinates of the "anchor" is 1 less than the cell indices of the back-left-lower cell

#ifndef MPI_DISABLE
	MPI_Barrier(MPI_CartSubComm);	//the file won't be closed (or deleted) until all nodes read from it
#endif

//cout << randfilename << endl;
//double a;
//cin >> a;

	//finally, delete the temporary random-material file
	if (rank==0)
	{
		if (unlink(randfilename.c_str())<0)
		{
			if (rank==0)
			{
				cout << "Warning: Cannot delete temporary file for random medium" << endl;
			}
		}
	}
}

template <typename RandMatType>
void PlaceWMCorrRandomBlock
(const string& constitutive_param_type, const double& mean, const double& std_dev, const double& corr_len, const double& m, const double& back_coord, const double& front_coord, const double& left_coord, const double& right_coord, const double& lower_coord, const double& upper_coord)
{//the version of PlaceWMCorrRandomBlock that creates its own random seed from the system time
	int random_seed = (unsigned int)time(NULL); //seed using the current time
	PlaceWMCorrRandomBlock<RandMatType>(constitutive_param_type,mean,std_dev,corr_len,m,back_coord,front_coord,left_coord,right_coord,lower_coord,upper_coord,random_seed);
}

#endif // RANDOM_H

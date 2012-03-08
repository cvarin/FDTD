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

//Definitions of the FFT wrappers (float type)

#include "headers.h"

//#include "fft.h"

// #include <gsl/gsl_fft_complex.h>
// #include <fftw3.h>

//KissFFT library header
#include "kiss_F_fft/kiss_F_fftnd.h"
#include "kiss_F_fft/kiss_F_fftndr.h"
// Although KissFFT does not *natively* support arbitrary FFT sizes, this can be partially alleviated by resizing the input to the desired FFT size and applying zero padding.


/** The internals of these functions are the same as in fft_d.cpp. **/
/** (otherwise the KissFFT library has to be c++'d with templates) **/
void perform_1D_FFT(const Array<complex<float>,1>& input, Array<complex<float>,1>& output, const int& N_FFT, const int& inverse_fft)
{//computes the 1D Fourier transform
	kiss_F_fft_cfg mycfg=kiss_F_fft_alloc(N_FFT,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0  NOTE:inverse FFT is NOT scaled by 1/(N_FFT)
	//Here, we assume that kiss_F_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_F_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fft(mycfg,input_kissfft,output_kissfft);
	free(mycfg);
}

void perform_1D_FFT_sym_i(const Array<complex<float>,1>& input, Array<float,1>& output, const int& N_FFT, const int& inverse_fft)
{//(Works only for even N_FFT)
 //computes the 1D Fourier transform for a conjugate-symmetric input array of length N_FFT/2+1
 //output will be of length N_FFT
	if (2*(N_FFT/2)!=N_FFT)
	{
		/** throw exception **/
		exit(-1); //remove later
	}
	kiss_F_fftr_cfg mycfg=kiss_F_fftr_alloc(N_FFT,inverse_fft,NULL,NULL);
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_scalar* output_kissfft = reinterpret_cast<kiss_F_fft_scalar*>(output.data());
	kiss_F_fftri(mycfg,input_kissfft,output_kissfft);
	/*
	 input has  nfft/2+1 complex points
	 output has nfft scalar points
	*/
	free(mycfg);
}

void perform_1D_FFT_sym_r(const Array<float,1>& input, Array<complex<float>,1>& output, const int& N_FFT, const int& inverse_fft)
{//(Works only for even N_FFT)
 //computes the 1D Fourier transform for a real input array of length N_FFT, returns one half of the conjugate-symmetric output array
 //output will be of length N_FFT/2+1
	if (2*(N_FFT/2)!=N_FFT)
	{
		/** throw exception **/
		exit(-1); //remove later
	}
	kiss_F_fftr_cfg mycfg=kiss_F_fftr_alloc(N_FFT,inverse_fft,NULL,NULL);
	const kiss_F_fft_scalar* input_kissfft = reinterpret_cast<const kiss_F_fft_scalar*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fftr(mycfg,input_kissfft,output_kissfft);
	/*
	 input has nfft scalar points
	 output has nfft/2+1 complex points
	*/
	free(mycfg);
}

void perform_2D_FFT(const Array<complex<float>,2>& input, Array<complex<float>,2>& output, const int& N1, const int& N2, const int& inverse_fft)
{//computes the 2D Fourier transform
	const int ndims = 2;
	const int dims[ndims] = {N1,N2};
	kiss_F_fftnd_cfg mycfg=kiss_F_fftnd_alloc(dims,ndims,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0   NOTE:inverse FFT is NOT scaled by 1/(N1*N2)
	//Here, we assume that kiss_F_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_F_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fftnd(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_2D_FFT_sym_r(const Array<float,2>& input, Array<complex<float>,2>& output, const int& N1, const int& N2, const int& inverse_fft)
{//(Works only for even N2)
 //computes the 2D Fourier transform for a conjugate-symmetric input array of dimensions N1 x (N2/2+1)
 //output will be of dimensions N1 x N2
	const int ndims = 2;
	const int dims[ndims] = {N1,N2};
	kiss_F_fftndr_cfg mycfg=kiss_F_fftndr_alloc(dims,ndims,inverse_fft,NULL,NULL);
	const kiss_F_fft_scalar* input_kissfft = reinterpret_cast<const kiss_F_fft_scalar*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fftndr(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_2D_FFT_sym_i(const Array<complex<float>,2>& input, Array<float,2>& output, const int& N1, const int& N2, const int& inverse_fft)
{//(Works only for even N2)
 //computes the 2D Fourier transform for a real input array of dimensions N1 x N2, returns one half of the conjugate-symmetric output array
 //output will be of dimensions N1 x (N2/2+1)
	const int ndims = 2;
	const int dims[ndims] = {N1,N2};
	kiss_F_fftndr_cfg mycfg=kiss_F_fftndr_alloc(dims,ndims,inverse_fft,NULL,NULL);
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_scalar* output_kissfft = reinterpret_cast<kiss_F_fft_scalar*>(output.data());
	kiss_F_fftndri(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_3D_FFT(const Array<complex<float>,3>& input, Array<complex<float>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft)
{//computes the 3D Fourier transform
	const int ndims = 3;
	const int dims[ndims] = {N1,N2,N3};
	kiss_F_fftnd_cfg mycfg=kiss_F_fftnd_alloc(dims,ndims,inverse_fft,NULL,NULL); //computes inverse FFT if inverse_fft is not 0   NOTE:inverse FFT is NOT scaled by 1/(N1*N2*N3)
	//Here, we assume that kiss_F_fft_cpx is bit-compatible to the C++ complex<T> datatype
	//kiss_F_fft_cpx is just a struct with two elements of type T, so it is sufficiently simple for this assumption
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fftnd(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_3D_FFT_sym_r(const Array<float,3>& input, Array<complex<float>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft)
{//(Works only for even N3)
 //computes the 3D Fourier transform for a conjugate-symmetric input array of dimensions N1 x N2 x (N3/2+1)
 //output will be of dimensions N1 x N2 x N3
	const int ndims = 3;
	const int dims[ndims] = {N1,N2,N3};
	kiss_F_fftndr_cfg mycfg=kiss_F_fftndr_alloc(dims,ndims,inverse_fft,NULL,NULL);
	const kiss_F_fft_scalar* input_kissfft = reinterpret_cast<const kiss_F_fft_scalar*>(input.data());
	kiss_F_fft_cpx* output_kissfft = reinterpret_cast<kiss_F_fft_cpx*>(output.data());
	kiss_F_fftndr(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

void perform_3D_FFT_sym_i(const Array<complex<float>,3>& input, Array<float,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft)
{//(Works only for even N3)
 //computes the 3D Fourier transform for a real input array of dimensions N1 x N2 x N3, returns one half of the conjugate-symmetric output array
 //output will be of dimensions N1 x N2 x (N3/2+1)
	const int ndims = 3;
	const int dims[ndims] = {N1,N2,N3};
	kiss_F_fftndr_cfg mycfg=kiss_F_fftndr_alloc(dims,ndims,inverse_fft,NULL,NULL);
	const kiss_F_fft_cpx* input_kissfft = reinterpret_cast<const kiss_F_fft_cpx*>(input.data());
	kiss_F_fft_scalar* output_kissfft = reinterpret_cast<kiss_F_fft_scalar*>(output.data());
	kiss_F_fftndri(mycfg,input_kissfft,output_kissfft);
	free(mycfg); //mycfg is allocated dynamically, so free it
}

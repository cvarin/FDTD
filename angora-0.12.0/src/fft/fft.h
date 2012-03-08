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

#ifndef FFT_H
#define FFT_H

//Declarations of the FFT wrappers

//Whatever FFT is used, it should have the following:
// 1) 2D FFT support (for imaging)
// 2) 3D FFT support (for random-medium generation)
// 3) complex<double> type support at the input
// 4) Different FFT size than the input size


/* 1D Fourier transform
  The blitz++ arrays input and output should have the same length (N_FFT), and should be allocated before calling this routine.
*/
//float and double versions are available
void perform_1D_FFT(const Array<complex<float>,1>& input, Array<complex<float>,1>& output, const int& N_FFT, const int& inverse_fft);
void perform_1D_FFT_sym_r(const Array<float,1>& input, Array<complex<float>,1>& output, const int& N_FFT, const int& inverse_fft);
void perform_1D_FFT_sym_i(const Array<complex<float>,1>& input, Array<float,1>& output, const int& N_FFT, const int& inverse_fft);
void perform_1D_FFT(const Array<complex<double>,1>& input, Array<complex<double>,1>& output, const int& N_FFT, const int& inverse_fft);
void perform_1D_FFT_sym_r(const Array<double,1>& input, Array<complex<double>,1>& output, const int& N_FFT, const int& inverse_fft);
void perform_1D_FFT_sym_i(const Array<complex<double>,1>& input, Array<double,1>& output, const int& N_FFT, const int& inverse_fft);

/* Rearranges the a 1D Fourier-spectrum array, shifting the zero-frequency component with index 0 to the center of the spectrum (similar to the "fftshift" operator of MATLAB)
For even L, the 0th element is shifted to become the (L/2)th element. (with 0-base ordering)
For odd L, the 0th element is shifted to become the ((L-1)/2)th element. (with 0-base ordering)
*/
template <class T>
void fftshift_1D(Array<T,1>& array)
{
    Array<T,1> temp(array.shape());
    int L = array.size();
	//new index of the zero-frequency element (for 0-base ordering)
    int center = L/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    temp(Range(center,L-1)) = array(Range(0,L-center-1));
    temp(Range(0,center-1)) = array(Range(L-center,L-1));
    array = temp;
}

/* Inverts the array resulting from fftshift_1D to its original state.
*/
template <class T>
void ifftshift_1D(Array<T,1>& fftshifted_array)
{
    Array<T,1> temp(fftshifted_array.shape());
    int L = fftshifted_array.size();
	//index of the zero-frequency element (for 0-base ordering)
    int center = L/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    temp(Range(0,L-center-1)) = fftshifted_array(Range(center,L-1));
    temp(Range(L-center,L-1)) = fftshifted_array(Range(0,center-1));
    fftshifted_array = temp;
}

/* 2D Fourier transform
  The blitz++ arrays input and output should have the same sizes (N1xN2), and should be allocated before calling this routine.
*/
//float and double versions are available
void perform_2D_FFT(const Array<complex<float>,2>& input, Array<complex<float>,2>& output, const int& N1, const int& N2, const int& inverse_fft);
void perform_2D_FFT_sym_r(const Array<float,2>& input, Array<complex<float>,2>& output, const int& N1, const int& N2, const int& inverse_fft);
void perform_2D_FFT_sym_i(const Array<complex<float>,2>& input, Array<float,2>& output, const int& N1, const int& N2, const int& inverse_fft);
void perform_2D_FFT(const Array<complex<double>,2>& input, Array<complex<double>,2>& output, const int& N1, const int& N2, const int& inverse_fft);
void perform_2D_FFT_sym_r(const Array<double,2>& input, Array<complex<double>,2>& output, const int& N1, const int& N2, const int& inverse_fft);
void perform_2D_FFT_sym_i(const Array<complex<double>,2>& input, Array<double,2>& output, const int& N1, const int& N2, const int& inverse_fft);

/* Rearranges the 4 quadrants of a 2D Fourier-spectrum array, shifting the zero-frequency component with index (0,0) to the center of the spectrum (similar to the "fftshift" operator of MATLAB)
For even L, the 0th element is shifted to become the (L/2)th element. (with 0-base ordering)
For odd L, the 0th element is shifted to become the ((L-1)/2)th element. (with 0-base ordering)
*/
template <class T>
void fftshift_2D(Array<T,2>& array)
{
    Array<T,2> temp(array.shape());
    int L1 = array.extent(firstDim);
    int L2 = array.extent(secondDim);
	//new indices of the zero-frequency elements (for 0-base ordering)
    int center1 = L1/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center2 = L2/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    //swap first dimension
    temp(Range(center1,L1-1),Range::all()) = array(Range(0,L1-center1-1),Range::all());
    temp(Range(0,center1-1),Range::all()) = array(Range(L1-center1,L1-1),Range::all());
    array = temp;
    //swap second dimension
    temp(Range::all(),Range(center2,L2-1)) = array(Range::all(),Range(0,L2-center2-1));
    temp(Range::all(),Range(0,center2-1)) = array(Range::all(),Range(L2-center2,L2-1));
    array = temp;
}

/* Inverts the array resulting from fftshift_1D to its original state.
*/
template <class T>
void ifftshift_2D(Array<T,2>& fftshifted_array)
{
    Array<T,2> temp(fftshifted_array.shape());
    int L1 = fftshifted_array.extent(firstDim);
    int L2 = fftshifted_array.extent(secondDim);
	//new indices of the zero-frequency elements (for 0-base ordering)
    int center1 = L1/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center2 = L2/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    //swap first dimension
    temp(Range(0,L1-center1-1),Range::all()) = fftshifted_array(Range(center1,L1-1),Range::all());
    temp(Range(L1-center1,L1-1),Range::all()) = fftshifted_array(Range(0,center1-1),Range::all());
    fftshifted_array = temp;
    //swap second dimension
    temp(Range::all(),Range(0,L2-center2-1)) = fftshifted_array(Range::all(),Range(center2,L2-1));
    temp(Range::all(),Range(L2-center2,L2-1)) = fftshifted_array(Range::all(),Range(0,center2-1));
    fftshifted_array = temp;
}

/* Computes the 3D Fourier transform
  The blitz++ arrays input and output should have the same sizes (N1xN2xN3), and should be allocated before calling this routine.
*/
//float and double versions are available
void perform_3D_FFT(const Array<complex<float>,3>& input, Array<complex<float>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);
void perform_3D_FFT_sym_r(const Array<float,3>& input, Array<complex<float>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);
void perform_3D_FFT_sym_i(const Array<complex<float>,3>& input, Array<float,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);
void perform_3D_FFT(const Array<complex<double>,3>& input, Array<complex<double>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);
void perform_3D_FFT_sym_r(const Array<double,3>& input, Array<complex<double>,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);
void perform_3D_FFT_sym_i(const Array<complex<double>,3>& input, Array<double,3>& output, const int& N1, const int& N2, const int& N3, const int& inverse_fft);

/* 3D version of the fftshift_2D above. (generalization is straightforward: see fftshift_2D for documentation)
*/
template <class T>
void fftshift_3D(Array<T,3>& array)
{
    Array<T,3> temp(array.shape());
    int L1 = array.extent(firstDim);
    int L2 = array.extent(secondDim);
    int L3 = array.extent(thirdDim);
	//new indices of the zero-frequency elements (for 0-base ordering)
    int center1 = L1/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center2 = L2/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center3 = L3/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    //swap first dimension
    temp(Range(center1,L1-1),Range::all(),Range::all()) = array(Range(0,L1-center1-1),Range::all(),Range::all());
    temp(Range(0,center1-1),Range::all(),Range::all()) = array(Range(L1-center1,L1-1),Range::all(),Range::all());
    array = temp;
    //swap second dimension
    temp(Range::all(),Range(center2,L2-1),Range::all()) = array(Range::all(),Range(0,L2-center2-1),Range::all());
    temp(Range::all(),Range(0,center2-1),Range::all()) = array(Range::all(),Range(L2-center2,L2-1),Range::all());
    array = temp;
    //swap third dimension
    temp(Range::all(),Range::all(),Range(center3,L3-1)) = array(Range::all(),Range::all(),Range(0,L3-center3-1));
    temp(Range::all(),Range::all(),Range(0,center3-1)) = array(Range::all(),Range::all(),Range(L3-center3,L3-1));
    array = temp;
}

/* Inverts the array resulting from fftshift_3D to its original state.
*/
template <class T>
void ifftshift_3D(Array<T,3>& fftshifted_array)
{
    Array<T,3> temp(fftshifted_array.shape());
    int L1 = fftshifted_array.extent(firstDim);
    int L2 = fftshifted_array.extent(secondDim);
    int L3 = fftshifted_array.extent(thirdDim);
	//new indices of the zero-frequency elements (for 0-base ordering)
    int center1 = L1/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center2 = L2/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L
    int center3 = L3/2; //integer division: (L/2) for even L, ((L-1)/2) for odd L

    //swap first dimension
    temp(Range(0,L1-center1-1),Range::all(),Range::all()) = fftshifted_array(Range(center1,L1-1),Range::all(),Range::all());
    temp(Range(L1-center1,L1-1),Range::all(),Range::all()) = fftshifted_array(Range(0,center1-1),Range::all(),Range::all());
    fftshifted_array = temp;
    //swap second dimension
    temp(Range::all(),Range(0,L2-center2-1),Range::all()) = fftshifted_array(Range::all(),Range(center2,L2-1),Range::all());
    temp(Range::all(),Range(L2-center2,L2-1),Range::all()) = fftshifted_array(Range::all(),Range(0,center2-1),Range::all());
    fftshifted_array = temp;
    //swap third dimension
    temp(Range::all(),Range::all(),Range(0,L3-center3-1)) = fftshifted_array(Range::all(),Range::all(),Range(center3,L3-1));
    temp(Range::all(),Range::all(),Range(L3-center3,L3-1)) = fftshifted_array(Range::all(),Range::all(),Range(0,center3-1));
    fftshifted_array = temp;
}

#endif

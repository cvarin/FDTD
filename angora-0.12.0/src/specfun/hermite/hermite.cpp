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

//calculation of the Hermite polynomials and polynomial coefficients

#include "headers.h"


double Hermite(const int& N, const double& x)
{
	double x1,x2,x3;

	x1 = 1;
	x2 = 2*x;

	if (N==0)
	{
		x3 = x1;
		return x3;
	}
	else if (N==1)
	{
		x3 = x2;
		return x3;
	}

	for(int i=2; i<=N; i++)
	{
		x3 = 2*x*x2 - 2*(i-1)*x1;
		x1 = x2;
		x2 = x3;
	}

	return x3;
}

void HermiteCoeff(const int& N, Array<int,1>& coeffs)
{
	/*  Hn(x) = coeff(0) + coeff(1)*x + ... + coeff(N)*x^N   */
	/*
 	H0(x) = 1
 	H1(x) = 2x
 	H2(x) = -2+4x^2
 	H3(x) = -12x+8x^3
 	H4(x) = 12-48x^2+16x^4
 	*/

    coeffs.resize(Range(0,N));
    for (int i=0;i<=N;i++)
    {
        coeffs(i) = 0;
    }
    coeffs(N) = (int)exp(N*log(2.0));
    for (int i=0;i<=N/2-1; i++)
    {
        coeffs(N-2*(i+1)) = -coeffs(N-2*i)*(N-2*i)*(N-2*i-1)/(4*(i+1));
    }
}

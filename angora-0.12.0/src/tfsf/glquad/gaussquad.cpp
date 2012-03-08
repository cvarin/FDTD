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

//Computes the Gaussian quadrature positions and weights between [-1,1]

#include "headers.h"

#include "gaussquad.h"


void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w)
{
	//uses Newton's rule with initial guess cos(pi*(4i+3)/(4n+2)), i=0..N-1
	// (Numerical Recipes, pp 152)

	//initialize blitz++ arrays
	x.resize(n);
	w.resize(n);
	//initialize to the values for n=1
	x=0;
	w=2;

	if (n>1)
	{
		//Roots r_i (i=0:n-1) are calculated backwards from 1 to -1
		//If i'th root is r, (n-1-i)th root is (-r)
		//Only half of them have to be calculated:
		// even n: n/2 roots are needed (i=0:n/2-1)
		// odd n: (n+1)/2 roots are needed (i=0:(n+1)/2-1)
		//        [although 0 is surely a root, calculation is required for computing the weight at 0]

		double r; //a root of the legendre polynomial
		double r_next; //next estimate in the Newton iteration
		double pj_2r,pj_1r,pjr; //recursion variables for computing P_n(r) [respectively, P_(j-2)(r),P_(j-1)(r),P_j(r)]
		double dpr; //derivative of P_n(r)

		for (int i=0; i<=(n+1)/2-1 /*amounts to n/2-1 for even n*/; i++)
		{
			//calculate the i'th root r_i
			//initial guess:
			r_next=cos(M_PI*(4*i+3.0)/(4*n+2.0));
			do
			{//start Newton iteration
				//calculate P_n(r) and the derivative P_n'(r) for Newton's algorithm
				r=r_next;
				pj_1r = 1;
				pjr = r;
				for (int j=2; j<=n ; j++)
				{//calculate P_n(r)
					pj_2r=pj_1r;
					pj_1r=pjr;
					pjr = (2*j-1.0)/j*r*pj_1r - (j-1.0)/j*pj_2r;
				}
				//now, pjr is P_n(r), pj_1r is P_(n-1)(r)
				dpr = (r*pjr - pj_1r)*n/(r*r-1); //recursion relation for P_n'(r)
				r_next = r - pjr/dpr; //next Newton estimate
			} while (abs(r_next-r)>LIBSTD_DBL_EPSILON*100); //7 decimal digits should be enough for double
			x(i) = r;
			//If i'th root is r, (n-1-i)th root is (-r)
			x(n-1-i) = -r;
			w(i) = 2.0/((1-r*r)*dpr*dpr);
			//weights are symmetric around 0
			w(n-1-i) = w(i);
		}
	}
// cout << x << endl;
// cout << w << endl;
}

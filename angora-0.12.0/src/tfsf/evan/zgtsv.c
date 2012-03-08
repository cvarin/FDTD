// Copyright (c) 1992-2010 The University of Tennessee and The University
//                         of Tennessee Research Foundation.  All rights
//                         reserved.
// Copyright (c) 2000-2010 The University of California Berkeley. All
//                         rights reserved.
// Copyright (c) 2006-2010 The University of Colorado Denver.  All rights
//                         reserved.
// 
// $COPYRIGHT$
// 
// Additional copyrights may follow
// 
// $HEADER$
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
// - Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
// 
// - Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer listed
//   in this license in the documentation and/or other materials
//   provided with the distribution.
// 
// - Neither the name of the copyright holders nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
// 
// The copyright holders provide no reassurances that the source code
// provided does not infringe any patent, copyright, or any other
// intellectual property rights of third parties.  The copyright holders
// disclaim any liability to any recipient for claims brought against
// recipient by any third party for infringement of that parties
// intellectual property rights.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*  (added by Ilker R. Capoglu)  This is the "modified BSD license".
*/

#include "f2c.h"

int zgtsv_(int *n, int *nrhs, 
	double *dl_r, double *dl_i, 
	double *d__r, double *d__i,
	double *du_r, double *du_i, 
	double *b_r, double *b_i, 
	int *ldb, int *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZGTSV  solves the equation   

       A*X = B,   

    where A is an N-by-N tridiagonal matrix, by Gaussian elimination with   
    partial pivoting.   

    Note that the equation  A'*X = B  may be solved by interchanging the   
    order of the arguments DU and DL.   

    Arguments   
    =========   

    N       (input) int   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) int   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input/output) COMPLEX*16 array, dimension (N-1)   
            On entry, DL must contain the (n-1) subdiagonal elements of   
            A.   
            On exit, DL is overwritten by the (n-2) elements of the   
            second superdiagonal of the upper triangular matrix U from   
            the LU factorization of A, in DL(1), ..., DL(n-2).   

    D       (input/output) COMPLEX*16 array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   
            On exit, D is overwritten by the n diagonal elements of U.   

    DU      (input/output) COMPLEX*16 array, dimension (N-1)   
            On entry, DU must contain the (n-1) superdiagonal elements   
            of A.   
            On exit, DU is overwritten by the (n-1) elements of the first   
            superdiagonal of U.   

    B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) int   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) int   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero, and the solution   
                  has not been computed.  The factorization has not been   
                  completed unless i = N.   

    =====================================================================   
*/

    /*   Parameter adjustments */
    /* System generated locals */
    static int b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    static double d__1, d__2, d__3, d__4;
    static double z__1_r, z__1_i, z__2_r, z__2_i, z__3_r, z__3_i, z__4_r, z__4_i, z__5_r, z__5_i;
    /* Builtin functions */
//    double d_imag(doublecomplex *);
    void z_div(double *, double *, double *, double *, double *, double *);
    /* Local variables */
    static double temp_r, temp_i, mult_r, mult_i;
    static int j, k;
//    extern /* Subroutine */ int xerbla_(char *, int *);

#define b_subscr(a_1,a_2) (a_2)*b_dim1 + a_1

	//Fortran arrays have a lower bound of 1, so adjust the pointers accordingly
	--dl_r; --dl_i;
    --d__r; --d__i;
    --du_r; --du_i;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b_r -= b_offset; b_i -= b_offset;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
	} else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
//	xerbla_("ZGTSV ", &i__1);
	return 0;
    }

    if (*n == 0) {
	return 0;
    }

    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if (dl_r[i__2] == 0. && dl_i[i__2] == 0.) {

/*           Subdiagonal is zero, no elimination is required. */

	    i__2 = k;
	    if (d__r[i__2] == 0. && d__i[i__2] == 0.) {

/*              Diagonal is zero: set INFO = K and return; a unique   
                solution can not be found. */

		*info = k;
		return 0;
	    }
	} else /* if(complicated condition) */ {
	    i__2 = k;
	    i__3 = k;
	    if ((d__1 = d__r[i__2], abs(d__1)) + (d__2 = d__i[k], 
		    abs(d__2)) >= (d__3 = dl_r[i__3], abs(d__3)) + (d__4 = 
		    dl_i[k], abs(d__4))) {

/*           No row interchange required */

		z_div(&z__1_r, &z__1_i, &dl_r[k], &dl_i[k], &d__r[k], &d__i[k]);
		mult_r = z__1_r, mult_i = z__1_i;
		i__2 = k + 1;
		i__3 = k + 1;
		i__4 = k;
		z__2_r = mult_r * du_r[i__4] - mult_i * du_i[i__4], z__2_i = 
			mult_r * du_i[i__4] + mult_i * du_r[i__4];
		z__1_r = d__r[i__3] - z__2_r, z__1_i = d__i[i__3] - z__2_i;
		d__r[i__2] = z__1_r, d__i[i__2] = z__1_i;
		i__2 = *nrhs;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = b_subscr(k + 1, j);
		    i__4 = b_subscr(k + 1, j);
		    i__5 = b_subscr(k, j);
		    z__2_r = mult_r * b_r[i__5] - mult_i * b_i[i__5], z__2_i =
			     mult_r * b_i[i__5] + mult_i * b_r[i__5];
		    z__1_r = b_r[i__4] - z__2_r, z__1_i = b_i[i__4] - z__2_i;
		    b_r[i__3] = z__1_r, b_i[i__3] = z__1_i;
/* L10: */
		}
		if (k < *n - 1) {
		    i__2 = k;
		    dl_r[i__2] = 0., dl_i[i__2] = 0.;
		}
	    } else {

/*           Interchange rows K and K+1 */

		z_div(&z__1_r, &z__1_i, &d__r[k], &d__i[k], &dl_r[k], &dl_i[k]);
		mult_r = z__1_r, mult_i = z__1_i;
		i__2 = k;
		i__3 = k;
		d__r[i__2] = dl_r[i__3], d__i[i__2] = dl_i[i__3];
		i__2 = k + 1;
		temp_r = d__r[i__2], temp_i = d__i[i__2];
		i__2 = k + 1;
		i__3 = k;
		z__2_r = mult_r * temp_r - mult_i * temp_i, z__2_i = mult_r * 
			temp_i + mult_i * temp_r;
		z__1_r = du_r[i__3] - z__2_r, z__1_i = du_i[i__3] - z__2_i;
		d__r[i__2] = z__1_r, d__i[i__2] = z__1_i;
		if (k < *n - 1) {
		    i__2 = k;
		    i__3 = k + 1;
		    dl_r[i__2] = du_r[i__3], dl_i[i__2] = du_i[i__3];
		    i__2 = k + 1;
		    z__2_r = -mult_r, z__2_i = -mult_i;
		    i__3 = k;
		    z__1_r = z__2_r * dl_r[i__3] - z__2_i * dl_i[i__3], 
			    z__1_i = z__2_r * dl_i[i__3] + z__2_i * dl_r[i__3]
			    ;
		    du_r[i__2] = z__1_r, du_i[i__2] = z__1_i;
		}
		i__2 = k;
		du_r[i__2] = temp_r, du_i[i__2] = temp_i;
		i__2 = *nrhs;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = b_subscr(k, j);
		    temp_r = b_r[i__3], temp_i = b_i[i__3];
		    i__3 = b_subscr(k, j);
		    i__4 = b_subscr(k + 1, j);
		    b_r[i__3] = b_r[i__4], b_i[i__3] = b_i[i__4];
		    i__3 = b_subscr(k + 1, j);
		    i__4 = b_subscr(k + 1, j);
		    z__2_r = mult_r * b_r[i__4] - mult_i * b_i[i__4], z__2_i =
			     mult_r * b_i[i__4] + mult_i * b_r[i__4];
		    z__1_r = temp_r - z__2_r, z__1_i = temp_i - z__2_i;
		    b_r[i__3] = z__1_r, b_i[i__3] = z__1_i;
/* L20: */
		}
	    }
	}
/* L30: */
    }
    i__1 = *n;
    if (d__r[i__1] == 0. && d__i[i__1] == 0.) {
	*info = *n;
	return 0;
    }

/*     Back solve with the matrix U from the factorization. */

    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
	i__2 = b_subscr(*n, j);
	z_div(&z__1_r, &z__1_i, &b_r[b_subscr(*n,j)], &b_i[b_subscr(*n,j)], &d__r[*n], &d__i[*n]);
	b_r[i__2] = z__1_r, b_i[i__2] = z__1_i;
	if (*n > 1) {
	    i__2 = b_subscr(*n - 1, j);
	    i__3 = b_subscr(*n - 1, j);
	    i__4 = *n - 1;
	    i__5 = b_subscr(*n, j);
	    z__3_r = du_r[i__4] * b_r[i__5] - du_i[i__4] * b_i[i__5], z__3_i =
		     du_r[i__4] * b_i[i__5] + du_i[i__4] * b_r[i__5];
	    z__2_r = b_r[i__3] - z__3_r, z__2_i = b_i[i__3] - z__3_i;
	    z_div(&z__1_r, &z__1_i, &z__2_r, &z__2_i, &d__r[*n - 1], &d__i[*n - 1]);
	    b_r[i__2] = z__1_r, b_i[i__2] = z__1_i;
	}
	for (k = *n - 2; k >= 1; --k) {
	    i__2 = b_subscr(k, j);
	    i__3 = b_subscr(k, j);
	    i__4 = k;
	    i__5 = b_subscr(k + 1, j);
	    z__4_r = du_r[i__4] * b_r[i__5] - du_i[i__4] * b_i[i__5], z__4_i =
		     du_r[i__4] * b_i[i__5] + du_i[i__4] * b_r[i__5];
	    z__3_r = b_r[i__3] - z__4_r, z__3_i = b_i[i__3] - z__4_i;
	    i__6 = k;
	    i__7 = b_subscr(k + 2, j);
	    z__5_r = dl_r[i__6] * b_r[i__7] - dl_i[i__6] * b_i[i__7], z__5_i =
		     dl_r[i__6] * b_i[i__7] + dl_i[i__6] * b_r[i__7];
	    z__2_r = z__3_r - z__5_r, z__2_i = z__3_i - z__5_i;
	    z_div(&z__1_r, &z__1_i, &z__2_r, &z__2_i, &d__r[k], &d__i[k]);
	    b_r[i__2] = z__1_r, b_i[i__2] = z__1_i;
/* L40: */
	}
/* L50: */
    }

    return 0;

/*     End of ZGTSV */

} /* zgtsv_ */

#undef b_subscr

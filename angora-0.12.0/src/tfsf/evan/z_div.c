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


/*  (added by Ilker R. Capoglu) This is the "modified BSD license".
*/

#include "f2c.h"

void z_div(double *c_r, double *c_i, double *a_r, double *a_i, double *b_r, double *b_i)
{
	static double ratio, den;
	static double abr, abi, cr;

	if( (abr = *b_r) < 0.)
		abr = - abr;
	if( (abi = *b_i) < 0.)
		abi = - abi;
	if( abr <= abi )
		{
		if(abi == 0){
			printf("complex division by zero\n");
			exit(0);
		}
		ratio = *b_r / *b_i ;
		den = *b_i * (1 + ratio*ratio);
		cr = (*a_r*ratio + *a_i) / den;
		*c_i = (*a_i*ratio - *a_r) / den;
		}

	else
		{
		ratio = *b_i / *b_r ;
		den = *b_r * (1 + ratio*ratio);
		cr = (*a_r + *a_i*ratio) / den;
		*c_i = (*a_i - *a_r*ratio) / den;
		}
	*c_r = cr;
}

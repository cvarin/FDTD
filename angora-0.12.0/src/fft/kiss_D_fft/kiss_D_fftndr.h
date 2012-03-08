#ifndef KISS_D_NDR_H
#define KISS_D_NDR_H

#include "kiss_D_fft.h"
#include "kiss_D_fftr.h"
#include "kiss_D_fftnd.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct kiss_D_fftndr_state *kiss_D_fftndr_cfg;


kiss_D_fftndr_cfg  kiss_D_fftndr_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem);
/*
 dims[0] must be even

 If you don't care to allocate space, use mem = lenmem = NULL
*/


void kiss_D_fftndr(
        kiss_D_fftndr_cfg cfg,
        const kiss_D_fft_scalar *timedata,
        kiss_D_fft_cpx *freqdata);
/*
 input timedata has dims[0] X dims[1] X ... X  dims[ndims-1] scalar points
 output freqdata has dims[0] X dims[1] X ... X  dims[ndims-1]/2+1 complex points
*/

void kiss_D_fftndri(
        kiss_D_fftndr_cfg cfg,
        const kiss_D_fft_cpx *freqdata,
        kiss_D_fft_scalar *timedata);
/*
 input and output dimensions are the exact opposite of kiss_D_fftndr
*/


#define kiss_D_fftr_free free

#ifdef __cplusplus
}
#endif

#endif

#ifndef KISS_D_FTR_H
#define KISS_D_FTR_H

#include "kiss_D_fft.h"
#ifdef __cplusplus
extern "C" {
#endif


/*

 Real optimized version can save about 45% cpu time vs. complex fft of a real seq.



 */

typedef struct kiss_D_fftr_state *kiss_D_fftr_cfg;


kiss_D_fftr_cfg kiss_D_fftr_alloc(int nfft,int inverse_fft,void * mem, size_t * lenmem);
/*
 nfft must be even

 If you don't care to allocate space, use mem = lenmem = NULL
*/


void kiss_D_fftr(kiss_D_fftr_cfg cfg,const kiss_D_fft_scalar *timedata,kiss_D_fft_cpx *freqdata);
/*
 input timedata has nfft scalar points
 output freqdata has nfft/2+1 complex points
*/

void kiss_D_fftri(kiss_D_fftr_cfg cfg,const kiss_D_fft_cpx *freqdata,kiss_D_fft_scalar *timedata);
/*
 input freqdata has  nfft/2+1 complex points
 output timedata has nfft scalar points
*/

#define kiss_D_fftr_free free

#ifdef __cplusplus
}
#endif
#endif

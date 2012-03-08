#ifndef KISS_D_FFTND_H
#define KISS_D_FFTND_H

#include "kiss_D_fft.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct kiss_D_fftnd_state * kiss_D_fftnd_cfg;

kiss_D_fftnd_cfg  kiss_D_fftnd_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem);
void kiss_D_fftnd(kiss_D_fftnd_cfg  cfg,const kiss_D_fft_cpx *fin,kiss_D_fft_cpx *fout);

#ifdef __cplusplus
}
#endif
#endif

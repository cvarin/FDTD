#ifndef KISS_F_FFTND_H
#define KISS_F_FFTND_H

#include "kiss_F_fft.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct kiss_F_fftnd_state * kiss_F_fftnd_cfg;

kiss_F_fftnd_cfg  kiss_F_fftnd_alloc(const int *dims,int ndims,int inverse_fft,void*mem,size_t*lenmem);
void kiss_F_fftnd(kiss_F_fftnd_cfg  cfg,const kiss_F_fft_cpx *fin,kiss_F_fft_cpx *fout);

#ifdef __cplusplus
}
#endif
#endif

#ifndef misc
#define misc

#include "mslib.h"

double wavelength( double kev );
void mat2grey(float *mat, int32_t *grey, uint32_t numPixel);
uint32_t max(float *a, int32_t n, int32_t i, int32_t j, int32_t k, uint32_t sizeBlock, uint32_t ref);
void heapdown(float *a, int32_t n, int32_t i, uint32_t sizeBlock, uint32_t ref);
void heapsort(float *a, int32_t n, uint32_t sizeBlock, uint32_t ref);
void displayParam(expPara *param);
void *malloc1D( int n, size_t size);
void **malloc2D( int nx, int ny, size_t size);
void ***malloc3D( int nx, int ny, int nz, size_t size);
void createFparams(double **fparams);

#endif // misc


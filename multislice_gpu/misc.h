#ifndef misc
#define misc

#include "mslib.h"

void displayParam(expPara *param);
void *malloc1D( int n, size_t size);
void **malloc2D( int nx, int ny, size_t size);
void ***malloc3D( int nx, int ny, int nz, size_t size);
void createFparams(double **fparams);

#endif // misc



#ifndef mslib
#define mslib

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define MAX_ELEMENT_Z 103
#define MAX_NUM_ELEMENT 5 // max number of element type in the material
#define NUM_RADIUS_SAMPLE 100 // sample point of the atomic potential look up table

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

struct expPara {
    uint32_t nx;//(pixel) size of output image
    uint32_t ny;//(pixel) size of output image
    uint32_t ncell_x;
    uint32_t ncell_y;
    uint32_t ncell_z;
    uint32_t numAtomUnitCell;
    float cell_a;
    float cell_b;
    float cell_c;
    float supercell_a;//(Angstrom) a size of the whole lattice
    float supercell_b;//(Angstrom) b size of the whole lattice
    float supercell_c;//(Angstrom) c size of the whole lattice
    float deltaZ;//(Angstrom) thickness of one slice
    float electronV;// electron source voltage
    float df0;//(Angstrom) defocus
    float sigmaf;
    float temperature;//(K) sample temperature
    float wavelen;
    uint32_t totalNumAtom;
    uint32_t numElement;
    float scale_x;
    float scale_y;
    float limit_x;
    float limit_y;
    float radMin_sq;
    float radMax_sq;
    uint32_t numSlice;
    cufftHandle fft_plan;

};

struct atomVZ_LUT {
    double splineCoeff[MAX_NUM_ELEMENT][3][NUM_RADIUS_SAMPLE];//
    double spline_x[NUM_RADIUS_SAMPLE];
    double spline_y[MAX_NUM_ELEMENT][NUM_RADIUS_SAMPLE];
};

struct atomVZ_LUT_all {
    double splineCoeff[MAX_ELEMENT_Z][3][NUM_RADIUS_SAMPLE];//
    double spline_x[NUM_RADIUS_SAMPLE];
    double spline_y[MAX_ELEMENT_Z][NUM_RADIUS_SAMPLE];
};

void constructSupercell(FILE *fpCell, float *atomSupercell, float *atomUnitCell, expPara *param, atomVZ_LUT_all *atomVZ_all, atomVZ_LUT *atomVZ);
__global__ void vzRealSpace(uint32_t numAtom, uint32_t startPoint, float *d_atomSupercell, float *d_projPotential);
void multislice_run(float *atomSupercell, expPara *param, atomVZ_LUT *elementVZ, uint32_t *startPointList, uint32_t *numAtomList);
void sliceSupercell(float *atomSupercell, expPara *param, uint32_t *startPointList, uint32_t *numAtomList);
void generateSplineCoeff(atomVZ_LUT_all *atomVZ_all);
float vzatom(uint32_t Z, float radius, double **fparams);
void splinh( double *x, double *y, double *b, double *c, double *d, uint32_t n);
double bessk0( double x );
double bessi0( double x );
__device__ double seval( double *x, double *y, double *b, double *c, double *d, int n, double x0 );
__global__ void initialiseComplex(uint32_t nx, uint32_t ny, float2 *array, float reIni, float imIni);
__global__ void createKArray_1D(uint32_t n, float dim, float *k);
__global__ void createKArray_2D(uint32_t nx, uint32_t ny, float *d_kx, float *d_ky, float *d_k2);
__global__ void createProp_1D(uint32_t n, float *k, float scale, float wavelen, float2 *prop);
__global__ void createProp_2D(uint32_t nx, float2 *d_propx, float2 *d_propy, float2 *d_propxy);


#endif


#ifndef mslib
#define mslib

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define MAX_ELEMENT_Z 103
#define NUM_RADIUS_SAMPLE 100

struct msPara {
    uint32_t nx;//(pixel) size of output image
    uint32_t ny;//(pixel) size of output image
    float supercell_a;//(Angstrom) a size of the whole lattice
    float supercell_b;//(Angstrom) b size of the whole lattice
    float deltaZ;//(Angstrom) thickness of one slice
    float electronV;// electron source voltage
    float df0;//(Angstrom) defocus
    float sigmaf;
    float temperature;//(K) sample temperature
    uint32_t totalNumAtoms;


};

struct sfInfo {
    char isElement[MAX_ELEMENT_Z];// store the existence of certain elements
    uint32_t nAtomType;
    double splineCoeff[MAX_ELEMENT_Z][3][NUM_RADIUS_SAMPLE];//
    double spline_x[NUM_RADIUS_SAMPLE];
    double spline_y[MAX_ELEMENT_Z][NUM_RADIUS_SAMPLE];
};

void constructSupercell(FILE *fpCell, float *atomSupercell, msPara *para, sfInfo *info);
__global__ void vzRealSpace(uint32_t numAtoms, sfInfo *info, float *atomSupercell, float *projPotential, uint32_t nx, uint32_t ny, float scale_x, float scale_y, uint32_t limit_x, uint32_t limit_y, float radMin_sq, float radMax_sq);
void calculatePhaseGrating(float *atomSupercell, msPara *para, sfInfo *info);
void generateSplineCoeff(sfInfo *info);
float vzatom(uint32_t Z, float radius, double **fparams);
void splinh( double *x, double *y, double *b, double *c, double *d, uint32_t n);
double bessk0( double x );
double bessi0( double x );
__device__ double seval( double *x, double *y, double *b, double *c,
         double *d, int n, double x0 );

#endif

#include "mslib.h"
#include "misc.h"
#include <cuda.h>
#include <cuda_profiler_api.h>


#define MAX_ATOM_RADIUS 3.0 //max atomic radius in the unit of Angstrom
#define MIN_ATOM_RADIUS 0.01 //min atomic radius in the unit of Angstrom
#define NUM_LORENZ 3
#define NUM_GAUSS 3
#define SMALL 1.0e-25
#define NUM_FPARAM 12 // number of parameters for calculating the scattering factor

void constructSupercell(FILE *fpCell, float *atomSupercell, float *atomUnitCell, expPara *param, atomVZ_LUT_all *atomVZ_all, atomVZ_LUT *atomVZ) {

    ///read the atomic number, x, y, z, occupation and thermal displacement
    uint32_t elementCount = 0;
    uint8_t elementList[MAX_NUM_ELEMENT] = {0};
    uint8_t currAtomZ = 255;//this is a impossible Z for existing atoms
    uint8_t elementExist;


    for (uint32_t readIndex = 0; readIndex < param->numAtomUnitCell; readIndex ++) {
        fscanf(fpCell, "%d %f %f %f %f %f", &currAtomZ, atomUnitCell+readIndex*6+1, atomUnitCell+readIndex*6+2, atomUnitCell+readIndex*6+3, atomUnitCell+readIndex*6+4, atomUnitCell+readIndex*6+5);
        printf("%2d %1.5f %1.5f %1.5f %1.1f %1.3f\n", currAtomZ, atomUnitCell[readIndex*6+1], atomUnitCell[readIndex*6+2], atomUnitCell[readIndex*6+3], atomUnitCell[readIndex*6+4], atomUnitCell[readIndex*6+5]);

        elementExist = 0;


        for (uint32_t i = 0; i < MAX_NUM_ELEMENT; i++){
            if (currAtomZ == elementList[i]) {
                elementExist = 1;
                atomUnitCell[readIndex*6] = (float)i;
            }

        }

        if (elementExist == 0) {

            elementList[elementCount] = currAtomZ;
            atomUnitCell[readIndex*6] = (float)elementCount;

            ///construct the look up table for atomic potentials of existing elements
            memcpy(atomVZ->splineCoeff[elementCount], atomVZ_all->splineCoeff[currAtomZ], 3*NUM_RADIUS_SAMPLE*sizeof(double));
            memcpy(atomVZ->spline_y[elementCount], atomVZ_all->spline_y[currAtomZ], NUM_RADIUS_SAMPLE*sizeof(double));

            elementCount++;
        }



    }


    memcpy(atomVZ->spline_x, atomVZ_all->spline_x, NUM_RADIUS_SAMPLE*sizeof(double));
    param->numElement = elementCount;

    printf("initialize the supercell\n");
    memcpy(atomSupercell, atomUnitCell, param->numAtomUnitCell*6*sizeof(float));

    ///replicate the unit cell on all three directions

    uint32_t currNAtoms = param->numAtomUnitCell;

    if( param->ncell_x > 1 ) {
        for (uint32_t i=1; i<param->ncell_x; i++) {
            for (uint32_t j=0; j<currNAtoms; j++) {
                *(atomSupercell + currNAtoms*i*6 + j*6 + 0) = *(atomSupercell + j*6 + 0);
                *(atomSupercell + currNAtoms*i*6 + j*6 + 1) = *(atomSupercell + j*6 + 1) + i*param->cell_a;
                *(atomSupercell + currNAtoms*i*6 + j*6 + 2) = *(atomSupercell + j*6 + 2);
                *(atomSupercell + currNAtoms*i*6 + j*6 + 3) = *(atomSupercell + j*6 + 3);
                *(atomSupercell + currNAtoms*i*6 + j*6 + 4) = *(atomSupercell + j*6 + 4);
                *(atomSupercell + currNAtoms*i*6 + j*6 + 5) = *(atomSupercell + j*6 + 5);

            }
        }
        currNAtoms = currNAtoms*param->ncell_x;
        //*ax = (*ax) * ncellx;
    }

    if( param->ncell_y > 1 ) {
        for(uint32_t i=1; i<param->ncell_y; i++)
        for(uint32_t j=0; j<currNAtoms; j++) {
            *(atomSupercell + currNAtoms*i*6 + j*6 + 0) = *(atomSupercell + j*6 + 0);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 1) = *(atomSupercell + j*6 + 1);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 2) = *(atomSupercell + j*6 + 2) + i*param->cell_b;
            *(atomSupercell + currNAtoms*i*6 + j*6 + 3) = *(atomSupercell + j*6 + 3);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 4) = *(atomSupercell + j*6 + 4);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 5) = *(atomSupercell + j*6 + 5);

        }
        currNAtoms = currNAtoms*param->ncell_y;
        //*by = (*by) * ncelly;
    }

    if( param->ncell_z > 1 ) {
        for(uint32_t i=1; i<param->ncell_z; i++)
        for(uint32_t j=0; j<currNAtoms; j++) {
            *(atomSupercell + currNAtoms*i*6 + j*6 + 0) = *(atomSupercell + j*6 + 0);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 1) = *(atomSupercell + j*6 + 1);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 2) = *(atomSupercell + j*6 + 2);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 3) = *(atomSupercell + j*6 + 3) + i*param->cell_c;
            *(atomSupercell + currNAtoms*i*6 + j*6 + 4) = *(atomSupercell + j*6 + 4);
            *(atomSupercell + currNAtoms*i*6 + j*6 + 5) = *(atomSupercell + j*6 + 5);
        }
        currNAtoms = currNAtoms*param->ncell_z;
        //*cz = (*cz) * ncellz;
    }

    ///check the constructed supercell
    printf("check the supercell.\n");
    for (uint32_t readIndex = 0; readIndex < param->totalNumAtom; readIndex ++) {
        for (uint32_t i = 0; i < 6; i++){
            printf("%f ", *(atomSupercell + readIndex*6 + i));
        }
        printf("\n");

    }

}

///store the experimental parameters and the look up table of elements in the constant memory of the device
__constant__ expPara d_param;
__constant__ atomVZ_LUT d_elementVZ;


__global__ void vzRealSpace(uint32_t numAtom, uint32_t startPoint, float *d_atomSupercell, float *d_projPotential) {

    uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y;
	uint32_t endPoint = startPoint + numAtom;

    for (uint32_t i = startPoint; i < endPoint; i++) {

        uint32_t elementIndex = (uint32_t)*(d_atomSupercell + i*6 + 0);
        float origin_x = *(d_atomSupercell + i*6 + 1);
        float origin_y = *(d_atomSupercell + i*6 + 2);
        float occupy = *(d_atomSupercell + i*6 + 4);

        float distance_x = (ix * d_param.scale_x - origin_x);
        float distance_y = (iy * d_param.scale_y - origin_y);
        float rsq = distance_x * distance_x + distance_y * distance_y;

        double *spline_x = d_elementVZ.spline_x;
        double *spline_y = d_elementVZ.spline_y[elementIndex];
        double *spline_b = d_elementVZ.splineCoeff[elementIndex][0];
        double *spline_c = d_elementVZ.splineCoeff[elementIndex][1];
        double *spline_d = d_elementVZ.splineCoeff[elementIndex][2];

        if( rsq < d_param.radMin_sq ) rsq = d_param.radMin_sq;

        if (rsq < d_param.radMax_sq) {
            float vz = occupy * seval(spline_x, spline_y, spline_b, spline_c, spline_d, NUM_RADIUS_SAMPLE,(double)rsq);
            *(d_projPotential + d_param.ny * ix + iy) += (float)vz;
        }


    }

}


void multislice_run(float *atomSupercell, expPara *param, atomVZ_LUT *elementVZ, uint32_t *startPointList, uint32_t *numAtomList) {

    ///initialize some useful parameters
    param->scale_x = param->supercell_a/param->nx;
    param->scale_y = param->supercell_b/param->ny;

    float radMin = 0.25 * sqrt(0.5*(param->scale_x*param->scale_x + param->scale_y*param->scale_y));
    param->radMin_sq = radMin*radMin;
    param->radMax_sq = MAX_ATOM_RADIUS*MAX_ATOM_RADIUS;

    param->limit_x = (uint32_t)(MAX_ATOM_RADIUS/param->scale_x) + 1;
    param->limit_y = (uint32_t)(MAX_ATOM_RADIUS/param->scale_y) + 1;

    float scale = 1;/// this the scaling factor used when calculating the propagation function

    ///initialize the device transmission and propagation functions
    cudaError_t err = cudaSuccess;
    dim3 dimBlock( 16, 16);
	dim3 dimGrid( param->nx / dimBlock.x, param->ny / dimBlock.y);

    float2 *d_trans, *d_wave;
    cudaMalloc((void **)&d_trans, param->nx*param->ny*sizeof(float2));
    cudaMalloc((void **)&d_wave, param->nx*param->ny*sizeof(float2));
    initialise<<<dimGrid, dimBlock>>>(param->nx, param->ny, d_trans, 0.0F, 0.0F);
    initialise<<<dimGrid, dimBlock>>>(param->nx, param->ny, d_wave, 1.0F, 0.0F);

    float *d_kx, *d_ky, *d_k2;
    cudaMalloc((void **)&d_kx, param->nx*sizeof(float2));
    cudaMalloc((void **)&d_ky, param->ny*sizeof(float2));
    cudaMalloc((void **)&d_k2, param->nx*param->ny*sizeof(float2));
    createKArray_1D<<<dimGrid.x, dimBlock.x>>>(param->nx, param->supercell_a, d_kx);
    createKArray_1D<<<dimGrid.y, dimBlock.y>>>(param->ny, param->supercell_b, d_ky);
    createKArray_2D<<<dimGrid, dimBlock>>>(param->nx, param->ny, d_kx, d_ky, d_k2);

    float2 *d_propx, *d_propy, *d_propxy;
    cudaMalloc((void **)&d_propx, param->nx*sizeof(float2));
    cudaMalloc((void **)&d_propy, param->ny*sizeof(float2));
    cudaMalloc((void **)&d_propxy, param->nx*param->ny*sizeof(float2));
    createProp_1D<<<dimGrid.x, dimBlock.x>>>(param->nx, d_kx, scale, param->wavelen, d_propx);
    createProp_1D<<<dimGrid.y, dimBlock.y>>>(param->ny, d_ky, scale, param->wavelen, d_propy);
    createProp_2D<<<dimGrid, dimBlock>>>(param->nx, d_propx, d_propy, d_propxy);

    ///copy the atom position information into the device
    float *d_atomSupercell = NULL;
    err = cudaMalloc((void **)&d_atomSupercell, param->totalNumAtom*sizeof(float)*6);
    err = cudaMemcpy(d_atomSupercell, atomSupercell, param->totalNumAtom*sizeof(float)*6, cudaMemcpyHostToDevice);



    ///initialize the projected potential
    float *h_projPotential = (float *)malloc(param->nx*param->ny*sizeof(float));
    float *d_projPotential = NULL;

     for (uint32_t i = 0; i < param->nx*param->ny; i++) {
        h_projPotential[i] = 0;
    }

    err = cudaMalloc((void **)&d_projPotential, param->nx*param->ny*sizeof(float));
    err = cudaMemcpy(d_projPotential, h_projPotential, param->nx*param->ny*sizeof(float), cudaMemcpyHostToDevice);

    ///copy the parameters and LUT into the constant memory of the device
    cudaMemcpyToSymbol(d_param, param, sizeof(expPara), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_elementVZ, elementVZ, sizeof(atomVZ_LUT), cudaMemcpyHostToDevice);

    ///perform the multislice calculation
    for (uint32_t i = 0; i < param->numSlice; i++) {
        ///number of atoms in one slice, may be changed
        printf("total number of atoms in this slice is %d.\n",numAtomList[i]);

        ///use the real space method to calculate the projected potential
        vzRealSpace<<<dimGrid, dimBlock>>>(numAtomList[i], startPointList[i], d_atomSupercell, d_projPotential);

        ///transmit and propagate the wave


    }



    err = cudaGetLastError();
    err = cudaFree(d_trans);
    err = cudaFree(d_wave);
    err = cudaFree(d_kx);
    err = cudaFree(d_ky);
    err = cudaFree(d_k2);
    err = cudaFree(d_propx);
    err = cudaFree(d_propy);
    err = cudaFree(d_propxy);
    err = cudaFree(d_atomSupercell);
    err = cudaFree(d_projPotential);

    free(h_projPotential);

}


void sliceSupercell(float *atomSupercell, expPara *param, uint32_t *startPointList, uint32_t *numAtomList) {

    float zslice = 0.95*param->deltaZ;//set the limit slightly lower than than the actual boundary for float number comparison
    uint32_t i_slice = 0;

    printf("zslice is %f.\n", zslice);
    for (uint32_t i = 0;i < param->totalNumAtom; i++) {
        printf("current z is :%f\n", atomSupercell[i*6+3]);
        if (atomSupercell[i*6+3] > zslice) {
            startPointList[i_slice] = i;
            if (i_slice > 0)
                numAtomList[i_slice] = i - startPointList[i_slice - 1];
            else
                numAtomList[i_slice] = i;
            i_slice++;
            zslice += param->deltaZ;
        }

    }

}


void generateSplineCoeff(atomVZ_LUT_all *atomVZ_all) {
    float dlnr = log(MAX_ATOM_RADIUS/MIN_ATOM_RADIUS)/(NUM_RADIUS_SAMPLE - 1);
    for(uint32_t i=0; i < NUM_RADIUS_SAMPLE; i++) {
         atomVZ_all->spline_x[i] = MIN_ATOM_RADIUS * exp( i * dlnr );
    }

    printf("spline x sampling array created.\n");

    double **fparams = (double **) malloc2D( MAX_ELEMENT_Z + 1, NUM_FPARAM, sizeof(double));
    createFparams(fparams);
    printf("scattering factor table loaded.\n");

    for (uint32_t i = 0; i < MAX_ELEMENT_Z; i++) {

        uint32_t Z = i + 1;// atomic number equals to index plus one.
        /// Use the look-up-table to calculate the atomic potential
        for(uint32_t j = 0; j < NUM_RADIUS_SAMPLE; j++) {
            atomVZ_all->spline_y[i][j] = vzatom(Z, atomVZ_all->spline_x[j], fparams);
        }

        /// Fit the spline
        splinh(atomVZ_all->spline_x, atomVZ_all->spline_y[i], atomVZ_all->splineCoeff[i][0], atomVZ_all->splineCoeff[i][1], atomVZ_all->splineCoeff[i][2], NUM_RADIUS_SAMPLE);
    }

    free(fparams);


}


float vzatom(uint32_t Z, float radius, double **fparamms) {

    double suml, sumg, x, r;

    /* Lorenzian, Gaussian constants */
    const double al=300.8242834, ag=150.4121417;
    const double pi=3.141592654;

    r = fabs( radius );
    if( r < 1.0e-10 ) r = 1.0e-10;  /* avoid singularity at r=0 */
    suml = sumg = 0.0;

    /* Lorenztians */
    x = 2.0*pi*r;
    for(uint32_t i=0; i<2*NUM_GAUSS; i+=2 )
        suml += fparamms[Z][i]* bessk0( x*sqrt(fparamms[Z][i+1]) );

    /* Gaussians */
    x = pi*r;
    x = x*x;
    for(uint32_t i=2*NUM_GAUSS; i<2*(NUM_LORENZ+NUM_GAUSS); i+=2 )
        sumg += fparamms[Z][i]*exp(-x/fparamms[Z][i+1]) / fparamms[Z][i+1];

    return( al*suml + ag*sumg );

}

/*------------------ splinh() -----------------------------*/
/*
    fit a quasi-Hermite  cubic spline

    [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
        'A New Method of Interpolation and Smooth
        Curve Fitting Based on Local Procedures'

    [2] H.Akima, Comm. ACM, 15(1972)p.914-918

    E. Kirkland 4-JUL-85
    changed zero test to be a small nonzero number 8-jul-85 ejk
    converted to C 24-jun-1995 ejk

    The inputs are:
        x[n] = array of x values in ascending order, each X(I) must
            be unique
        y[n] = array of y values corresponding to X(N)
        n  = number of data points must be 2 or greater

    The outputs are (with z=x-x(i)):
        b[n] = array of spline coefficients for (x-x[i])
        c[n] = array of spline coefficients for (x-x[i])**2
        d[n] = array of spline coefficients for (x-x[i])**3
        ( x[i] <= x <= x[i+1] )
    To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

    The coefficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
    interval. NOTE that the last set of coefficients,
    b[n-1], c[n-1], d[n-1] are meaningless.
*/
void splinh( double *x, double *y, double *b, double *c, double *d, uint32_t n)
{


    uint32_t i, nm1, nm4;
    double m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

    if( n < 4) return;

    /* Do the first end point (special case),
       and get starting values */

    m5 = ( y[3] - y[2] ) / ( x[3] - x[2] ); /* mx = slope at pt x */
    m4 = ( y[2] - y[1] ) / ( x[2] - x[1] );
    m3 = ( y[1] - y[0] ) / ( x[1] - x[0] );

    m2 = m3 + m3 - m4;  /* eq. (9) of reference [1] */
    m1 = m2 + m2 - m3;

    m54 = fabs( m5 - m4);
    m43 = fabs( m4 - m3);
    m32 = fabs( m3 - m2);
    m21 = fabs( m2 - m1);

    if ( (m43+m21) > SMALL )
        t1 = ( m43*m2 + m21*m3 ) / ( m43 + m21 );
    else
        t1 = 0.5 * ( m2 + m3 );

    /*  Do everything up to the last end points */

    nm1 = n-1;
    nm4 = n-4;

    for( i=0; i<nm1; i++) {

        if( (m54+m32) > SMALL )
            t2= (m54*m3 + m32*m4) / (m54 + m32);
        else
            t2 = 0.5* ( m3 + m4 );

        x43 = x[i+1] - x[i];
        b[i] = t1;
        c[i] = ( 3.0*m3 - t1 - t1 - t2 ) /x43;
        d[i] = ( t1 + t2 - m3 - m3 ) / ( x43*x43 );

        m1 = m2;
        m2 = m3;
        m3 = m4;
        m4 = m5;
        if( i < nm4 ) {
            m5 = ( y[i+4] - y[i+3] ) / ( x[i+4] - x[i+3] );
        } else {
            m5 = m4 + m4 - m3;
        }

        m21 = m32;
        m32 = m43;
        m43 = m54;
        m54 = fabs( m5 - m4 );
        t1 = t2;
    }

    return;

} /* end splinh() */

/*-------------------- bessk0() ---------------*/
/*
    modified Bessel function K0(x)
    see Abramowitz and Stegun page 380

    Note: K0(0) is not define and this function
    returns 1E20

    x = (double) real arguments

    this routine calls bessi0() = Bessel function I0(x)

    12-feb-1997 E. Kirkland
 */
 double bessk0( double x )
 {
    double bessi0(double);

    int i;
    double ax, x2, sum;
    double k0a[] = { -0.57721566, 0.42278420, 0.23069756,
         0.03488590, 0.00262698, 0.00010750, 0.00000740};

    double k0b[] = { 1.25331414, -0.07832358, 0.02189568,
         -0.01062446, 0.00587872, -0.00251540, 0.00053208};

    ax = fabs( x );
    if( (ax > 0.0)  && ( ax <=  2.0 ) ){
        x2 = ax/2.0;
        x2 = x2 * x2;
        sum = k0a[6];
        for( i=5; i>=0; i--) sum = sum*x2 + k0a[i];
        sum = -log(ax/2.0) * bessi0(x) + sum;
    } else if( ax > 2.0 ) {
        x2 = 2.0/ax;
        sum = k0b[6];
        for( i=5; i>=0; i--) sum = sum*x2 + k0b[i];
        sum = exp( -ax ) * sum / sqrt( ax );
    } else sum = 1.0e20;
    return ( sum );

}  /* end bessk0() */


/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x )
 {
    int i;
    double ax, sum, t;

    double i0a[] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
        0.2659732, 0.0360768, 0.0045813 };

    double i0b[] = { 0.39894228, 0.01328592, 0.00225319,
        -0.00157565, 0.00916281, -0.02057706, 0.02635537,
        -0.01647633, 0.00392377};

    ax = fabs( x );
    if( ax <= 3.75 ) {
        t = x / 3.75;
        t = t * t;
        sum = i0a[6];
        for( i=5; i>=0; i--) sum = sum*t + i0a[i];
    } else {
        t = 3.75 / ax;
        sum = i0b[8];
        for( i=7; i>=0; i--) sum = sum*t + i0b[i];
        sum = exp( ax ) * sum / sqrt( ax );
    }
    return( sum );

}  /* end bessi0() */


/*----------------------- seval() ----------------------*/
/*
    Interpolate from cubic spline coefficients

    E. Kirkland 4-JUL-85
    modified to do a binary search for efficiency 13-Oct-1994 ejk
    converted to C 26-jun-1995 ejk
    fixed problem on end-of-range 16-July-1995 ejk

    The inputs are:
        x[n] = array of x values in ascending order, each x[i] must
            be unique
        y[n] = array of y values corresponding to x[n]
        b[n] = array of spline coefficients for (x-x[i])
        c[n] = array of spline coefficients for (x-x[i])**2
        d[n] = array of spline coefficients for (x-x[i])**3
        n  = number of data points
        x0  = the x value to interpolate at
        (x[i] <= x <= x[i+1]) and all inputs remain unchanged

    The value returned is the interpolated y value.

    The coefficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
    interval. NOTE that the last set of coefficients,
    b[n-1], c[n-1], d[n-1] are meaningless.
*/
__device__ double seval( double *x, double *y, double *b, double *c,
         double *d, int n, double x0 )
{
    int i, j, k;
    double z, seval1;

    /*  exit if x0 is outside the spline range */
    if( x0 <= x[0] ) i = 0;
    else if( x0 >= x[n-2] ) i = n-2;
    else {
        i = 0;
        j = n;
        do{ k = ( i + j ) / 2 ;
            if( x0 < x[k] )  j = k;
            else if( x0 >= x[k] ) i = k;
        } while ( (j-i) > 1 );
    }

    z = x0 - x[i];
    seval1 = y[i] + ( b[i] + ( c[i] + d[i] *z ) *z) *z;

    return( seval1 );

} /* end seval() */


__global__ void initialise(uint32_t nx, uint32_t ny, float2 *array, float reIni, float imIni) {
    uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y;

	array[(ix * ny)+iy].x = reIni; //real space coordinate
	array[(ix * ny)+iy].y = imIni; //real space coordinate

}

__global__ void createKArray_1D(uint32_t n, float dim, float *k) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t mid = n/2;

    if (i > mid) {
        k[i] = ((float)(i - n))/dim;
    }
    else {
        k[i] = ((float)i)/dim;
    }
}

__global__ void createKArray_2D(uint32_t nx, uint32_t ny, float *d_kx, float *d_ky, float *d_k2) {
    uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y;

	d_k2[(iy * nx) + ix] = d_kx[ix] * d_kx[ix] + d_ky[iy] * d_ky[iy]; //k space coordinate

}

__global__ void createProp_1D(uint32_t n, float *k, float scale, float wavelen, float2 *prop) {
    uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    float t;
    t = scale*k[i]*k[i]*wavelen;
    prop[i].x = acosf(t);
    prop[i].y = asinf(t);

}

__global__ void createProp_2D(uint32_t nx, float2 *d_propx, float2 *d_propy, float2 *d_propxy) {
    uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y;

    d_propxy[(iy * nx) + ix].x = d_propx[ix].x * d_propy[iy].x - d_propx[ix].y * d_propy[iy].y;
    d_propxy[(iy * nx) + ix].y = d_propx[ix].x * d_propy[iy].y + d_propx[ix].y * d_propy[iy].x;


}



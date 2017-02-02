#include "mslib.h"
#include "misc.h"
#include <cuda.h>
#include <cuda_profiler_api.h>

#define MAX_DSCRP_LENGTH 100
#define MAX_ATOM_RADIUS 3.0 //max atomic radius in the unit of Angstrom
#define MIN_ATOM_RADIUS 0.01 //min atomic radius in the unit of Angstrom
#define NUM_LORENZ 3
#define NUM_GAUSS 3
#define SMALL 1.0e-25
#define NUM_PARAMS 12 // number of parameters for calculating the scattering factor

void constructSupercell(FILE *fpCell, float *atomSupercell, msPara *para, sfInfo *info) {

    char currLine[MAX_DSCRP_LENGTH];
    uint32_t ncellx, ncelly, ncellz, nAtomsUnitCell;
    float cell_a, cell_b, cell_c;
    ///read the first description line
    fgets(currLine, MAX_DSCRP_LENGTH, fpCell);
    printf("%s", currLine);

    ///read the number of cells expanding in three directions
    fscanf(fpCell, "%u %u %u", &ncellx, &ncelly, &ncellz);
    printf("%u, %u, %u\n", ncellx, ncelly, ncellz);

    ///read the unit cell size
    fscanf(fpCell, "%f %f %f", &cell_a, &cell_b, &cell_c);
    printf("%f, %f, %f\n", cell_a, cell_b, cell_c);

    ///read the number of atoms in one unit cell
    fscanf(fpCell, "%u\n", &nAtomsUnitCell);
    float atomInfoOneCell[nAtomsUnitCell][6];
    printf("# of atoms in the unit cell is: %u\n", nAtomsUnitCell);

    ///allocate the memory to store the atom information in the whole lattice
    uint32_t totalNumAtoms = ncellx*ncelly*ncellz*nAtomsUnitCell;
    para->totalNumAtoms = totalNumAtoms;
    para->supercell_a = cell_a*ncellx;
    para->supercell_b = cell_b*ncelly;
    atomSupercell = (float *)malloc(totalNumAtoms*6*sizeof(float));

    ///read the atoms' infomation: atomic number, x, y, z, occupation and thermal displacement
    for (uint32_t readIndex = 0; readIndex < nAtomsUnitCell; readIndex ++) {
        fscanf(fpCell, "%f %f %f %f %f %f", atomInfoOneCell[readIndex], atomInfoOneCell[readIndex] + 1, atomInfoOneCell[readIndex] + 2, atomInfoOneCell[readIndex] + 3, atomInfoOneCell[readIndex] + 4, atomInfoOneCell[readIndex] + 5);
        printf("%2f %1.5f %1.5f %1.5f %1.1f %1.3f\n", atomInfoOneCell[readIndex][0], atomInfoOneCell[readIndex][1], atomInfoOneCell[readIndex][2], atomInfoOneCell[readIndex][3], atomInfoOneCell[readIndex][4], atomInfoOneCell[readIndex][5]);

        char currAtomZ = (char)atomInfoOneCell[readIndex][0];
        if ((currAtomZ > 0) && (currAtomZ <= MAX_ELEMENT_Z)) {
            if (info->isElement[currAtomZ - 1] == 0) {
                info->nAtomType++;
                info->isElement[currAtomZ - 1] = 1;
            }
        }
    }

    printf("initialize the supercell\n");
    memcpy(atomSupercell, atomInfoOneCell, nAtomsUnitCell*6*sizeof(float));
    for (uint32_t readIndex = 0; readIndex < nAtomsUnitCell; readIndex ++) {
        for (uint32_t i = 0; i < 6; i++){
            printf("%2f ", *(atomSupercell + readIndex*6 + i));
        }
        printf("\n");

    }

    ///replicate the unit cell on all three directions

    uint32_t currNAtoms = nAtomsUnitCell;
    if( ncellx > 1 ) {
        for(uint32_t i=1; i<ncellx; i++)
        for(uint32_t j=0; j<currNAtoms; j++) {
            *(atomSupercell + currNAtoms*i*6 + j + 0) = *(atomSupercell + j + 0);
            *(atomSupercell + currNAtoms*i*6 + j + 1) = *(atomSupercell + j + 1) + i*cell_a;
            *(atomSupercell + currNAtoms*i*6 + j + 2) = *(atomSupercell + j + 2);
            *(atomSupercell + currNAtoms*i*6 + j + 3) = *(atomSupercell + j + 3);
            *(atomSupercell + currNAtoms*i*6 + j + 4) = *(atomSupercell + j + 4);
            *(atomSupercell + currNAtoms*i*6 + j + 5) = *(atomSupercell + j + 5);

        }
        currNAtoms = currNAtoms*ncellx;
        //*ax = (*ax) * ncellx;
    }

    if( ncelly > 1 ) {
        for(uint32_t i=1; i<ncelly; i++)
        for(uint32_t j=0; j<currNAtoms; j++) {
            *(atomSupercell + currNAtoms*i*6 + j + 0) = *(atomSupercell + j + 0);
            *(atomSupercell + currNAtoms*i*6 + j + 1) = *(atomSupercell + j + 1);
            *(atomSupercell + currNAtoms*i*6 + j + 2) = *(atomSupercell + j + 2) + i*cell_b;
            *(atomSupercell + currNAtoms*i*6 + j + 3) = *(atomSupercell + j + 3);
            *(atomSupercell + currNAtoms*i*6 + j + 4) = *(atomSupercell + j + 4);
            *(atomSupercell + currNAtoms*i*6 + j + 5) = *(atomSupercell + j + 5);

        }
        currNAtoms = currNAtoms*ncelly;
        //*by = (*by) * ncelly;
    }

    if( ncellz > 1 ) {
        for(uint32_t i=1; i<ncellz; i++)
        for(uint32_t j=0; j<currNAtoms; j++) {
            *(atomSupercell + currNAtoms*i*6 + j + 0) = *(atomSupercell + j + 0);
            *(atomSupercell + currNAtoms*i*6 + j + 1) = *(atomSupercell + j + 1);
            *(atomSupercell + currNAtoms*i*6 + j + 2) = *(atomSupercell + j + 2);
            *(atomSupercell + currNAtoms*i*6 + j + 3) = *(atomSupercell + j + 3) + i*cell_c;
            *(atomSupercell + currNAtoms*i*6 + j + 4) = *(atomSupercell + j + 4);
            *(atomSupercell + currNAtoms*i*6 + j + 5) = *(atomSupercell + j + 5);
        }
        currNAtoms = currNAtoms*ncellz;
        //*cz = (*cz) * ncellz;
    }

}

__global__ void vzRealSpace(uint32_t numAtoms, sfInfo *info, float *atomSupercell, float *projPotential, uint32_t nx, uint32_t ny, float scale_x, float scale_y, uint32_t limit_x, uint32_t limit_y, float radMin_sq, float radMax_sq) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numAtoms) {
        uint32_t atomZ = (uint32_t)*(atomSupercell + i*6 + 0);
        float origin_x = *(atomSupercell + i*6 + 1);
        float origin_y = *(atomSupercell + i*6 + 2);
        float occupy = *(atomSupercell + i*6 + 4);

        double *spline_x = info->spline_x;
        double *spline_y = info->spline_y[atomZ];
        double *spline_b = info->splineCoeff[atomZ][0];
        double *spline_c = info->splineCoeff[atomZ][1];
        double *spline_d = info->splineCoeff[atomZ][2];

        int32_t start_x = (uint32_t)origin_x - limit_x;
        int32_t end_x = (uint32_t)origin_x + limit_x;
        int32_t start_y = (uint32_t)origin_y - limit_y;
        int32_t end_y = (uint32_t)origin_y + limit_y;

        float distance_x, distance_y, distance_x_sq, distance_y_sq, rsq, vz;
        int32_t ixw, iyw;

        for (int32_t ix = start_x; ix < end_x; ix++) {
            distance_x = origin_x - (float)ix*scale_x;
            ixw = ix;
            distance_x_sq = distance_x*distance_x;
            while( ixw < 0 ) ixw = ixw + nx;
            ixw = ixw % nx;
            for (int32_t iy = start_y; iy < end_y; iy++) {
                    distance_y = origin_y - (float)iy*scale_y;
                    distance_y_sq = distance_y*distance_y;
                    rsq = distance_x_sq + distance_y_sq;
                    if (rsq < radMax_sq) {
                        iyw = iy;
                        while( iyw < 0 ) iyw = iyw + ny;
                        iyw = iyw % ny;
                        if( rsq < radMin_sq ) rsq = radMin_sq;
                        vz = occupy * seval(spline_x, spline_y, spline_b, spline_c, spline_d, NUM_RADIUS_SAMPLE,(double)rsq);
                        atomicAdd(projPotential+ny*ixw + iyw, (float)vz);

                    }



            }

        }
    }

}

void calculatePhaseGrating(float *atomSupercell, msPara *para, sfInfo *info) {
    float scale_x = para->supercell_a/para->nx;
    float scale_y = para->supercell_b/para->ny;

    float radMin = 0.25 * sqrt(0.5*(scale_x*scale_x + scale_y*scale_y));
    float radMin_sq = radMin*radMin;
    float radMax_sq = MAX_ATOM_RADIUS*MAX_ATOM_RADIUS;

    uint32_t limit_x = (uint32_t)(MAX_ATOM_RADIUS/scale_x) + 1;
    uint32_t limit_y = (uint32_t)(MAX_ATOM_RADIUS/scale_y) + 1;

    // For now, assume only one slice
    // parallel the following part, use a compact spline interpolation coefficient table

    uint32_t numAtoms = para->totalNumAtoms;//number of atoms in one slice
    printf("total number of atoms is %d.\n",numAtoms);

    cudaError_t err = cudaSuccess;

    float *h_projPotential = (float *)malloc(para->nx*para->ny*sizeof(float));
    sfInfo *h_info = info;
    float *d_projPotential = NULL;

     for (uint32_t i = 0; i < para->nx*para->ny; i++) {
        h_projPotential[i] = 0;
    }

    sfInfo *d_info = NULL;
    float *d_atomSupercell = NULL;

    err = cudaMalloc((void **)&d_projPotential, para->nx*para->ny*sizeof(float));
    err = cudaMalloc((void **)&d_info, sizeof(sfInfo));
    err = cudaMalloc((void **)&d_atomSupercell, numAtoms*sizeof(float)*6);

    err = cudaMemcpy(d_info, h_info, sizeof(sfInfo), cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_projPotential, h_projPotential, para->nx*para->ny*sizeof(float), cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_atomSupercell, atomSupercell, numAtoms*sizeof(float)*6, cudaMemcpyHostToDevice);


    uint32_t threadsPerBlock = 256;
    uint32_t blocksPerGrid = (numAtoms + threadsPerBlock - 1)/threadsPerBlock;


    vzRealSpace<<<blocksPerGrid, threadsPerBlock>>>(numAtoms, d_info, d_atomSupercell, d_projPotential, para->nx, para->ny, scale_x, scale_y, limit_x, limit_y, radMin_sq, radMax_sq);


    err = cudaMemcpy(h_projPotential, d_projPotential, para->nx*para->ny*sizeof(float), cudaMemcpyDeviceToHost);

    err = cudaGetLastError();
    err = cudaFree(d_info);
    err = cudaFree(d_atomSupercell);
    err = cudaFree(d_projPotential);

    free(h_projPotential);


}


void generateSplineCoeff(sfInfo *info) {
    float dlnr = log(MAX_ATOM_RADIUS/MIN_ATOM_RADIUS)/(NUM_RADIUS_SAMPLE - 1);
    for(uint32_t i=0; i < NUM_RADIUS_SAMPLE; i++) {
         info->spline_x[i] = MIN_ATOM_RADIUS * exp( i * dlnr );
    }

    printf("spline X created");

    double **fparams = (double **) malloc2D( MAX_ELEMENT_Z + 1, NUM_PARAMS, sizeof(double));
    createFparams(fparams);
    printf("sf table loaded.");

    for (uint32_t i = 0; i < MAX_ELEMENT_Z; i++) {

        uint32_t Z = i + 1;// atomic number equals to index plus one.
        /// Use the look-up-table to calculate the atomic potential
        for(uint32_t j = 0; j < NUM_RADIUS_SAMPLE; j++) {
            info->spline_y[i][j] = vzatom(Z, info->spline_x[j], fparams);
        }

        /// Fit the spline
        splinh(info->spline_x, info->spline_y[i], info->splineCoeff[i][0], info->splineCoeff[i][1], info->splineCoeff[i][2], NUM_RADIUS_SAMPLE);
    }

    free(fparams);


}


float vzatom(uint32_t Z, float radius, double **fparams) {

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
        suml += fparams[Z][i]* bessk0( x*sqrt(fparams[Z][i+1]) );

    /* Gaussians */
    x = pi*r;
    x = x*x;
    for(uint32_t i=2*NUM_GAUSS; i<2*(NUM_LORENZ+NUM_GAUSS); i+=2 )
        sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / fparams[Z][i+1];

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

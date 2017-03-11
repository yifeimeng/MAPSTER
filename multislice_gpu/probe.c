
/*------------------------ chi() ---------------------*/
/*  return the aberration function
    must multiply this result by 2pi/wavelength before using
    C10= -df is defocus with the opposite sign

    put in a subroutine so I only have to get it right once

    now includes thru 5th order 4-jul-2010 ejk
    convert to param array p[] 2-may-2011 ejk
    add multiMode 8-may-2011 ejk
    fix bug in astigmatism calculation 4-jul-2012 ejk

input:
    p = param[] array with aberration coeff. of chi (in Angstroms)
    alx, aly = x,y components of alpha in radians
    multiMode = if not 0 then use multipole aberrations
           set to 0 to speed up the calculation (no C34a etc.)

*/


double chi(float p[], double alx, double aly, int multiMode) {
    double theta, al, al2, aln;
    double w, ct, st, c2t, s2t, c3t, s3t, c4t, s4t, c5t, s5t,
        c6t, s6t;

    aln = al2 = alx*alx + aly*aly;    /*  alpha squared */

    /*  just rotationally symm. aberrations (faster) */
    w = ( ( al2*p[pCS5]/6.0 +  0.25*p[pCS] )*al2 - 0.5*p[pDEFOCUS] )*al2;

    if( multiMode != 0 ) {
            /* ---- first order ----- */
            theta = atan2( aly, alx );
            ct = cos( theta );
            st = sin( theta );
            c2t = ct*ct - st*st;    /*  cos/sin of 2*theta */
            s2t = 2.0*ct*st;
            w += al2*(  p[pC12a]*c2t + p[pC12b]*s2t )/2.0;

            al = sqrt( al2 );

            /* ---- second order ----- */
            /*   generate the other theta's recursively to reduce CPU time */
            aln = al2*al;  /* alpha^3 */
            c3t = ct*c2t - st*s2t;    /*  cos/sin of 3*theta */
            s3t = ct*s2t + st*c2t;
            w += aln*( p[pC21a]*ct + p[pC21b]*st + p[pC23a]*c3t + p[pC23b]*s3t )/3.0;

            /* ---- third order ----- */
            aln = al2*al2;  /* alpha^4 */
            c4t = ct*c3t - st*s3t;    /*  cos/sin of 4*theta */
            s4t = ct*s3t + st*c3t;
            w += aln*( p[pC32a]*c2t + p[pC32b]*s2t + p[pC34a]*c4t
                     + p[pC34b]*s4t )/4.0;

             /* ---- fourth order ----- */
            aln = aln*al;  /* alpha^5 */
            c5t = ct*c4t - st*s4t;    /*  cos/sin of 5*theta */
            s5t = ct*s4t + st*c4t;
            w += aln*( p[pC41a]*ct + p[pC41b]*st +  p[pC43a]*c3t + p[pC43b]*s3t
                     + p[pC45a]*c5t + p[pC45b]*s5t )/5.0;

              /* ---- fifth order ----- */
            aln = aln*al;  /* alpha^6 */
            c6t = ct*c5t - st*s5t;    /*  cos/sin of 5*theta */
            s6t = ct*s5t + st*c5t;
            w += aln*(  p[pC52a]*c2t + p[pC52b]*s2t +  p[pC54a]*c4t + p[pC54b]*s4t
                      + p[pC56a]*c6t + p[pC56b]*s6t )/6.0;


    }

    return( w );
}


__global__ void calculateProbeWave(float2 *probeWave, float *pAberr, float k2max, uint8_t multiMode, float *d_kx, float *d_ky, float *d_k2, float xp, float yp) {

    uint32_t ix = blockIdx.x * blockDim.x + threadIdx.x;
	uint32_t iy = blockIdx.y * blockDim.y + threadIdx.y;
	uint32_t i = (iy * nx) + ix;

    float k2 = d_k2[i];
    float chi0;
	if (k2 <= k2max) {
        float aly = param->wavelen*d_ky[iy];
        float alx = param->wavelen*d_kx[ix];
        float dx2p = 2.0*PI*xp;
        float dy2p = 2.0*PI*yp;

        chi0 = (2.0*PI/d_param->wavelen)*chi(pAberr, alx, aly, multiMode) - ((dx2p*d_kx[ix])+(dy2p*d_ky[iy]));
        probeWave[i].x = _cosf(chi0);
        probeWave[i].y = _sinf(chi0);

	}
	else {
        probeWave[i].x = 0.0;
        probeWave[i].y = 0.0;
	}





}



void makeProbe(float2 *h_probe, expPara *param, float *p) {

    double xGH[]={ 3.190993201781528, 2.266580584531843, 1.468553289216668,
        0.723551018752838, 0.000000000000000, -0.723551018752838,
        -1.468553289216668,-2.266580584531843,-3.190993201781528};
    double wGH[]={3.960697726326e-005, 4.943624275537e-003 ,8.847452739438e-002,
        4.326515590026e-001, 7.202352156061e-001, 4.326515590026e-001,
        8.847452739438e-002, 4.943624275537e-003, 3.960697726326e-005};

    param->wavelen = wavelength(param->electronV); // calculate the wavelength
    sigmae = sigma( keV )/ 1000.0;

    float2 *d_probeWave;
    cudaMalloc((void **)&d_probeWave, param->nx*param->ny*sizeof(float2));

    float *d_kx, *d_ky, *d_k2;
    cudaMalloc((void **)&d_kx, param->nx*sizeof(float2));
    cudaMalloc((void **)&d_ky, param->ny*sizeof(float2));
    cudaMalloc((void **)&d_k2, param->nx*param->ny*sizeof(float2));
    createKArray_1D<<<dimGrid.x, dimBlock.x>>>(param->nx, param->supercell_a, d_kx);
    createKArray_1D<<<dimGrid.y, dimBlock.y>>>(param->ny, param->supercell_b, d_ky);
    createKArray_2D<<<dimGrid, dimBlock>>>(param->nx, param->ny, d_kx, d_ky, d_k2);

    float k2max = (p[pOAPERT]/param->wavelen)*(p[pOAPERT]/param->wavelen);

    float ddf2, df;
    uint32_t ndf;
    if (p[pDDF] > 1.0) {
        ndf = NGH;
        ddf2 = sqrt(log(2.0)/(p[pDDF]*p[pDDF]/4.0));
    }
    else {
        ndf = 1;
        ddf2 = 0;
    }

    for (uint32_t idf = 0; idf < ndf; idf ++) {
        if (ndf > 1) {

        df = param->df0 + xGH[idf]/ddf2;
        weight = (float) wGH[idf];
        }
        else {
            df = param->df0;
            weight = 1.0;
        }

        dim3 dimBlock( 16, 16);
        dim3 dimGrid( param->nx / dimBlock.x, param->ny / dimBlock.y);
        calculateProbeWave<<<dimGrid, dimBlock>>>(d_probeWave, p, k2max, multiMode, d_kx, d_ky, d_k2, xp, yp);

        cufftExecC2C(d_param->fft_plan, (cufftComplex *)d_probeWave, (cufftComplex *)d_probeWave, CUFFT_INVERSE);

        // calculate probe intensity at each pixel
        float *pixsq;
        module<<<dimGrid, dimBlock>>>(d_probeWave, d_pixsq);
        // normalize the probe
        float d_sum;
        float *ctf;
        reduceArray<<<1, param->nx*param->ny>>>(d_pixsq, &d_sum, param->nx*param->ny);
        scale = weight * ((float)sqrt( 1.0 / d_sum ) );

        scaleReal<<<dimGrid, dimBlock>>>(d_pixsq, d_pixsq, float scale);
        addReal<<<dimGrid, dimBlock>>>(d_ctf, d_pixsq, d_ctf);


    }

    /// convolve with the probe
    if (dsource > 0.001) {
    float2 *temp;
    real2Complex<<<dimGrid, dimBlock>>>(d_ctf, d_temp, param->nx*param->ny);
    cufftExecC2C(d_param->fft_plan, (cufftComplex *)d_temp, (cufftComplex *)d_temp, CUFFT_FORWARD);
    dsource = 0.5* dsource;   // convert diameter to radius
    ds = pi*pi * dsource*dsource/log(2.0);  // source size factor- convert to FWHM
    scaleComplex<<<dimGrid, dimBlock>>>(d_temp, d_temp, ds);
    cufftExecC2C(d_param->fft_plan, (cufftComplex *)d_temp, (cufftComplex *)d_temp, CUFFT_INVERSE);
    complex2real<<<dimGrid, dimBlock>>>(d_temp, d_ctf, param->nx*param->ny);

    }

    reduceArray<<<1, param->nx*param->ny>>>(d_ctf, &d_sum, param->nx*param->ny);
    scale = weight * ((float)sqrt( 1.0 / d_sum ) );
    scaleReal<<<dimGrid, dimBlock>>>(d_ctf, d_ctf, float scale);

    real2complex<<<dimGrid, dimBlock>>>(d_ctf, d_probe, param->nx*param->ny);

    cudaMemcpy(h_probe, d_probe, param->nx*param->ny*sizeof(float2), cudaMemcpyDeviceToHost);

}

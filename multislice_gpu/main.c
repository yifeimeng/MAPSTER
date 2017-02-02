#include <stdio.h>
#include <stdint.h>
#include "mslib.h"
#include "misc.h"

int main() {


    FILE *fpCell;
    fpCell = fopen("STO.xyz", "rb");
    if (fpCell == NULL) {
        printf("File open failed.");
        return(0);
    }

    float *atomSupercell;//Pointer to all atoms in the whole lattice.
    msPara *para = (msPara *)malloc(sizeof(msPara));
    sfInfo *info = (sfInfo *)malloc(sizeof(sfInfo));
    para->nx = 512;
    para->ny = 512;


    generateSplineCoeff(info);
    printf("generate spline coeff successful.\n");
    constructSupercell(fpCell, atomSupercell, para, info);
    printf("construct supercell successful.\n");

    calculatePhaseGrating(atomSupercell, para, info);
    printf("calculate phase grating successful.\n");



    fclose(fpCell);
    free(para);
    free(info);
    free(atomSupercell);
}


#include <stdio.h>
#include <stdint.h>
#include "mslib.h"
#include "misc.h"

#define MAX_STRING_LENGTH 100

int main() {

    ///get the LUT for atomic potential
    expPara *param = (expPara *)malloc(sizeof(expPara));
    atomVZ_LUT_all *atomVZ_all = (atomVZ_LUT_all *)malloc(sizeof(atomVZ_LUT_all));
    atomVZ_LUT *atomVZ = (atomVZ_LUT *)malloc(sizeof(atomVZ_LUT));
    generateSplineCoeff(atomVZ_all);
    printf("generate spline coeff successful.\n");

    ///initialize some parameters manually
    param->nx = 256;
    param->ny = 256;
    param->deltaZ = 1.95255;

    ///load the unit cell file
    FILE *fpCell;
    fpCell = fopen("STO.xyz", "rb");
    if (fpCell == NULL) {
        printf("File open failed.");
        return(0);
    }

    char currLine[MAX_STRING_LENGTH];
    uint32_t ncellx, ncelly, ncellz, numAtomUnitCell;
    float cell_a, cell_b, cell_c;
    ///read the first description line
    fgets(currLine, MAX_STRING_LENGTH, fpCell);
    printf("%s", currLine);

    ///read the number of cells expanding in three directions
    fscanf(fpCell, "%d %d %d", &ncellx, &ncelly, &ncellz);
    printf("%d, %d, %d\n", ncellx, ncelly, ncellz);
    param->ncell_x = ncellx;
    param->ncell_y = ncelly;
    param->ncell_z = ncellz;

    ///read the unit cell size
    fscanf(fpCell, "%f %f %f", &cell_a, &cell_b, &cell_c);
    printf("%f, %f, %f\n", cell_a, cell_b, cell_c);
    param->cell_a = cell_a;
    param->cell_b = cell_b;
    param->cell_c = cell_c;

    ///read the number of atoms in one unit cell
    fscanf(fpCell, "%u\n", &numAtomUnitCell);
    printf("# of atoms in the unit cell is: %u\n", numAtomUnitCell);
    param->numAtomUnitCell = numAtomUnitCell;
    float *atomUnitCell = (float *)malloc(numAtomUnitCell*6*sizeof(float));

    ///update some parameters
    param->totalNumAtom = ncellx*ncelly*ncellz*numAtomUnitCell;
    param->supercell_a = cell_a*(float)ncellx;
    param->supercell_b = cell_b*(float)ncelly;
    param->supercell_c = cell_c*(float)ncellz;

    ///construct the supercell (not sorted)
    float *atomSupercell = (float *)malloc(param->totalNumAtom*6*sizeof(float));
    constructSupercell(fpCell, atomSupercell, atomUnitCell, param, atomVZ_all, atomVZ);
    printf("construct supercell successful.\n");

    ///check the atomVZ look up table
    for (uint32_t i = 0; i < 100; i ++) {

        printf("%f ",atomVZ->spline_y[0][i]);

    }
    printf("\n");

    ///sort the supercell according to the z position
    float *atomSupercell_sorted = (float *)malloc(param->totalNumAtom*6*sizeof(float));
    memcpy(atomSupercell_sorted, atomSupercell, param->totalNumAtom*6*sizeof(float));
    heapsort(atomSupercell_sorted, param->totalNumAtom, 6, 4);// the z position is the 4th data in a block

    ///check the sorted supercell
    printf("check the sorted supercell.\n");
    for (uint32_t i = 0; i < param->totalNumAtom; i ++) {
        for (uint32_t j = 0; j < 6; j ++){
            printf("%f ", *(atomSupercell_sorted + i*6 + j));
        }
        printf("\n");

    }

    ///slice the sorted supercell
    displayParam(param);
    param->numSlice = (uint32_t)(param->supercell_c/param->deltaZ + 0.5);//round the float number to integer
    printf("the number of slice is %d.\n", param->numSlice);

    uint32_t *startPointList, *numAtomList;
    startPointList = (uint32_t *)malloc(param->numSlice*sizeof(uint32_t));
    numAtomList = (uint32_t *)malloc(param->numSlice*sizeof(uint32_t));
    sliceSupercell(atomSupercell_sorted, param, startPointList, numAtomList);
    printf("check slices:\n");
    for (uint32_t i = 0; i < param->numSlice; i ++) {
        printf("slice %d: start point %d, number of atom %d\n", i, startPointList[i], numAtomList[i]);
    }


    ///run the multislice, the phase grating is calculated on the fly
    multislice_run(atomSupercell_sorted, param, atomVZ, startPointList, numAtomList);

    fclose(fpCell);
    free(param);
    free(atomVZ);
    free(atomVZ_all);
    free(atomSupercell);
    free(atomUnitCell);
    free(startPointList);
    free(numAtomList);


}


#include "crystal.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_NUM_ATOM 80000
#define SMALL 0.001

int main() {
    
    FILE *fp, *fp_o, *fp_o2;
    float *unitCell;
    cell_param *cp, *cp_super;
    uint32_t numAtom;
    
    float rMat[9] = {0};
    
    rMat[0] = 0.70946004;
    rMat[1] = 0;
    rMat[2] = 0;
    rMat[3] = 0;
    rMat[4] = 1;
    rMat[5] = 0;
    rMat[6] = -0.70474566;
    rMat[7] = 0;
    rMat[8] = 1;
    
    
    unitCell = (float *)malloc(sizeof(float)*4*MAX_NUM_ATOM);
    cp = (cell_param *)malloc(sizeof(cell_param));
    fp = fopen("SrRuO3-62-[110]-supercell-113.cif", "r");
    if (fp == NULL) {
        printf("File not found!\n");
        return 0;
    }
    
    numAtom = readCIF(fp, unitCell, cp);
    printf("%d atoms read.\n", numAtom);
    
    float tVector[3] = {0};
    
    tVector[0] = 0;
    tVector[1] = 0;
    tVector[2] = cp->a * 0.70474566;
    
    cell2Angstrom(unitCell, numAtom, cp);
    rotateCell(unitCell, unitCell, numAtom, rMat);
    translateCell(unitCell, unitCell, numAtom, tVector);
    
    cp_super = (cell_param *)malloc(sizeof(cell_param));
    cp_super->a = cp->a * 0.70946004;
    cp_super->b = cp->b;
    cp_super->c = cp->c + cp->a * 0.70474566;
    cp_super->alpha = 90;
    cp_super->beta = 90;
    cp_super->gamma = 90;
    
    fp_o = fopen("SrRuO3-62-[110]-supercell-113.xyz", "w");
    printCellXYZ(fp_o, unitCell, numAtom);
    
    fp_o2 = fopen("SrRuO3-62-[110]-supercell-113-muSTEM.xtl", "w");
    angstrom2Cell(unitCell, numAtom, cp_super);
    printCellXTL_muSTEM(fp_o2, unitCell, numAtom, cp_super);
    

    
    fclose(fp);
    fclose(fp_o);
    free(unitCell);
    free(cp);
    free(cp_super);
    
    /*
    vector zone1, zone2;
    zone1.x = 0;
    zone1.y = 0;
    zone1.z = 1;

    zone2.x = 1;
    zone2.y = 1;
    zone2.z = 0;

    float rMat[9] = {0};

    //rotMatrix(zone1, zone2, rMat);
    
    rMat[0] = -0.704745;
    rMat[1] = 0.709460;
    rMat[2] = 0;
    rMat[3] = 0;
    rMat[4] = 0;
    rMat[5] = 1;
    rMat[6] = 0.709460;
    rMat[7] = 0.704745;
    rMat[8] = 0;
    
    printf("the rotation matrix is:\n");
    printf("%f %f %f\n", rMat[0], rMat[1], rMat[2]);
    printf("%f %f %f\n", rMat[3], rMat[4], rMat[5]);
    printf("%f %f %f\n", rMat[6], rMat[7], rMat[8]);
    
    FILE *fp, *fp_o, *fp_o2;
    float *unitCell, *superCell, *newCell, *croppedCell;
    cell_param *cp;
    
    uint32_t ncell_x = 1, ncell_y = 1, ncell_z = 1;
    
    cp = (cell_param *)malloc(sizeof(cell_param));
    fp = fopen("SrRuO3-62-[001]-supercell.xtl", "r");
    if (fp == NULL) {
        printf("File not found!\n");
        return 0;
    }
    
    unitCell = (float *)malloc(sizeof(float)*4*MAX_NUM_ATOM);
    superCell = (float *)malloc(sizeof(float)*4*ncell_x*ncell_y*ncell_z*MAX_NUM_ATOM);
    newCell = (float *)malloc(sizeof(float)*4*ncell_x*ncell_y*ncell_z*MAX_NUM_ATOM);
    croppedCell = (float *)malloc(sizeof(float)*4*ncell_x*ncell_y*ncell_z*MAX_NUM_ATOM);
    
    uint32_t numAtom;
    numAtom = readXTL_VESTA(fp, unitCell, cp);
    
    printf("num of atom is %d\n", numAtom);
    
   
    replicateCell(unitCell, cp, numAtom, superCell, ncell_x, ncell_y, ncell_z);
    
    numAtom = numAtom*ncell_x*ncell_y*ncell_z;
    
    cell2Angstrom(superCell, numAtom, cp);
    
    rotateCell(superCell, newCell, numAtom, rMat);
    
    float crop_param[6];
    
    crop_param[0] = 25; // x upper limit
    crop_param[1] = -25;    // x lower limit
    crop_param[2] = 50;    // y upper limit
    crop_param[3] = 0;    // y lower limit
    crop_param[4] = 250;    // z upper limit
    crop_param[5] = 0;    // z lower limit
    
   
    
    cell_param *cp_super;
    cp_super = (cell_param *)malloc(sizeof(cell_param));
    cp_super->a = 50;
    cp_super->b = 50;
    cp_super->c = 250;
    cp_super->alpha = 90;
    cp_super->beta = 90;
    cp_super->gamma = 90;
    
    fp_o = fopen("SrRuO3-62-[110]-supercell.xyz", "w");
    printCellXYZ(fp_o, newCell, numAtom);
    
    fp_o2 = fopen("SrRuO3-62-[110]-supercell-muSTEM.xtl", "w");
    numAtom = cropCell(newCell, croppedCell, numAtom, crop_param);
    angstrom2Cell(croppedCell, numAtom, cp_super);
    printCellXTL_muSTEM(fp_o2, croppedCell, numAtom, cp_super); // muSTEM only take fractional coordinates
    
    
    
    
    fclose(fp);
    fclose(fp_o);
    fclose(fp_o2);
    free(cp);
    free(cp_super);
    free(unitCell);
    free(newCell);
    free(superCell);
    free(croppedCell);
     
     */
    
    return 1;
    
    
}

#ifndef crystal
#define crystal

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

struct element {
    char name[3];
    uint8_t Z;

};

struct vector {
    float x;
    float y;
    float z;

};

struct cell_param {
    float a;
    float b;
    float c;
    float alpha;
    float beta;
    float gamma;
    
};

typedef struct vector vector;
typedef struct cell_param cell_param;

uint32_t cropCell(float *cell, float *newCell, uint32_t numAtom, float *crop_param);
void angstrom2Cell(float *cell, uint32_t numAtom, cell_param *cp);
void cell2Angstrom(float *cell, uint32_t numAtom, cell_param *cp);
void replicateCell(float *unitCell, cell_param *cp, uint32_t numAtom, float *newCell, uint32_t ncell_x, uint32_t ncell_y, uint32_t ncell_z);
uint32_t name2Z(char *name);
void Z2name(uint32_t Z, char *name);
vector rotateVector(vector u, float *R);
void rotateCell(float *oldCell, float *newCell, uint32_t numAtom, float *R);
void translateCell(float *oldCell, float *newCell, uint32_t numAtom, float *T);
float magnitude(vector u);
vector crossProduct(vector u, vector v);
float dotProduct(vector u, vector v);
float angleTwoVector(vector u, vector v);
void rotMatrix(vector u, vector v, float *R);
uint32_t readXTL_VESTA(FILE *fp, float *unitCell, cell_param *p);
void printCellXYZ(FILE *fp, float *cell, uint32_t numAtom);
void printCellXTL_muSTEM(FILE *fp, float *cell, uint32_t numAtom, cell_param *cp);
uint32_t readCIF(FILE *fp, float *unitCell, cell_param *cp);

#endif // crystal

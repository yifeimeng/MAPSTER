#include "crystal.h"
#define MAX_STRING_LENGTH 200
#define NUM_ELEMENT 3


uint32_t name2Z(char *name) {
    
    if (strcmp(name, "Sr") == 0) return 38;
    if (strcmp(name, "Ru") == 0) return 44;
    if (strcmp(name, "O") == 0) return 8;
    if (strcmp(name, "Ti") == 0) return 22;
    
    return 0;
}

void Z2name(uint32_t Z, char *name) {
    
    switch (Z) {
        case 38:
            strcpy(name, "Sr");
            break;
        case 44:
            strcpy(name, "Ru");
            break;
        case 8:
            strcpy(name, "O");
            break;
        case 22:
            strcpy(name, "Ti");
            break;
        default:
            strcpy(name, "NA");
    }
}

vector rotateVector(vector u, float *R) {
    vector v;
    
    v.x = R[0]*u.x + R[1]*u.y + R[2]*u.z;
    v.y = R[3]*u.x + R[4]*u.y + R[5]*u.z;
    v.z = R[6]*u.x + R[7]*u.y + R[8]*u.z;
    
    return v;
    
}

void rotateCell(float *oldCell, float *newCell, uint32_t numAtom, float *R) {
    
    vector newV, oldV;
    
    for (uint32_t i = 0; i < numAtom; i++) {
        oldV.x = oldCell[i*4 + 1];
        oldV.y = oldCell[i*4 + 2];
        oldV.z = oldCell[i*4 + 3];
        
        newV = rotateVector(oldV, R);
        
        newCell[i*4 + 0] = oldCell[i*4 + 0];
        newCell[i*4 + 1] = newV.x;
        newCell[i*4 + 2] = newV.y;
        newCell[i*4 + 3] = newV.z;
    }
    
    printf("cell rotation successful!\n");
    
}

void translateCell(float *oldCell, float *newCell, uint32_t numAtom, float *T) {
    
    for (uint32_t i = 0; i < numAtom; i++) {
        
        newCell[i*4 + 0] = oldCell[i*4 + 0];
        newCell[i*4 + 1] = oldCell[i*4 + 1] + T[0];
        newCell[i*4 + 2] = oldCell[i*4 + 2] + T[1];
        newCell[i*4 + 3] = oldCell[i*4 + 3] + T[2];
    }
    
    printf("cell translation successful!\n");
    
}

uint32_t cropCell(float *cell, float *newCell, uint32_t numAtom, float *crop_param) {

    float x, y, z, xs, ys, zs;
    
    float xlim_u, xlim_l, ylim_u, ylim_l, zlim_u, zlim_l;
    float xshift, yshift, zshift;
    
    xlim_u = crop_param[0];
    xlim_l = crop_param[1];
    ylim_u = crop_param[2];
    ylim_l = crop_param[3];
    zlim_u = crop_param[4];
    zlim_l = crop_param[5];
    
    xshift = -1 * xlim_l;
    yshift = -1 * ylim_l;
    zshift = -1 * zlim_l;
    
    
    uint32_t count = 0;
    
    /// crop along the 111 direction. this function is related to the rotation matrix
    for (uint32_t i = 0; i < numAtom; i ++) {
        
        x = cell[i*4 + 1];
        y = cell[i*4 + 2];
        z = cell[i*4 + 3];
        
        if ((x < xlim_u) && (x > xlim_l) && (y < ylim_u) && (y > ylim_l) && (z < zlim_u) && (z > zlim_l)) {
            
            printf("point %f %f %f include.\n", x, y, z);
            
            xs = x + xshift;
            ys = y + yshift;
            zs = z + zshift;
            
            printf("%f %f %f\n", xs, ys, zs);
            
            newCell[count*4 + 0] = cell[i*4 + 0];
            newCell[count*4 + 1] = xs;
            newCell[count*4 + 2] = ys;
            newCell[count*4 + 3] = zs;
            
            count ++;
        
        }
    
    }
    
    printf("%d atoms left after cropping.\n",count);
    return count;
    

}

void angstrom2Cell(float *cell, uint32_t numAtom, cell_param *cp) {

    for (uint32_t i = 0; i < numAtom; i++) {
        
        cell[i * 4 + 1] = cell[i * 4 + 1] / cp->a;
        cell[i * 4 + 2] = cell[i * 4 + 2] / cp->b;
        cell[i * 4 + 3] = cell[i * 4 + 3] / cp->c;
        
    }

    printf("Angstrom coordinations converted to fractions.\n");

}

void cell2Angstrom(float *cell, uint32_t numAtom, cell_param *cp) {
    
    for (uint32_t i = 0; i < numAtom; i ++) {
  
        //printf("%d\n", Z);
        //printf("%s\n", element);
        cell[i * 4 + 1] = cell[i * 4 + 1] * cp->a;
        cell[i * 4 + 2] = cell[i * 4 + 2] * cp->b;
        cell[i * 4 + 3] = cell[i * 4 + 3] * cp->c;
        
    }
    
    printf("Fractional coordinations converted to angstrom.\n");
    
}

void printCellXTL_muSTEM(FILE *fp, float *cell, uint32_t numAtom, cell_param *cp) {
    
    float *element[NUM_ELEMENT];
    
    for (uint32_t i = 0; i < NUM_ELEMENT; i ++) {
        element[i] = (float *)malloc(sizeof(float)*4*numAtom);
    
    }
    
    uint32_t elementList[NUM_ELEMENT] = {0};
    uint32_t elementCount[NUM_ELEMENT] = {0};
    uint32_t elementIndex = 0;
    
    for (uint32_t i = 0; i < numAtom; i ++) {
        
        char found = 0;
        uint32_t Z = (uint32_t)(cell[i*4] + 0.5);
        
        for (uint32_t j = 0; j < NUM_ELEMENT; j ++) {
            
            if (Z == elementList[j]) {
                found = 1;
               
                *(element[j] + elementCount[j]*4 + 0) = cell[i*4 + 0];
                *(element[j] + elementCount[j]*4 + 1) = cell[i*4 + 1];
                *(element[j] + elementCount[j]*4 + 2) = cell[i*4 + 2];
                *(element[j] + elementCount[j]*4 + 3) = cell[i*4 + 3];
                elementCount[j] ++;
                
            }
    
        }
        
        if (found == 0) {
            
            *(element[elementIndex] + 0) = cell[i*4 + 0];
            *(element[elementIndex] + 1) = cell[i*4 + 1];
            *(element[elementIndex] + 2) = cell[i*4 + 2];
            *(element[elementIndex] + 3) = cell[i*4 + 3];
            
            elementList[elementIndex] = Z;
            elementCount[elementIndex] ++;
            elementIndex ++;
        }
        
    }
    
    char *elementName = (char *)malloc(sizeof(char)*3);

    printf("sort atom list successful.\n");
    
    for (uint32_t i = 0; i < NUM_ELEMENT; i ++) {
        
        uint32_t Z = elementList[i];
        Z2name(Z, elementName);
        printf("%d %s (%d) atoms\n", elementCount[i], elementName, elementList[i]);
        
    }
    
    float occupy = 1.00;
    float B = 0.00634;
    
    
    fprintf(fp, "new structure\n");
    fprintf(fp, "%f %f %f %f %f %f\n", cp->a, cp->b, cp->c, cp->alpha, cp->beta, cp->gamma);
    fprintf(fp, "%f\n", 300.00);
    fprintf(fp, "%d\n", NUM_ELEMENT);
    
    
    for (uint32_t i = 0; i < NUM_ELEMENT; i ++) {
        
        uint32_t Z = elementList[i];
        Z2name(Z, elementName);
        fprintf(fp, "%s\n", elementName);
        fprintf(fp, "%d %d %f %f\n", elementCount[i], Z, occupy, B);
        
        
        for (uint32_t j = 0; j < elementCount[i]; j ++) {
            
            fprintf(fp, "%f %f %f\n", (*(element[i] + j*4 + 1)), (*(element[i] + j*4 + 2)), (*(element[i] + j*4 + 3)));
            
        }
        
        
    }
    
    for (uint32_t i = 0; i < NUM_ELEMENT; i ++) {
        free(element[i]);
        
    }
    
    free(elementName);
    
    fprintf(fp, "Orientation\n");
    fprintf(fp, "0 0 1\n");
    fprintf(fp, "1 0 0\n");
    fprintf(fp, "0 1 0\n");

}

void printCellXYZ(FILE *fp, float *cell, uint32_t numAtom) {
    
    char *element = (char *)malloc(sizeof(char)*3);
    uint32_t Z;
    
    fprintf(fp, "%d\n", numAtom);
    fprintf(fp, "%s\n", "new structure");
    
    for (uint32_t i = 0; i < numAtom; i ++) {
        
        Z = (uint32_t)(cell[i * 4] + 0.5);
        Z2name(Z, element);
        //printf("%d\n", Z);
        //printf("%s\n", element);
        fprintf(fp, "%s %f %f %f\n", element, cell[i * 4 + 1], cell[i * 4 + 2], cell[i * 4 + 3]);
        
    }
    
    free(element);
    
}

void replicateCell(float *unitCell, cell_param *cp, uint32_t numAtom, float *newCell, uint32_t ncell_x, uint32_t ncell_y, uint32_t ncell_z) {
    
    
    memcpy(newCell, unitCell, numAtom*4*sizeof(float));
    
    ///replicate the unit cell on all three directions
    
    uint32_t currNAtoms = numAtom;
    
    if( ncell_x > 1 ) {
        for (uint32_t i=1; i<ncell_x; i++) {
            for (uint32_t j=0; j<currNAtoms; j++) {
                *(newCell + currNAtoms*i*4 + j*4 + 0) = *(newCell + j*4 + 0);
                *(newCell + currNAtoms*i*4 + j*4 + 1) = *(newCell + j*4 + 1) + i;
                *(newCell + currNAtoms*i*4 + j*4 + 2) = *(newCell + j*4 + 2);
                *(newCell + currNAtoms*i*4 + j*4 + 3) = *(newCell + j*4 + 3);
                
            }
        }
        currNAtoms = currNAtoms*ncell_x;
        //*ax = (*ax) * ncellx;
    }
    
    if( ncell_y > 1 ) {
        for(uint32_t i=1; i<ncell_y; i++)
            for(uint32_t j=0; j<currNAtoms; j++) {
                *(newCell + currNAtoms*i*4 + j*4 + 0) = *(newCell + j*4 + 0);
                *(newCell + currNAtoms*i*4 + j*4 + 1) = *(newCell + j*4 + 1);
                *(newCell + currNAtoms*i*4 + j*4 + 2) = *(newCell + j*4 + 2) + i;
                *(newCell + currNAtoms*i*4 + j*4 + 3) = *(newCell + j*4 + 3);
                
            }
        currNAtoms = currNAtoms*ncell_y;
        //*by = (*by) * ncelly;
    }
    
    if( ncell_z > 1 ) {
        for(uint32_t i=1; i<ncell_z; i++)
            for(uint32_t j=0; j<currNAtoms; j++) {
                *(newCell + currNAtoms*i*4 + j*4 + 0) = *(newCell + j*4 + 0);
                *(newCell + currNAtoms*i*4 + j*4 + 1) = *(newCell + j*4 + 1);
                *(newCell + currNAtoms*i*4 + j*4 + 2) = *(newCell + j*4 + 2);
                *(newCell + currNAtoms*i*4 + j*4 + 3) = *(newCell + j*4 + 3) + i;
                
            }
        currNAtoms = currNAtoms*ncell_z;
        //*cz = (*cz) * ncellz;
    }
    
    
    printf("replicate cell successful!\n");
    
    
    
    
}

uint32_t readCIF(FILE *fp, float *unitCell, cell_param *cp) {
    
    char currLine[MAX_STRING_LENGTH];
    char str[MAX_STRING_LENGTH], currElement[3];
    char iniChar, *r;
    float value;
    uint32_t count = 0, numLoop = 0;
    
    r = fgets(currLine, MAX_STRING_LENGTH, fp);
    
    while (r != NULL) {
        
        iniChar = currLine[0];
        if ((iniChar != '#') && (iniChar != '\n')) {
            if (iniChar == '_') {
                sscanf(currLine, "%s %f", str, &value);
    
                if (strcmp("cell_length_a", str + 1) == 0) {
                    cp->a = value;
                }
                
                if (strcmp("cell_length_b", str + 1) == 0) {
                    cp->b = value;
                }
                
                if (strcmp("cell_length_c", str + 1) == 0) {
                    cp->c = value;
                }
                
            }
            
            if (currLine[0] == 'l' && currLine[1] == 'o' && currLine[2] == 'o' && currLine[3] == 'p' && currLine[4] == '_') {
                
                numLoop ++;
                
                
                if (numLoop == 2) {
                    float occupy, x, y, z, B;
                    uint32_t Z;
                    
                    r = fgets(currLine, MAX_STRING_LENGTH, fp);
                    
                 
                    while (currLine[0] == '_') {
                        r= fgets(currLine, MAX_STRING_LENGTH, fp);
                        
                    }
                
                    while ((r != NULL) && (currLine[0] != '\n')) {
                    
                        sscanf(currLine, " %s %s %f %f %f %f %f", str, currElement, &occupy, &x, &y, &z, &B);
                        
                        printf("%s %s %f %f %f %f %f\n", str, currElement, occupy, x, y, z, B);
                        
                        Z = name2Z(currElement);
                        unitCell[count*4 + 0] = (float)Z;
                        unitCell[count*4 + 1] = x ;
                        unitCell[count*4 + 2] = y ;
                        unitCell[count*4 + 3] = z ;

                        count ++;
                        r = fgets(currLine, MAX_STRING_LENGTH, fp);
                    }
                    
                }
                
               
                
            }
            
            
        }
        
        r= fgets(currLine, MAX_STRING_LENGTH, fp);
        
    }
    
    return count;
    
}

uint32_t readXTL_VESTA(FILE *fp, float *unitCell, cell_param *cp) {
    char currLine[MAX_STRING_LENGTH];
    fgets(currLine, MAX_STRING_LENGTH, fp);
    fgets(currLine, MAX_STRING_LENGTH, fp);
    fscanf(fp, "%f %f %f %f %f %f\n", &(cp->a), &(cp->b), &(cp->c), &(cp->alpha), &(cp->beta), &(cp->gamma));
    fgets(currLine, MAX_STRING_LENGTH, fp);
    fgets(currLine, MAX_STRING_LENGTH, fp);
    fgets(currLine, MAX_STRING_LENGTH, fp);
    fgets(currLine, MAX_STRING_LENGTH, fp);
    
    printf("%f %f %f %f %f %f\n", cp->a, cp->b, cp->c, cp->alpha, cp->beta, cp->gamma);
    
    char currElement[3];
    currElement[2] = '\0';
    float x, y, z;
    uint32_t count = 0, Z;
    
    while (strcmp(fgets(currLine, MAX_STRING_LENGTH, fp), "EOF") != 0) {
        
     
        sscanf(currLine, "%s %f %f %f\n", currElement, &x, &y, &z);
        printf("%s %f %f %f\n", currElement, x, y, z);
        Z = name2Z(currElement);
        //printf("%d\n", Z);
        unitCell[count*4 + 0] = (float)Z;
        unitCell[count*4 + 1] = x ;
        unitCell[count*4 + 2] = y ;
        unitCell[count*4 + 3] = z ;
        
        count++;
    }
    
    return count;

}


/* calculate the cross product of two 3d vectors */
float magnitude(vector u) {
    float result;
    result = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    return result;

}

vector crossProduct(vector u, vector v) {
    vector cp;
    cp.x = u.y * v.z - u.z * v.y;
    cp.y = u.z * v.x - u.x * v.z;
    cp.z = u.x * v.y - u.y * v.x;

    return cp;

}

float dotProduct(vector u, vector v) {
    float dp;
    dp = u.x * v.x + u.y * v.y + u.z*v.z;
    return dp;
}


float angleTwoVector(vector u, vector v) {
    float dp_normalized, angle;
    dp_normalized = dotProduct(u ,v)/(magnitude(u) * magnitude(v));
    angle = acos(dp_normalized); // the angle unit is rad
    return angle;

}

void rotMatrix(vector u, vector v, float *R) {
    float angle = angleTwoVector(u ,v);
    vector rotAxis = crossProduct(u, v);
    R[0] = cos(angle) + rotAxis.x * rotAxis.x * (1 - cos(angle));
    R[1] = rotAxis.x * rotAxis.y * (1 - cos(angle)) - rotAxis.z * sin(angle);
    R[2] = rotAxis.x * rotAxis.z * (1 - cos(angle)) + rotAxis.y * sin(angle);
    R[3] = rotAxis.y * rotAxis.x  * (1 - cos(angle)) + rotAxis.z * sin(angle);
    R[4] = cos(angle) + rotAxis.y * rotAxis.y * (1 - cos(angle));
    R[5] = rotAxis.y * rotAxis.z * (1 - cos(angle)) - rotAxis.x * sin(angle);
    R[6] = rotAxis.z * rotAxis.x * (1 - cos(angle)) - rotAxis.y * sin(angle);
    R[7] = rotAxis.z * rotAxis.y * (1 - cos(angle)) + rotAxis.x * sin(angle);
    R[8] = cos(angle) + rotAxis.z * rotAxis.z * (1 - cos(angle));

}                                                                                                 

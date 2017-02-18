#include "crystal.h"

void main() {
    vector3D zone1, zone2;
    zone1.x = 0;
    zone1.y = 0;
    zone1.z = 1;

    zone2.x = 1;
    zone2.y = 1;
    zone2.z = 0;

    float rMat[9] = {0};

    rotMatrix(zone1, zone2, rMat);

    printf("the rotation matrix is:\n");
    printf("%f %f %f\n", rMat[0], rMat[1], rMat[2]);
    printf("%f %f %f\n", rMat[3], rMat[4], rMat[5]);
    printf("%f %f %f\n", rMat[6], rMat[7], rMat[8]);

}

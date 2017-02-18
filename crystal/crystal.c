#include "crystal.h"
#include <math.h>

/* calculate the cross product of two 3d vectors */
float magnitude(vector3D u) {
    float result;
    result = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    return result;

}

vector3D crossProduct(vector3D u, vector3D v) {
    vector3D cp;
    cp.x = u.y * v.z - u.z * v.y;
    cp.y = u.z * v.x - u.x * v.z;
    cp.z = u.x * v.y - u.y * v.x;

    return cp;

}

float dotProduct(vector3D u, vector3D v) {
    float dp;
    dp = u.x * v.x + u.y * v.y + u.z*v.z;
    return dp;
}


float angleTwoVector(vector3D u, vector3D v) {
    float dp_normalized, angle;
    dp_normalized = dotProduct(u ,v)/(magnitude(u) * magnitude(v));
    angle = acos(dp_normalized); // the angle unit is rad
    return angle;

}

void rotMatrix(vector3D u, vector3D v, float *R) {
    float angle = angleTwoVector(u ,v);
    vector3D rotAxis = dotProduct(u, v);
    R[0] = cos(angle) + rotAxis.x * rotAxis.x * (1 - cos(angle));
    R[1] = rotAxis.x * rotAxis.y * (1 - cos(angle)) - rotAxis.z * sin(angle);
    R[2] = rotAxis.x * rotAxis.z * (1 - cos(angle)) + rotAxis.y * sin(angle);
    R[3] = rotAxis.y * rotAxis.z  * (1 - cos(angle)) + rotAxis.z * sin(angle);
    R[4] = cos(angle) + rotAxis.y * rotAxis.y * (1 - cos(angle));
    R[5] = rotAxis.y * rotAxis.z * (1 - cos(angle)) - rotAxis.x * sin(angle);
    R[6] = rotAxis.z * rotAxis.x * (1 - cos(angle)) - rotAxis.y * sin(angle);
    R[7] = rotAxis.z * rotAxis.y * (1 - cos(angle)) + rotAxis.x * sin(angle);
    R[8] = cos(angle) + rotAxis.z * rotAxis.z * (1 - cos(angle));

}                                                                                                 

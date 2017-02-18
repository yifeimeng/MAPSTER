#ifndef crystal
#define crystal

struct vector3D {
    float x;
    float y;
    float z;

};

float magnitude(vector3D u);
vector3D crossProduct(vector3D u, vector3D v);
float dotProduct(vector3D u, vector3D v);
float angleTwoVector(vector3D u, vector3D v);
void rotMatrix(vector3D u, vector3D v, float *R);

#endif // crystal

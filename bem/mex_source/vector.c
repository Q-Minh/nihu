#include "vector.h"

/* Compute dot product of two 3D vectors */
double dot(const double *v1,
           const double *v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/* Compute dot product of two 2D vectors */
double dot2D(const double *v1,
           const double *v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1];
}

/* Compute cross product of two vectors */
double * cross(const double *v1,
               const double *v2,
               double *w)
{
    w[0] = v1[1]*v2[2]-v1[2]*v2[1];
    w[1] = v1[2]*v2[0]-v1[0]*v2[2];
    w[2] = v1[0]*v2[1]-v1[1]*v2[0];

    return w;
}

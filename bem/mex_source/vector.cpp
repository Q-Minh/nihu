/**
* \file vector.cpp
* \brief Basic vector arithmetics
* \author Peter Fiala fiala@hit.bme.hu
*/

#include "vector.h"

/**
* \brief dot product of 3D vectors
* \param v1 left vector
* \param v2 right vector
* \returns dot product
*/
double dot(double const v1[],
           double const v2[])
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/**
* \brief dot product of 2D vectors
* \param v1 left vector
* \param v2 right vector
* \returns dot product
*/
double dot2D(double const v1[],
           double const v2[])
{
    return v1[0]*v2[0] + v1[1]*v2[1];
}

/**
* \brief cross product of 3D vectors
* \param v1 left vector
* \param v2 right vector
* \param w cross product
*/
double *cross(double const v1[],
              double const v2[],
              double w[])
{
    w[0] = v1[1]*v2[2]-v1[2]*v2[1];
    w[1] = v1[2]*v2[0]-v1[0]*v2[2];
    w[2] = v1[0]*v2[1]-v1[1]*v2[0];

    return w;
}


/**
* \file vector.h
* \brief header file of vector.cpp
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef VECTOR_H
#define VECTOR_H

double dot(double const v1[],
           double const v2[]);
		   
double dot2D(double const v1[],
             double const v2[]);
		   
double *cross(double const v1[],
              double const v2[],
              double w[]);

#endif


#ifndef VECTOR_H
#define VECTOR_H

/* Compute dot product of two 3D vectors */
double dot(const double *,
           const double *);
		   
/* Compute dot product of two 2D vectors */
double dot2D(const double *,
           const double *);
		   
/* Compute cross product of two 3D vectors */
double * cross(const double *,
               const double *,
               double *);
#endif


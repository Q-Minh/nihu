#ifndef MESH_H
#define MESH_H

#include "types.h"

/* ------------------------------------------------------------------------ */
/* Compute the center node of each element */
void init_accelerators(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       accelerator_t *accelerators);

/* ------------------------------------------------------------------------ */
/* Compute the center node of each element */
void init_accelerators2D(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       accelerator2D_t *accelerators);

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division(const double *q,
                   const double *elemcenter,
                   const double *dist);

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division2D(const double *q,
                   const double *elemcenter,
                   const double *dist);


#endif

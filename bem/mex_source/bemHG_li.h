#ifndef BEMHG_LI_H
#define BEMHG_LI_H
#include "types.h"

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
                     double *Ar,
                     double *Ai,
                     double *Br,
                     double *Bi);

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_lin2D(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss2D_t *g,
                     const double *dist,
                     double k,
                     double *Ar,
                     double *Ai,
                     double *Br,
                     double *Bi);

/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
void matrix_field_lin(int nnodes,
                      const double *nodes,
                      int nelements,
                      const double *elements,
                      int npoints,
                      const double *points,
                      const gauss_t *g3,
                      const gauss_t *g4,
                      const double *dist,
                      double k,
                      double *Ar,
                      double *Ai,
                      double *Br,
                      double *Bi);

/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
void matrix_field_lin2D(int nnodes,
                      const double *nodes,
                      int nelements,
                      const double *elements,
                      int npoints,
                      const double *points,
                      const gauss_t *g,
                      const double *dist,
                      double k,
                      double *Ar,
                      double *Ai,
                      double *Br,
                      double *Bi);

#endif

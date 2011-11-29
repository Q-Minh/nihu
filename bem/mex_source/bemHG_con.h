#ifndef BEMHG_CON_H
#define BEMHG_CON_H
#include "types.h"

void matrix_surf_const(int nnodes,
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

void matrix_field_const(int nnodes,
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

void matrix_surf_const_sparse(int nnodes,
                              const double *nodes,
                              int nelements,
                              const double *elements,
                              int npairs,
                              const double *pairs,
                              const gauss_t *g3,
                              const gauss_t *g4,
                              const double *dist,
                              double k,
                              double *Ar,
                              double *Ai,
                              double *Br,
                              double *Bi);

void matrix_field_const_sparse(int nnodes,
                               const double *nodes,
                               int nelements,
                               const double *elements,
                               int npoints,
                               const double *points,
                               int npairs,
                               const double *pairs,
                               const gauss_t *g3,
                               const gauss_t *g4,
                               const double *dist,
                               double k,
                               double *Ar,
                               double *Ai,
                               double *Br,
                               double *Bi);
#endif

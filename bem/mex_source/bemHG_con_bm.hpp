#include "bemHG_con_bm.h"

#include "types.h"
#include "green.hpp"

void matrix_surf_const_bm(int nnodes,
                          const double *nodes,
                          int nelements,
                          const double *elements,
                          const gauss_t *g3,
                          const gauss_t *g4,
                          const double *dist,
                          const complex_scalar &k,
                          const complex_scalar &alpha,
                          double *Ar,
                          double *Ai,
                          double *Br,
                          double *Bi);


void matrix_surf_const_bm_sparse(int nnodes,
                                 const double *nodes,
                                 int nelements,
                                 const double *elements,
                                 int npairs,
                                 const double *pairs,
                                 const gauss_t *g3,
                                 const gauss_t *g4,
                                 const double *dist,
                                 const complex_scalar & k,
                                 const complex_scalar &alpha,
                                 double *Ar,
                                 double *Ai,
                                 double *Br,
                                 double *Bi);

#ifndef IBEMD_CON_H
#define IBEMD_CON_H
#include "types.h"

void iBEM_D_surf_const(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       const gauss_t *g3,
                       const gauss_t *g4,
                       const double *dist,
                       double k,
                       double *Br,
                       double *Bi);

#endif

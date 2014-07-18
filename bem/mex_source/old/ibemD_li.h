#ifndef IBEMD_LI_H
#define IBEMD_LI_H
#include "types.h"

void iBEM_D_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
                     double *Dr,
                     double *Di);

#endif

#ifndef INTEGRAL_DIRECT_H
#define INTEGRAL_DIRECT_H
#include "types.h"

/* Regular integral over a linear LINE element */
void int_line_lin(const gauss2D_t *g,
                 const double *nodes,
                 const accelerator2D_t *accelerator,
                 const double *q,
                 double k,
                 double *ar,
                 double *ai,
                 double *br,
                 double *bi);

/* Singular integral over a linear LINE element */
void int_line_lin_sing(const gauss2D_t *g,
                 const double *nodes,
                 const accelerator2D_t *accelerator,
                 int corner,
                 double k,
                 double *ar,
                 double *ai,
                 double *br,
                 double *bi);

#endif

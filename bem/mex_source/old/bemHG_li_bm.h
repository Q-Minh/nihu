#ifndef BEMHG_LI_H
#define BEMHG_LI_H
#include "types.h"

/* Peter Rucz, 2012-06-20 */

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model - Burton Miller*/
void matrix_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
					 double alphar,
					 double alphai,
                     double *Ar,
                     double *Ai,
                     double *Br,
                     double *Bi);

#endif

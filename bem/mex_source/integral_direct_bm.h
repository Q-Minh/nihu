#ifndef INTEGRAL_DIRECT_H
#define INTEGRAL_DIRECT_H
#include "types.h"

/* Burton - Miller integral over linear QUAD elements */

/* Regular integral over a linear QUAD element */
/* TODO: give alpha parameter */
void int_quad_lin_bm(const gauss_t *g,
			const double *nodes,
			const accelerator_t *accelerator,
			const double *q,
			const double *nq,
			double k,
			double alphar,
			double alphai,
			double *ar,
			double *ai,
			double *br,
			double *bi);

/* Singular integral over a linear QUAD element */
void int_quad_lin_sing_bm(const gauss_t *g,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       int corner,
					   double *nq, 					/* Normal vector at source point */
                       double k,
					   double alphar, 
					   double alphai,
                       double *ar,
                       double *ai,
                       double *br,
                       double *bi);

#endif
#ifndef INTEGRAL_INDIRECT_H
#define INTEGRAL_INDIRECT_H

#include "types.h"

/* Regular double B integral over a constant TRIA element */
void dint_tritri_const(const gauss_t *gx,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi);

/* Regular double B integral over a constant QUAD element */
void dint_quadquad_const(const gauss_t *gx,
                         const double *nodesx,
                         const accelerator_t *acceleratorx,
                         const double *nodesy,
                         const accelerator_t *acceleratory,
                         double k,
                         double *br,
                         double *bi);

/* Regular double B integral over a constant TRIA and a QUAD element */
void dint_triquad_const(const gauss_t *gx,
                        const double *nodesx,
                        const accelerator_t *acceleratorx,
                        const gauss_t *gy,
                        const double *nodesy,
                        const accelerator_t *acceleratory,
                        double k,
                        double *br,
                        double *bi);

/* Singular double B integral over a constant QUAD element */
void dint_quad_const_sing(const gauss_t *g4,
                          const double *nodes,
                          const accelerator_t *accelerator,
                          double k,
                          double *br,
                          double *bi);

/* Singular double B integral over a constant TRIA element */
void dint_tri_const_sing(const gauss_t *g3,
                         const double *nodes,
                         const accelerator_t *accelerator,
                         const gauss_t *g4,
                         double k,
                         double *br,
                         double *bi);

/* Regular double B integral over linear TRIA elements */
void dint_tritri_lin(const gauss_t *g3,
                     const double *nodesx,
                     const accelerator_t *acceleratorx,
                     const double *nodesy,
                     const accelerator_t *acceleratory,
                     double k,
                     double *br,
                     double *bi);

/* Regular double B integral over linear TRIA elements */
void dint_tritri_lin_fast(const gauss_t *g3,
                     const double *nodesx,
                     const accelerator_t *acceleratorx,
                     const double *nodesy,
                     const accelerator_t *acceleratory,
                     double k,
                     double *br,
                     double *bi);

/* Regular double B integral over constant QUAD elements */
void dint_quadquad_lin(const gauss_t *g4,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi);

/* Regular double B integral over a lin TRIA and a QUAD element */
void dint_triquad_lin(const gauss_t *g3,
                      const double *nodesx,
                      const accelerator_t *acceleratorx,
                      const gauss_t *g4,
                      const double *nodesy,
                      const accelerator_t *acceleratory,
                      double k,
                      double *br,
                      double *bi);

/* Singular double B integral over a lin TRIA element */
void dint_tri_lin_sing(const gauss_t *g3,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       const gauss_t *g4,
                       double k,
                       double *br,
                       double *bi);

/* Singular double B integral over a lin QUAD element */
void dint_quad_lin_sing(const gauss_t *g4,
                        const double *nodes,
                        const accelerator_t *accelerator,
                        double k,
                        double *br,
                        double *bi);
						
						
/* Regular double D integral over a constant TRIA element */
void dintD_tritri_const(const gauss_t *gx,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi);

/* Regular double D integral over a constant QUAD element */
void dintD_quadquad_const(const gauss_t *gx,
                         const double *nodesx,
                         const accelerator_t *acceleratorx,
                         const double *nodesy,
                         const accelerator_t *acceleratory,
                         double k,
                         double *br,
                         double *bi);

/* Regular double D integral over a constant TRIA and a QUAD element */
void dintD_triquad_const(const gauss_t *gx,
                        const double *nodesx,
                        const accelerator_t *acceleratorx,
                        const gauss_t *gy,
                        const double *nodesy,
                        const accelerator_t *acceleratory,
                        double k,
                        double *br,
                        double *bi);

/* Singular double D integral over a constant QUAD element */
void dintD_quad_const_sing(const gauss_t *g4,
                          const double *nodes,
                          const accelerator_t *accelerator,
                          double k,
                          double *br,
                          double *bi);

/* Singular double D integral over a constant TRIA element */
void dintD_tri_const_sing(const gauss_t *g3,
                         const double *nodes,
                         const accelerator_t *accelerator,
                         const gauss_t *g4,
                         double k,
                         double *br,
                         double *bi);

						

#endif

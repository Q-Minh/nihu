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

/* Regular integral over a linear TRIA element */
void int_tri_lin(const gauss_t *g,
                 const double *nodes,
                 const accelerator_t *accelerator,
                 const double *q,
                 double k,
                 double *ar,
                 double *ai,
                 double *br,
                 double *bi);

/* Singular integral over a linear TRIA element */
void int_tri_lin_sing(const gauss_t *g,
                      const double *nodes,
                      const accelerator_t *accelerator,
                      int corner,
                      double k,
                      double *ar,
                      double *ai,
                      double *br,
                      double *bi);

/* Regular integral over a linear QUAD element */
void int_quad_lin(const gauss_t *g,
                  const double *nodes,
                  const accelerator_t *accelerator,
                  const double *q,
                  double k,
                  double *ar,
                  double *ai,
                  double *br,
                  double *bi);

/* Singular integral over a linear QUAD element */
void int_quad_lin_sing(const gauss_t *g,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       int corner,
                       double k,
                       double *ar,
                       double *ai,
                       double *br,
                       double *bi);

/* Regular integral over a constant LINE element */
void int_line_const(gauss2D_t g,
                   const double *nodes,
                   const accelerator2D_t *accelerator,
                   const double *q,
                   double k,
                   double *ar,
                   double *ai,
                   double *br,
                   double *bi);

/* Regular integral over a constant TRIA element */
void int_tri_const(gauss_t g,
                   const double *nodes,
                   const accelerator_t *accelerator,
                   const double *q,
                   double k,
                   double *ar,
                   double *ai,
                   double *br,
                   double *bi);

/* Singular integral over a constant TRIA element */
void int_tri_const_sing(const gauss_t *g,
                        const double *nodes,
                        const accelerator_t *accelerator,
                        const double *q,
                        double k,
                        double *ar,
                        double *ai,
                        double *br,
                        double *bi);

/* Regular integral over a constant QUAD element */
void int_quad_const(gauss_t g,
                    const double *nodes,
                    const accelerator_t *accelerator,
                    const double *q,
                    double k,
                    double *ar,
                    double *ai,
                    double *br,
                    double *bi);

/* Singular integral over a constant QUAD element */
void int_quad_const_sing(const gauss_t *g,
                         const double *nodes,
                         const accelerator_t *accelerator,
                         const double *q,
                         double k,
                         double *ar,
                         double *ai,
                         double *br,
                         double *bi);
#endif

#ifndef GREEN_H
#define GREEN_H

#include <complex>
typedef std::complex<double> complex_scalar;

void green(const double *r,
           const complex_scalar &k,
           complex_scalar &g,
           const double *n,
           complex_scalar &dg);

void green2D(const double *r,
             const complex_scalar &k,
             complex_scalar &g,
             const double *n,
             complex_scalar &dg);

/* 3D Green's function and its derivatives  */
void ddgreen(const double *r, /* y - x */
             const complex_scalar &k,
             const double *nx,
             const double *ny,
             complex_scalar &g, 	/* Green's function */
             complex_scalar &dgx, 	/* normal derivative with respect to x */
             complex_scalar &dgy,	/* normal derivative with respect to y */
             complex_scalar &ddg); 	/* double normal derivative */

/* 3D static Green's function and its derivatives  */
void ddgreen0(const double *r, /* y - x */
              const double *nx,
              const double *ny,
              double &g, 	/* Green's function */
              double &dgx, 	/* normal derivative with respect to x */
              double &dgy,	/* normal derivative with respect to y */
              double &ddg); 	/* double normal derivative */

#endif

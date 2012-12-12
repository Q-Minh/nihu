#ifndef GREEN_HPP
#define GREEN_HPP

#include "vector.h"

#include <complex>
typedef std::complex<double> complex_scalar;

const complex_scalar compJ = complex_scalar(0.0, 1.0);

template <typename kType>
complex_scalar greenr(const double ar,
                      const kType &k)
{
    return exp(-compJ * k*ar) / (4.0 * M_PI * k);
}


template <typename kType>
complex_scalar green(const double *r,
                     const kType &k)
{
    double ar = sqrt(dot(r, r));
    return exp(-compJ*k*ar) / (4.0 * M_PI * ar);
}


template <typename kType>
void green(const double *r,
           const kType &k,
           complex_scalar &g,
           const double *n,
           complex_scalar &dg)
{
    double ar = sqrt(dot(r, r));
    g = exp(-compJ*k*ar) / (4.0 * M_PI * ar);
    double rn = dot(r, n) / ar;
    dg = -g * (1.0 + compJ*k*ar)/ar * rn;
}

template <typename kType>
void ddgreen(const double *r, /* y - x */
             const kType &k,
             const double *nx, 	/* unit normal at x */
             const double *ny,	/* unit normal at y */
             complex_scalar &g, 		/* Green's function */
             complex_scalar &dgx, 		/* normal derivative with respect to x */
             complex_scalar &dgy,		/* normal derivative with respect to y */
             complex_scalar &ddg)
{
    /* Green's function */
    double ar2 = dot(r,r); 			/* square norm of distance vector */
    double ar = sqrt(ar2); 			/* norm of distance vector */
    g = exp(-compJ*k*ar) / (4.0 * M_PI * ar);

    /* Normal derivatives */
    double rnx = -dot(r, nx) / ar;	/* normal derivative of distance wrt x */
    double rny = dot(r, ny) / ar;		/* normal derivative of distance wrt y */

    dgx = -g * (1.0 + compJ*k*ar)/ar * rnx;
    dgy = -g * (1.0 + compJ*k*ar)/ar * rny;

    /* Double normal derivative */
    ddg = g * ((3.0/ar2 + 3.0*compJ*k/ar-k*k)*rnx*rny + 1.0/ar2*(1.0+compJ*k*ar)*dot(nx, ny));
}

#endif // GREEN_HPP

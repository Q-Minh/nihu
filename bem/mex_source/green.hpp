/**
* \file green.hpp
* \brief Green's functions
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef GREEN_HPP
#define GREEN_HPP

#include "vector.h"

#include <complex>

/**
* \brief Complex Scalar tyle
*/
typedef std::complex<double> complex_scalar;

/**
* \brief Imaginary Unit constant
*/
complex_scalar const compJ = complex_scalar(0.0, 1.0);

/**
* \brief Scalar green's function
* \param ar absolute value of r
* \param k wave number
* \returns the Green's function \f$\exp(-jkr)/4\pi r\f$
*/
template <typename kType>
complex_scalar green(double ar,
                      kType const &k)
{
    return exp(-compJ * k*ar) / (4.0 * M_PI * k);
}


/**
* \brief Scalar green's function
* \param r distance vector
* \param k wave number
* \returns the Green's function \f$\exp(-jkr)/4\pi r\f$
*/
template <typename kType>
complex_scalar green(double const r[],
                     kType const &k)
{
    double ar = sqrt(dot(r, r));
    return exp(-compJ*k*ar) / (4.0 * M_PI * ar);
}


/**
* \brief Scalar green's function and its normal derivative
* \param r distance vector
* \param k wave number
* \param g Green's function \f$\exp(-jkr)/4\pi r\f$
* \param n normal vector
* \param dg Green's function's normal derivative  \f$-g \frac{1+jkr}{r} \frac{{\bf r} \cdot {\bf n}}{r}\f$
*/
template <typename kType>
void green(double const r[],
           kType const &k,
           complex_scalar &g,
           double const n[],
           complex_scalar &dg)
{
    double ar = sqrt(dot(r, r));
    g = exp(-compJ*k*ar) / (4.0 * M_PI * ar);
    double rn = dot(r, n) / ar;
    dg = -g * (1.0 + compJ*k*ar)/ar * rn;
}

/**
* \brief Scalar green's function and its first and second order normal derivatives
* \param r distance vector \f${\bf y} - {\bf x}\f$
* \param k wave number
* \param nx normal vector
* \param ny normal vector
* \param g Green's function \f$\exp(-jkr)/4\pi r\f$
* \param dgx Green's function's normal derivative with respect to nx \f$g \frac{1+jkr}{r} \frac{{\bf r} \cdot {\bf n}_x}{r}\f$
* \param dgy Green's function's normal derivative with respect to ny \f$-g \frac{1+jkr}{r} \frac{{\bf r} \cdot {\bf n}_y}{r}\f$
* \param ddg Green's function's double normal derivative with respect to nx and ny
*/
template <typename kType>
void ddgreen(double const r[], /* y - x */
             kType const &k,
             double const nx[], 	/* unit normal at x */
             double const ny[],	/* unit normal at y */
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


/**
* \file integral_direct_bm.hpp
* \brief Numerical integration of Green's functions over BEM elements with Burton and Miller
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef INTEGRAL_DIRECT_BM_HPP
#define INTEGRAL_DIRECT_BM_HPP

#include "green.hpp"
#include "types.h"
#include "vector.h"

#include "elem_traits.hpp"

/**
* \brief Numerical integration of Green's functions over regular BEM element with Burton and Miller
* \param gau Gaussian integration structure
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param q reference node coordinates
* \param nq reference normal vector
* \param k wave number
* \param alpha coupling coefficient
* \param a A elem integral
* \param b B elem integral
*/
template <class ElemType, typename kType>
void int_const_bm(const gauss_t &gau,
                  double const nodes[],
                  accelerator_t const &accelerator,
                  double const q[],
                  double const nq[],
                  kType const &k, complex_scalar const &alpha,
                  complex_scalar &a, complex_scalar &b)
{
    enum {NDIM = 3};
    bool const isLinear = elem_traits<ElemType>::isLinear;
    unsigned const nNodes = elem_traits<ElemType>::nNodes;

    double norm[NDIM], jac;

    /* Initialize result to zero */
    a = b = 0.0;

    if (isLinear)
    {
        /* Jacobian and surface normal calculation  */
        jac = sqrt(dot(accelerator.n0, accelerator.n0));
        for (int j = 0; j < NDIM; j++)
            norm[j] = accelerator.n0[j]/jac;
    }

    /* for each gaussian integration point */
    for (int i = 0; i < gau.num; i++)
    {
        double r[NDIM], rxi[NDIM], reta[NDIM];

        /* computing integration location */
        for (int j = 0; j < NDIM; j++) 				/* for all directions */
        {
            r[j] = -q[j];
            if (!isLinear)
                rxi[j] = reta[j] = 0.0;
            for (int s = 0; s < nNodes; s++) 		/* for all vertices */
            {
                r[j] += gau.N[i+s*gau.num]*nodes[s+nNodes*j];
                if (!isLinear)
                {
                    rxi[j] += gau.Nxi[i+s*gau.num]*nodes[s+nNodes*j];
                    reta[j] += gau.Neta[i+s*gau.num]*nodes[s+nNodes*j];
                }
            }
        }

        double w;
        if (!isLinear)
        {
            /* surface normal and jacobian */
            cross(rxi, reta, norm);
            w = sqrt(dot(norm, norm));
            for (int j = 0; j < NDIM; j++)
                norm[j] /= w;
            w *= gau.w[i];
        }
        else
            w = gau.w[i];

        /* Evaluate Green function and its derivatives */
        complex_scalar g, dgx, dgy, ddg;
        ddgreen(r, k, nq, norm, g, dgx, dgy, ddg);

        /* Evaluate matrix elements */
        a += (dgy + alpha * ddg)*w;
        b += (g + alpha * dgx)*w;
    }

    if (isLinear)
    {
        /* Finally, multiply with jacobian */
        a *= jac;
        b *= jac;
    }
}


double gauss_xi_bm_sing[] =
{
    -0.995187219997022,    -0.974728555971310,    -0.938274552002733,    -0.886415527004401,
    -0.820001985973903,    -0.740124191578554,    -0.648093651936975,    -0.545421471388839,
    -0.433793507626045,    -0.315042679696164,    -0.191118867473616,    -0.064056892862606,
    0.064056892862605,    0.191118867473617,    0.315042679696164,    0.433793507626045,
    0.545421471388839,    0.648093651936975,    0.740124191578555,    0.820001985973903,
    0.886415527004401,    0.938274552002733,    0.974728555971310,    0.995187219997021
};
double gauss_w_bm_sing[] =
{
    0.012341229799987,    0.028531388628934,    0.044277438817420,    0.059298584915436,
    0.073346481411081,    0.086190161531953,    0.097618652104114,    0.107444270115965,
    0.115505668053726,    0.121670472927803,    0.125837456346828,    0.127938195346752,
    0.127938195346752,    0.125837456346828,    0.121670472927803,    0.115505668053726,
    0.107444270115965,    0.097618652104114,    0.086190161531953,    0.073346481411080,
    0.059298584915437,    0.044277438817420,    0.028531388628933,    0.012341229799988
};

enum {GNUM = sizeof(gauss_w_bm_sing)/sizeof(gauss_w_bm_sing[0])};

/**
* \brief Numerical integration of Green's functions over singular BEM element with Burton and Miller
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param q reference node coordinates
* \param nq reference normal vector
* \param k wave number
* \param alpha coupling coefficient
* \param a A elem integral
* \param b B elem integral
*/
template <class ElemType, typename kType>
void int_const_sing_bm(double const nodes[],
                       accelerator_t const &accelerator,
                       double const q[], double const nq[], 		    /* Source normal */
                       kType const &k, 	complex_scalar const &alpha,	/* coupling constant */
                       complex_scalar &a, complex_scalar &b)
{
    enum {NDIM = 3};
    bool const isLinear = elem_traits<ElemType>::isLinear;
    unsigned const nNodes = elem_traits<ElemType>::nNodes;

    /* Initialize the result */
    a = b = 0.0;

    /* Go through all three sides */
    for (int i = 0; i < nNodes; ++i)
    {
        double d[NDIM];
        /* Initialize nodes */
        int n1 = i;
        int n2 = (n1+1)%nNodes; 			/* Calculate between ith and i+1th nodes */
        /* Obtain distance vector */
        for (int j = 0; j < NDIM; ++j) /* For all dimensions */
            d[j] = nodes[n2+nNodes*j] - nodes[n1+nNodes*j];
        /* Calculate element length */
        double L = sqrt(dot(d,d));

        /* Go through all integration points */
        for (int ig = 0; ig < GNUM; ig++)
        {
            /* Calculate actual location x(\xi)-x_q*/
            double r[NDIM];
            for (int j=0; j < NDIM; ++j)
                r[j] = 0.5*(1.0-gauss_xi_bm_sing[ig])*nodes[n1+nNodes*j]
                       + 0.5*(1.0+gauss_xi_bm_sing[ig])*nodes[n2+nNodes*j] - q[j];
            /* Absolute value of distance */
            double lr = sqrt(dot(r,r));
            /* Jacobian is sin(beta)/ar*L/2 */
            /* Weight is also part of jacobian */
            double tmp = dot(r,d)/(lr*L);
            double jac = gauss_w_bm_sing[ig] * sqrt(1.0 - tmp*tmp) / lr * L / 2.0;

            /* Here the calculation of the integrand should be performed */
            /* Matrix H: the negative of the simple green function should be evaluated */
            a -= (green(r, k) * alpha)*jac;
            /* Matrix G: the reduced Green is evaluated */
            b += compJ * green(lr, k) * jac;
        }
    }
}

#endif


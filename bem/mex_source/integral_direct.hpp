/**
* \file integral_direct.hpp
* \brief Numerical integration of Green's functions over BEM elements
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef INTEGRAL_DIRECT_HPP
#define INTEGRAL_DIRECT_HPP

#include "elem_traits.hpp"

#include "green.hpp"
#include "quadrature.hpp"
#include "element.hpp"
#include "types.h"
#include "vector.h"


/**
* \brief Regular integral of Green's functions over constant BEM element
* \param gau Gaussian integration structure
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param q reference node coordinates
* \param k wave number
* \param a A elem integral
* \param b B elem integral
*/
template <class ElemType, typename kType>
void int_const(const gauss_t &gau,
        double const nodes[],
        accelerator_t const &accelerator,
        double const q[],
        kType const &k,
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
        
        complex_scalar g, dg;
        green(r, k, g, norm, dg);
        
        /* Evaluate matrix elements */
        a += dg*w;
        b += g*w;
    }
    
    if (isLinear)
    {
        /* Finally, multiply with jacobian */
        a *= jac;
        b *= jac;
    }
}

/**
* \brief Singular integral of Green's functions over constant BEM element
* \param gau Gaussian integration structure
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param q reference node coordinates
* \param k wave number
* \param a A elem integral
* \param b B elem integral
*/
template <class ElemType, typename kType>
void int_const_sing(const gauss_t &gau,
        double const nodes[],
        accelerator_t const &accelerator,
        double const q[],
        kType const &k,
        complex_scalar &a, complex_scalar &b)
{
    enum {NDIM = 3};
    bool const isLinear = elem_traits<ElemType>::isLinear;
    unsigned const nNodes = elem_traits<ElemType>::nNodes;
    
    double norm[NDIM], jac;
    
    a = b = 0.0;
    
    if (isLinear)
    {
        jac = sqrt(dot(accelerator.n0, accelerator.n0));
        for (int j = 0; j < NDIM; j++)
            norm[j] = accelerator.n0[j] / jac;
    }
    
    double *xiprime = new double [nNodes * gau.num];
    double *etaprime = new double [nNodes * gau.num];
    double *wprime = new double [nNodes * gau.num];
    
    int gnum;
    if (nNodes == 3)
        sing_quadr_face<ElemType>(gau, 1.0/3.0, 1.0/3.0, &gnum, xiprime, etaprime, wprime);
    else if (nNodes == 4)
        sing_quadr_face<ElemType>(gau, 0.0, 0.0, &gnum, xiprime, etaprime, wprime);
    
    /* for each gaussian integration point */
    for (int i = 0; i < gnum; i++)
    {
        double N[nNodes], r[NDIM];
        shapefun<ElemType>(xiprime[i], etaprime[i], N);
        
        /* for each coordinate direction */
        for (int j = 0; j < NDIM; j++)
        {
            r[j] = -q[j];
            /* computing integration location and its derivatives */
            for (int s = 0; s < nNodes; s++)
                r[j] += N[s]*nodes[s+nNodes*j];
            if (!isLinear)
                norm[j] = accelerator.n0[j] + accelerator.nxi[j] * xiprime[i] + accelerator.neta[j]*etaprime[i];
        }
        
        double w;
        if (!isLinear)
        {
            w = sqrt(dot(norm, norm));
            for (int j = 0; j < NDIM; j++)
                norm[j] /= w;
            w *= wprime[i];
        }
        else
            w = wprime[i];
        
        complex_scalar g, dg;
        green(r, k, g, norm, dg);
        
        a += dg*w;
        b += g*w;
    }
    
    if (isLinear)
    {
        a *= jac;
        b *= jac;
    }
    
    delete [] xiprime;
    delete [] etaprime;
    delete [] wprime;
}


/**
* \brief Regular integral of Green's functions over a linear BEM element
* \param gau Gaussian integration structure
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param q reference node coordinates
* \param k wave number
* \param a A elem integrals
* \param b B elem integrals
*/
template <class ElemType, typename kType>
void int_lin(gauss_t const &gau,
        double const *nodes,
        accelerator_t const &accelerator,
        double const *q,
        kType const &k,
        complex_scalar a[], complex_scalar b[])
{
    enum {NDIM = 3};
    bool const isLinear = elem_traits<ElemType>::isLinear;
    unsigned const nNodes = elem_traits<ElemType>::nNodes;
    
    for (int s = 0; s < nNodes; s++)
        a[s] = b[s] = 0.0;
    
    /* for each gaussian integration point */
    for (int i = 0; i < gau.num; i++)
    {
        double r[NDIM], rxi[NDIM], reta[NDIM], norm[NDIM];
        
        /* for each coordinate direction */
        for (int j = 0; j < NDIM; j++)
        {
            r[j] = -q[j];
            rxi[j] = reta[j] = 0.0;
            /* computing integration location and its derivatives */
            for (int s = 0; s < nNodes; s++)
            {
                r[j] += gau.N[i+s*gau.num]*nodes[s+nNodes*j];
                rxi[j] += gau.Nxi[i+s*gau.num]*nodes[s+nNodes*j];
                reta[j] += gau.Neta[i+s*gau.num]*nodes[s+nNodes*j];
            }
        }
        /* surface normal and jacobian */
        cross(rxi, reta, norm);
        double jac = sqrt(dot(norm, norm));
        for (int j = 0; j < NDIM; j++)
            norm[j] /= jac;
        jac *= gau.w[i];
        
        complex_scalar g, dg;
        green(r, k, g, norm, dg);
        
        g *= jac;
        dg *= jac;
        
        for (int s = 0; s < nNodes; s++)
        {
            a[s] += gau.N[i+s*gau.num]*dg;
            b[s] += gau.N[i+s*gau.num]*g;
        }
    }
}

/**
* \brief Singular integral of Green's functions over a linear BEM element
* \param gau Gaussian integration structure
* \param nodes node coordinates
* \param accelerator integration accelerator structure
* \param corner singular corner index
* \param k wave number
* \param a A elem integrals
* \param b B elem integrals
*/
template <class ElemType, typename kType>
void int_lin_sing(gauss_t const &gau,
        double const nodes[],
        accelerator_t const &accelerator,
        int corner,
        kType const &k,
        complex_scalar a[], complex_scalar b[])
{
    enum {NDIM = 3};
    bool const isLinear = elem_traits<ElemType>::isLinear;
    unsigned const nNodes = elem_traits<ElemType>::nNodes;
    
    for (int s = 0; s < nNodes; s++)
        a[s] = b[s] = 0.0;
    
    double *xiprime = new double [2*gau.num];
    double *etaprime = new double [2*gau.num];
    double *wprime = new double [2*gau.num];
    
    int gnum;
    sing_quadr_corner<ElemType>(gau, corner, &gnum, xiprime, etaprime, wprime);
    
    /* for each gaussian integration point */
    for (int i = 0; i < gnum; i++)
    {
        double r[NDIM], norm[NDIM];
        double N[nNodes];
        
        shapefun<ElemType>(xiprime[i], etaprime[i], N);
        
        /* computing integration location */
        for (int j = 0; j < NDIM; j++)
        {
            r[j] = -nodes[corner+nNodes*j];
            for (int s = 0; s < nNodes; s++)
                r[j] += N[s]*nodes[s+nNodes*j];
            
            norm[j] = accelerator.n0[j] + accelerator.nxi[j] * xiprime[i] + accelerator.neta[j] * etaprime[i];
        }
        
        double jac = sqrt(dot(norm, norm));
        for (int j = 0; j < NDIM; j++)
            norm[j] /= jac;
        jac *= wprime[i];
        
        complex_scalar g, dg;
        green(r, k, g, norm, dg);
        
        g *= jac;
        dg *= jac;
        
        for (int s = 0; s < nNodes; s++)
        {
            a[s] += N[s]*dg;
            b[s] += N[s]*g;
        }
    }
    
    delete [] xiprime;
    delete [] etaprime;
    delete [] wprime;
}


#endif


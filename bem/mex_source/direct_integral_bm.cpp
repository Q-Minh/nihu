#include "green.hpp"
#include "types.h"
#include "vector.h"

template <class elemType>
struct elem_traits;

class TriElem;
class QuadElem;

template<>
struct elem_traits<TriElem>
{
    enum
    {
        nNodes = 3,
        isLinear = true
    };
};


template<>
struct elem_traits<QuadElem>
{
    enum
    {
        nNodes = 4,
        isLinear = false
    };
};

template <class ElemType>
void int_const_bm(const gauss_t *gau,
                  const double *nodes,
                  const accelerator_t *accelerator,
                  const double *q,					/* source location */
                  const double *nq, 					/* source normal vector */
                  const complex_scalar &k, 							/* Wave number */
                  const complex_scalar &alpha,						/* Coupling real */
                  complex_scalar &a,
                  complex_scalar &b)
{
    enum {NDIM = 3};
    const bool isLinear = elem_traits<ElemType>::isLinear;
    const unsigned nNodes = elem_traits<ElemType>::nNodes;

    double norm[NDIM], jac;

    /* Initialize result to zero */
    a = b = 0.0;

    if (isLinear)
    {
        /* Jacobian and surface normal calculation  */
        jac = sqrt(dot(accelerator->n0, accelerator->n0));
        for (int j = 0; j < NDIM; j++)
            norm[j] = accelerator->n0[j]/jac;
    }

    /* for each gaussian integration point */
    for (int i = 0; i < gau->num; i++)
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
                r[j] += gau->N[i+s*gau->num]*nodes[s+nNodes*j];
                if (!isLinear)
                {
                    rxi[j] += gau->Nxi[i+s*gau->num]*nodes[s+nNodes*j];
                    reta[j] += gau->Neta[i+s*gau->num]*nodes[s+nNodes*j];
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
            w *= gau->w[i];
        }
        else
            w = gau->w[i];

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

/**
* \file bemHG_const_bm.hpp
* \brief Surface and Field point system matrices of a constant bem model with Burton Miller
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef BEMHG_COM_BM_HPP
#define BEMHG_COM_BM_HPP

#include "mesh.h"
#include "integral_direct_bm.hpp"
#include "vector.h"

/**
* \brief Surface system matrix of a constant bem model
* \param nnodes number of nodes
* \param nodes node coordinates
* \param nelements number of elements
* \param elements element node indices
* \param g3 Gaussian structures for Tria elements
* \param g4 Gaussian structures for Quad elements
* \param dist distance vector
* \param k wave number
* \param alpha coupling coefficient
* \param Ar real part of matrix A
* \param Ai imaginary part of matrix A
* \param Br real part of matrix B
* \param Bi imaginary part of matrix B
*/
template <typename kType>
void matrix_surf_const_bm(int nnodes, const double nodes[],
                          int nelements, const double elements[],
                          const gauss_t g3[], const gauss_t g4[],
                          const double dist[],
                          const kType &k, const complex_scalar &alpha,
                          double Ar[], double Ai[], double Br[], double Bi[])
{
    enum {MAXNVERT = 4, NDIM = 3};

    /* Compute element centres */
    accelerator_t *accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (int n = 0; n < nelements; n++)
    {
        /* Source location is element center */
        const double *q = accelerators[n].center;

        /* Source normal calculation */
        double nq[NDIM];
        double jac = sqrt(dot(accelerators[n].n0, accelerators[n].n0));
        for (int j = 0; j < NDIM; j++)
            nq[j] = accelerators[n].n0[j]/jac;

        /* Integration for each element */
        for (int e = 0; e < nelements; e++)
        {
            int elem[MAXNVERT];
            double nod[MAXNVERT*NDIM];

            /* number of vertices */
            int nvert = (int)elements[e];
            /* Collect element vertex nodes and coordinates */
            for (int s = 0; s < nvert; s++)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (int j = 0; j < NDIM; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }

            complex_scalar a, b;
            /* Singular element */
            if (e == n)
            {
                switch(nvert)
                {
                case 4:
                    int_const_sing_bm<QuadElem>(nod, accelerators[e], q, nq, k, alpha, a, b);
                    break;
                case 3:
                    int_const_sing_bm<TriaElem>(nod, accelerators[e], q, nq, k, alpha, a, b);
                    break;
                }
            }
            /* Non-singular element */
            else
            {
                int gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                case 4:
                    int_const_bm<QuadElem>(g4[gs], nod, accelerators[e], q, nq, k, alpha, a, b);
                    break;
                case 3:
                    int_const_bm<TriaElem>(g3[gs], nod, accelerators[e], q, nq, k, alpha, a, b);
                    break;
                }
            }
            Ar[n+nelements*e] = a.real();
            Ai[n+nelements*e] = a.imag();
            Br[n+nelements*e] = b.real();
            Bi[n+nelements*e] = b.imag();
        }
    }

    delete[] accelerators;
}


/**
* \brief Sparse surface system matrix of a constant bem model
* \param nnodes number of nodes
* \param nodes node coordinates
* \param nelements number of elements
* \param elements element node indices
* \param npairs number of element-element pairs
* \param pairs pair indices
* \param g3 Gaussian structures for Tria elements
* \param g4 Gaussian structures for Quad elements
* \param dist distance vector
* \param k wave number
* \param alpha coupling coefficient
* \param Ar real part of matrix A
* \param Ai imaginary part of matrix A
* \param Br real part of matrix B
* \param Bi imaginary part of matrix B
*/
template <typename kType>
void matrix_surf_const_bm_sparse(int nnodes, const double nodes[],
                                 int nelements, const double elements[],
                                 int npairs, const double pairs[],
                                 const gauss_t g3[], const gauss_t g4[],
                                 const double dist[],
                                 const kType& k, const complex_scalar &alpha,
                                 double Ar[], double Ai[], double Br[], double Bi[])
{
    enum {MAXNVERT = 4, NDIM = 3};

    /* Allocate space for element centres */
    accelerator_t *accelerators = new accelerator_t [nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (int p = 0; p < npairs; p++)
    {
        int n = (int)pairs[p];
        int e = (int)pairs[p+npairs];

        /* reference location */
        const double *q = accelerators[n].center;
        /* Source normal calculation */
        double nq[NDIM];
        double jac = sqrt(dot(accelerators[n].n0, accelerators[n].n0));
        for (int j = 0; j < NDIM; j++)
            nq[j] = accelerators[n].n0[j]/jac;

        int nvert = (int)elements[e];

        int elem[MAXNVERT];
        double nod[MAXNVERT*NDIM];
        for (int s = 0; s < nvert; s++)
        {
            elem[s] = (int)elements[e+(s+1)*nelements];
            for (int j = 0; j < NDIM; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }

        complex_scalar a, b;
        if (e == n) /* singular case */
        {
            switch(nvert)
            {
            case 4:
                int_const_sing_bm<QuadElem>(nod, accelerators[e], q, nq, k, alpha, a, b);
                break;
            case 3:
                int_const_sing_bm<TriaElem>(nod, accelerators[e], q, nq, k, alpha, a, b);
                break;
            }
        }
        else /* regular case */
        {
            int gs = gauss_division(q, accelerators[e].center, dist);
            switch(nvert)
            {
            case 4:
                int_const_bm<QuadElem>(g4[gs], nod, accelerators[e], q, nq, k, alpha, a, b);
                break;
            case 3:
                int_const_bm<TriaElem>(g3[gs], nod, accelerators[e], q, nq, k, alpha, a, b);
                break;
            }
        }

        Ar[p] = a.real();
        Ai[p] = a.imag();
        Br[p] = b.real();
        Bi[p] = b.imag();
    }
    delete[] accelerators;
}

#endif // BEM_HG_CON_BM_HPP


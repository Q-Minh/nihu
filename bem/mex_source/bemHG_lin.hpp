#ifndef BEMHG_LIN_HPP
#define BEMHG_LIN_HPP

#include "elem_traits.hpp"

#include "integral_direct.hpp"
#include "vector.h"
#include "mesh.h"

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
template <typename kType>
        void matrix_surf_lin(int nnodes,
        double const *nodes,
        int nelements,
        double const *elements,
        gauss_t const *g3,
        gauss_t const *g4,
        double const *dist,
        kType const &k,
        double *Ar,
        double *Ai,
        double *Br,
        double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};
    
    /* Compute element centres */
    accelerator_t *accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);
    
    /* Clear output matrices */
    for (int n = 0; n < nnodes; ++n)
        for (int j = 0; j < nnodes; ++j)
            Ar[n+nnodes*j] = Ai[n+nnodes*j] = Br[n+nnodes*j] = Bi[n+nnodes*j] = 0.0;
    
    /* Integration for each node as reference point */
    for (int n = 0; n < nnodes; ++n)
    {
        /* reference location */
        double q[NDIM];
        for (int j = 0; j < 3; ++j)
            q[j] = nodes[n+j*nnodes];
        
        /* Integration for each element */
        for (int e = 0; e < nelements; ++e)
        {
            int elem[MAXNVERT];
            double nod[MAXNVERT*NDIM];
            
            int nvert = (int)elements[e];
            
            /* Collect element vertex nodes and coordinates */
            bool sing = false;
            int corner;
            complex_scalar a[MAXNVERT], b[MAXNVERT];
            for (int s = 0; s < nvert; ++s)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                if (elem[s] == n)
                {
                    sing = true;
                    corner = s;
                }
            }
            for (int s = 0; s < nvert; ++s)
                for (int j = 0; j < NDIM; ++j)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            if (sing)
            {
                switch(nvert)
                {
                    case 3:
                        int_lin_sing<TriaElem, kType>(g4[0], nod, accelerators[e], corner, k, a, b);
                        break;
                    case 4:
                        int_lin_sing<QuadElem, kType>(g4[0], nod, accelerators[e], corner, k, a, b);
                        break;
                }
            }
            else
            {
                int gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                    case 3:
                        int_lin<TriaElem, kType>(g3[gs], nod, accelerators[e], q, k, a, b);
                        break;
                    case 4:
                        int_lin<QuadElem, kType>(g4[gs], nod, accelerators[e], q, k, a, b);
                        break;
                }
            }
            for (int s = 0; s < nvert; ++s)
            {
                Ar[n+nnodes*elem[s]] += a[s].real();
                Ai[n+nnodes*elem[s]] += a[s].imag();
                Br[n+nnodes*elem[s]] += b[s].real();
                Bi[n+nnodes*elem[s]] += b[s].imag();
            }
        }
    }
    delete [] accelerators;
}


/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
template <typename kType>
        void matrix_field_lin(int nnodes,
        double const *nodes,
        int nelements,
        double const *elements,
        int npoints,
        double const *points,
        const gauss_t *g3,
        const gauss_t *g4,
        double const *dist,
        kType const &k,
        double *Ar,
        double *Ai,
        double *Br,
        double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};
    
    /* Compute element centres */
    accelerator_t *accelerators = new accelerator_t [nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);
    
    /* Clear output matrices */
    for (int n = 0; n < npoints; ++n)
        for (int j = 0; j < nnodes; ++j)
            Ar[n+npoints*j] = Ai[n+npoints*j] = Br[n+npoints*j] = Bi[n+npoints*j] = 0.0;
    
    /* Integration for each node as reference point */
    for (int n = 0; n < npoints; ++n)
    {
        double q[NDIM];
        /* reference location */
        for (int j = 0; j < 3; ++j)
            q[j] = points[n+j*npoints];
        
        /* Integration for each element */
        for (int e = 0; e < nelements; ++e)
        {
            int elem[MAXNVERT];
            double nod[MAXNVERT*NDIM];
            
            int nvert = (int)elements[e];
            
            /* Collect element vertex nodes and coordinates */
            for (int s = 0; s < nvert; ++s)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (int j = 0; j < 3; ++j)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
            int gs = gauss_division(q, accelerators[e].center, dist);
            complex_scalar a[MAXNVERT], b[MAXNVERT];
            switch(nvert)
            {
                case 3:
                    int_lin<TriaElem, kType>(g3[gs], nod, accelerators[e], q, k, a, b);
                    break;
                case 4:
                    int_lin<QuadElem, kType>(g4[gs], nod, accelerators[e], q, k, a, b);
                    break;
            }
            
            for (int s = 0; s < nvert; ++s)
            {
                Ar[n+npoints*elem[s]] += a[s].real();
                Ai[n+npoints*elem[s]] += a[s].imag();
                Br[n+npoints*elem[s]] += b[s].real();
                Bi[n+npoints*elem[s]] += b[s].imag();
            }
        }
    }
    
    delete [] accelerators;
}

#endif


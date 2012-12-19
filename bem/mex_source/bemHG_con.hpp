#ifndef BEMHG_CON_HPP
#define BEMHG_CON_HPP

#include "elem_traits.hpp"

#include "integral_direct.hpp"
#include "vector.h"
#include "mesh.h"

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
template <typename kType>
        void matrix_surf_const(int nnodes,
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
    
    /* Integration for each node as reference point */
    for (int n = 0; n < nelements; n++)
    {
        double const *q = accelerators[n].center;
        
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
            if (e == n)
            {
                switch(nvert)
                {
                    case 3:
                        int_const_sing<TriaElem>(g4[0], nod, accelerators[e], q, k, a, b);
                        break;
                    case 4:
                        int_const_sing<QuadElem>(g4[0], nod, accelerators[e], q, k, a, b);
                        break;
                }
            }
            else
            {
                int gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                    case 3:
                        int_const<TriaElem>(g3[gs], nod, accelerators[e], q, k, a, b);
                        break;
                    case 4:
                        int_const<QuadElem>(g4[gs], nod, accelerators[e], q, k, a, b);
                        break;
                }
            }
            Ar[n+nelements*e] = a.real();
            Ai[n+nelements*e] = a.imag();
            Br[n+nelements*e] = b.real();
            Bi[n+nelements*e] = b.imag();
        }
    }
    delete [] accelerators;
}


/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
template <typename kType>
        void matrix_field_const(int nnodes,
        double const *nodes,
        int nelements,
        double const *elements,
        int npoints,
        double const *points,
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
    
    /* Allocate space for element centres */
    accelerator_t *accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);
    
    /* Integration for each node as reference point */
    for (int n = 0; n < npoints; n++)
    {
        double q[NDIM];
        /* reference location */
        for (int j = 0; j < NDIM; j++)
            q[j] = points[n+j*npoints];
        
        /* Integration for each element */
        for (int e = 0; e < nelements; e++)
        {
            int elem[MAXNVERT];
            double nod[NDIM*MAXNVERT];
            /* number of vertices */
            int nvert = (int)elements[e];
            /* collect element vertices and coordinates */
            for (int s = 0; s < nvert; s++)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (int j = 0; j < NDIM; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            
            /* perform the integration */
            int gs = gauss_division(q, accelerators[e].center, dist);
            complex_scalar a, b;
            switch(nvert)
            {
                case 4:
                    int_const<QuadElem>(g4[gs], nod, accelerators[e], q, k, a, b);
                    break;
                case 3:
                    int_const<TriaElem>(g3[gs], nod, accelerators[e], q, k, a, b);
                    break;
            }
            Ar[n+npoints*e] = a.real();
            Ai[n+npoints*e] = a.imag();
            Br[n+npoints*e] = b.real();
            Bi[n+npoints*e] = b.imag();
        }
    }
    delete [] accelerators;
}


/* ------------------------------------------------------------------------ */
/* Compute surface sparse system matrices of a bem model */
template <typename kType>
        void matrix_surf_const_sparse(int nnodes,
        double const *nodes,
        int nelements,
        double const *elements,
        int npairs,
        double const *pairs,
        const gauss_t *g3,
        const gauss_t *g4,
        double const *dist,
        const kType &k,
        double *Ar,
        double *Ai,
        double *Br,
        double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};
    
    /* Allocate space for element centres */
    accelerator_t *accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);
    
    /* Integration for each node as reference point */
    for (int p = 0; p < npairs; p++)
    {
        int elem[MAXNVERT];
        double nod[MAXNVERT*NDIM];
        int n = (int)pairs[p];
        int e = (int)pairs[p+npairs];
        
        /* reference location */
        double const *q = accelerators[n].center;
        
        int nvert = (int)elements[e];
        
        for (int s = 0; s < nvert; s++)
        {
            elem[s] = (int)elements[e+(s+1)*nelements];
            for (int j = 0; j < NDIM; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }
        
        complex_scalar a, b;;
        if (e == n) /* singular case */
        {
            switch(nvert)
            {
                case 4:
                    int_const_sing<QuadElem>(g4[0], nod, accelerators[e], q, k, a, b);
                    break;
                case 3:
                    int_const_sing<TriaElem>(g4[0], nod, accelerators[e], q, k, a, b);
                    break;
            }
        }
        else /* regular case */
        {
            int gs = gauss_division(q, accelerators[e].center, dist);
            switch(nvert)
            {
                case 4:
                    int_quad_const(&g4[gs], nod, &accelerators[e], q, k, a, b);
                    break;
                case 3:
                    int_tri_const(&g3[gs], nod, &accelerators[e], q, k, a, b);
                    break;
            }
        }
        
        Ar[p] = a.real();
        Ai[p] = a.imag();
        Br[p] = b.real();
        Bi[p] = b.imag();
    }
    delete [] accelerators;
}

/* ------------------------------------------------------------------------ */
/* Compute far field sparse system matrices of a bem model */
template <typename kType>
        void matrix_field_const_sparse(int nnodes,
        double const *nodes,
        int nelements,
        double const *elements,
        int npoints,
        double const *points,
        int npairs,
        double const *pairs,
        const gauss_t *g3,
        const gauss_t *g4,
        double const *dist,
        double k,
        double *Ar,
        double *Ai,
        double *Br,
        double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};
    
    double q[NDIM];
    double nod[MAXNVERT*NDIM];
    int j, e, s, n, gs, p;
    int elem[MAXNVERT], nvert;
    complex_scalar a, b;;
    
    /* Allocate space for element centres */
    accelerator_t *accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);
    
    /* Integration for each node as reference point */
    for (p = 0; p < npairs; p++)
    {
        n = (int)pairs[p];
        e = (int)pairs[p+npairs];
        
        /* reference location */
        for (j = 0; j < NDIM; j++)
            q[j] = points[n+j*npoints];
        
        nvert = (int)elements[e];
        
        for (s = 0; s < nvert; s++)
        {
            elem[s] = (int)elements[e+(s+1)*nelements];
            for (j = 0; j < NDIM; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }
        
        gs = gauss_division(q, accelerators[e].center, dist);
        switch(nvert)
        {
            case 4:
                int_quad_const(&g4[gs], nod, &accelerators[e], q, k, a, b);
                break;
            case 3:
                int_tri_const(&g3[gs], nod, &accelerators[e], q, k, a, b);
                break;
        }
        Ar[p] = a.real();
        Ai[p] = a.imag();
        Br[p] = b.imag();
        Bi[p] = b.real();
    }
    
    delete [] accelerators;
}

#endif

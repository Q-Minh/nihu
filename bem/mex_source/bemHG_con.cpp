#include "bemHG_con.hpp"

#include "integral_direct.h"
#include "vector.h"
#include "mesh.h"

#include <cmath>

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_const(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       const gauss_t *g3,
                       const gauss_t *g4,
                       const double *dist,
                       const complex<double> &k,
                       double *Ar,
                       double *Ai,
                       double *Br,
                       double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};

    accelerator_t *accelerators;
    const double *q;
    int j, e, s, n, gs;
    int elem[MAXNVERT], nvert;
    double nod[MAXNVERT*NDIM];
    complex<double> a, b;

    /* Compute element centres */
    accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (n = 0; n < nelements; n++)
    {
        q = accelerators[n].center;

        /* Integration for each element */
        for (e = 0; e < nelements; e++)
        {
            /* number of vertices */
            nvert = (int)elements[e];
            /* Collect element vertex nodes and coordinates */
            for (s = 0; s < nvert; s++)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (j = 0; j < NDIM; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            if (e == n)
            {
                switch(nvert)
                {
                case 4:
                    int_quad_const_sing(&g4[0], nod, &accelerators[e], q, k, a, b);
                    break;
                case 3:
                    int_tri_const_sing(&g4[0], nod, &accelerators[e], q, k, a, b);
                    break;
                }
            }
            else
            {
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
/* Compute surface system matrices of a bem model */
void matrix_surf_const2D(int nnodes,
                         const double *nodes,
                         int nelements,
                         const double *elements,
                         const gauss2D_t *g,
                         const double *dist,
                         const complex<double> &k,
                         double *Ar,
                         double *Ai,
                         double *Br,
                         double *Bi)
{
    enum {NDIM = 2, NVERT = 2};

    accelerator2D_t *accelerators;
    const double *q;
    int j, e, s, n, gs;
    int elem[NVERT];
    double nod[NVERT*NDIM];
    complex<double> a, b;

    /* Compute element centres */
    accelerators = new accelerator2D_t[nelements];
    init_accelerators2D(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (n = 0; n < nelements; n++)
    {
        q = accelerators[n].center;

        /* Integration for each element */
        for (e = 0; e < nelements; e++)
        {
            /* Collect element vertex nodes and coordinates */
            for (s = 0; s < NVERT; s++)
            {
                elem[s] = (int)elements[e+s*nelements];
                for (j = 0; j < NDIM; j++)
                    nod[s+NVERT*j] = nodes[elem[s]+j*nnodes];
            }
            if (e == n)
                int_line_const_sing(&g[0], nod, &accelerators[e], q, k, a, b);
            else
            {
                gs = gauss_division2D(q, accelerators[e].center, dist);
                int_line_const(&g[gs], nod, &accelerators[e], q, k, a, b);
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
void matrix_field_const(int nnodes,
                        const double *nodes,
                        int nelements,
                        const double *elements,
                        int npoints,
                        const double *points,
                        const gauss_t *g3,
                        const gauss_t *g4,
                        const double *dist,
                        const complex<double> &k,
                        double *Ar,
                        double *Ai,
                        double *Br,
                        double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};

    accelerator_t *accelerators;
    double q[NDIM];
    double nod[MAXNVERT*NDIM];
    int j, e, s, n, gs;
    int elem[MAXNVERT], nvert;
    complex<double> a, b;

    /* Allocate space for element centres */
    accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (n = 0; n < npoints; n++)
    {
        /* reference location */
        for (j = 0; j < NDIM; j++)
            q[j] = points[n+j*npoints];

        /* Integration for each element */
        for (e = 0; e < nelements; e++)
        {
            /* number of vertices */
            nvert = (int)elements[e];
            /* collect element vertices and coordinates */
            for (s = 0; s < nvert; s++)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (j = 0; j < NDIM; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
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
            Ar[n+npoints*e] = a.real();
            Ai[n+npoints*e] = a.imag();
            Br[n+npoints*e] = b.real();
            Bi[n+npoints*e] = b.imag();
        }
    }
    delete [] accelerators;
}

/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
void matrix_field_const2D(int nnodes,
                          const double *nodes,
                          int nelements,
                          const double *elements,
                          int npoints,
                          const double *points,
                          const gauss2D_t *g,
                          const double *dist,
                          const complex<double> &k,
                          double *Ar,
                          double *Ai,
                          double *Br,
                          double *Bi)
{
    enum {NDIM = 2, NVERT = 2};

    accelerator2D_t *accelerators;
    double q[NDIM];
    double nod[NVERT*NDIM];
    int j, e, s, n, gs;
    int elem[NVERT];
    complex<double> a, b;

    /* Allocate space for element centres */
    accelerators = new accelerator2D_t[nelements];
    init_accelerators2D(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (n = 0; n < npoints; n++)
    {
        /* reference location */
        for (j = 0; j < NDIM; j++)
            q[j] = points[n+j*npoints];

        /* Integration for each element */
        for (e = 0; e < nelements; e++)
        {
            /* collect element vertices and coordinates */
            for (s = 0; s < NVERT; s++)
            {
                elem[s] = (int)elements[e+s*nelements];
                for (j = 0; j < NDIM; j++)
                    nod[s+NVERT*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
            gs = gauss_division2D(q, accelerators[e].center, dist);
            int_line_const(&g[gs], nod, &accelerators[e], q, k, a, b);
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
void matrix_surf_const_sparse(int nnodes,
                              const double *nodes,
                              int nelements,
                              const double *elements,
                              int npairs,
                              const double *pairs,
                              const gauss_t *g3,
                              const gauss_t *g4,
                              const double *dist,
                              const complex<double> &k,
                              double *Ar,
                              double *Ai,
                              double *Br,
                              double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};

    accelerator_t *accelerators;
    const double *q;
    int j, e, s, n, gs, p;
    int elem[MAXNVERT], nvert;
    double nod[MAXNVERT*NDIM];
    complex<double> a, b;

    /* Allocate space for element centres */
    accelerators = new accelerator_t[nelements];
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (p = 0; p < npairs; p++)
    {
        n = (int)pairs[p];
        e = (int)pairs[p+npairs];

        /* reference location */
        q = accelerators[n].center;

        nvert = (int)elements[e];

        for (s = 0; s < nvert; s++)
        {
            elem[s] = (int)elements[e+(s+1)*nelements];
            for (j = 0; j < NDIM; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }

        if (e == n) /* singular case */
        {
            switch(nvert)
            {
            case 4:
                int_quad_const_sing(&g4[0], nod, &accelerators[e], q, k, a, b);
                break;
            case 3:
                int_tri_const_sing(&g4[0], nod, &accelerators[e], q, k, a, b);
                break;
            }
        }
        else /* regular case */
        {
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
void matrix_field_const_sparse(int nnodes,
                               const double *nodes,
                               int nelements,
                               const double *elements,
                               int npoints,
                               const double *points,
                               int npairs,
                               const double *pairs,
                               const gauss_t *g3,
                               const gauss_t *g4,
                               const double *dist,
                               const complex<double> &k,
                               double *Ar,
                               double *Ai,
                               double *Br,
                               double *Bi)
{
    enum {MAXNVERT = 4, NDIM = 3};

    accelerator_t *accelerators;
    double q[NDIM];
    double nod[MAXNVERT*NDIM];
    int j, e, s, n, gs, p;
    int elem[MAXNVERT], nvert;
    complex<double> a, b;

    /* Allocate space for element centres */
    accelerators = new accelerator_t[nelements];
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
        Br[p] = b.real();
        Bi[p] = b.imag();
    }

    delete [] accelerators;
}

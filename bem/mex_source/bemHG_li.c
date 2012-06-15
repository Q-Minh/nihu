#include "bemHG_li.h"

#include "integral_direct.h"
#include "mesh.h"

#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
                     double *Ar,
                     double *Ai,
                     double *Br,
                     double *Bi)
{
    double q[3];
    int j, e, s, n;
    int elem[4], nvert;
    double nod[12];
    double ai[4], bi[4], ar[4], br[4];
    accelerator_t* accelerators;
    boolean sing;
    int corner, gs;

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (n = 0; n < nnodes; ++n)
        for (j = 0; j < nnodes; ++j)
            Ar[n+nnodes*j] = Ai[n+nnodes*j] = Br[n+nnodes*j] = Bi[n+nnodes*j] = 0.0;

    /* Integration for each node as reference point */
    for (n = 0; n < nnodes; ++n)
    {
        /* reference location */
        for (j = 0; j < 3; ++j)
            q[j] = nodes[n+j*nnodes];

        /* Integration for each element */
        for (e = 0; e < nelements; ++e)
        {
            nvert = (int)elements[e];

            /* Collect element vertex nodes and coordinates */
            sing = 0;
            for (s = 0; s < nvert; ++s)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                if (elem[s] == n)
                {
                    sing = 1;
                    corner = s;
                }
            }
            for (s = 0; s < nvert; ++s)
                for (j = 0; j < 3; ++j)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            if (sing)
            {
                switch(nvert)
                {
                case 3:
                    int_tri_lin_sing(&g4[0], nod, &accelerators[e], corner, k, ar, ai, br, bi);
                    break;
                case 4:
                    int_quad_lin_sing(&g4[0], nod, &accelerators[e], corner, k, ar, ai, br, bi);
                    break;
                }
            }
            else
            {
                gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                case 3:
                    int_tri_lin(&g3[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
                    break;
                case 4:
                    int_quad_lin(&g4[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
                    break;
                }
            }
            for (s = 0; s < nvert; ++s)
            {
                Ar[n+nnodes*elem[s]] += ar[s];
                Ai[n+nnodes*elem[s]] += ai[s];
                Br[n+nnodes*elem[s]] += br[s];
                Bi[n+nnodes*elem[s]] += bi[s];
            }
        }
    }
    free(accelerators);
}


/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_lin2D(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g,
                     const double *dist,
                     double k,
                     double *Ar,
                     double *Ai,
                     double *Br,
                     double *Bi)
{
#define NDIM 2
#define NVERT 2
    double q[NDIM];
    int j, e, s, n;
    int elem[NVERT];
    double nod[NVERT*NDIM];
    double ai[NVERT], bi[NVERT], ar[NVERT], br[NVERT];
    accelerator2D_t* accelerators;
    boolean sing;
    int corner, gs;

    /* Compute element centres */
    accelerators = (accelerator2D_t *)calloc(nelements, sizeof(accelerator2D_t));
    init_accelerators2D(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (n = 0; n < nnodes; ++n)
        for (j = 0; j < nnodes; ++j)
            Ar[n+nnodes*j] = Ai[n+nnodes*j] = Br[n+nnodes*j] = Bi[n+nnodes*j] = 0.0;

    /* Integration for each node as reference point */
    for (n = 0; n < nnodes; ++n)
    {
        /* reference location */
        for (j = 0; j < NDIM; ++j)
            q[j] = nodes[n+j*nnodes];

        /* Integration for each element */
        for (e = 0; e < nelements; ++e)
        {
            /* Collect element vertex nodes and coordinates */
            sing = 0;
            for (s = 0; s < NVERT; ++s)
            {
                elem[s] = (int)elements[e+s*nelements];
                if (elem[s] == n)
                {
                    sing = 1;
                    corner = s;
                }
            }
            for (s = 0; s < NVERT; ++s)
                for (j = 0; j < NDIM; ++j)
                    nod[s+NVERT*j] = nodes[elem[s]+j*nnodes];
            if (sing)
                int_line_lin_sing(&g[0], nod, &accelerators[e], corner, k, ar, ai, br, bi);
            else
            {
                gs = gauss_division2D(q, accelerators[e].center, dist);
                int_line_lin(&g[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
            }
            for (s = 0; s < NVERT; ++s)
            {
                Ar[n+nnodes*elem[s]] += ar[s];
                Ai[n+nnodes*elem[s]] += ai[s];
                Br[n+nnodes*elem[s]] += br[s];
                Bi[n+nnodes*elem[s]] += bi[s];
            }
        }
    }
    free(accelerators);
}


/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
void matrix_field_lin(int nnodes,
                      const double *nodes,
                      int nelements,
                      const double *elements,
                      int npoints,
                      const double *points,
                      const gauss_t *g3,
                      const gauss_t *g4,
                      const double *dist,
                      double k,
                      double *Ar,
                      double *Ai,
                      double *Br,
                      double *Bi)
{
    double q[3];
    int j, e, s, n, gs;
    int elem[4], nvert;
    double nod[12];
    double ai[4], bi[4], ar[4], br[4];
    accelerator_t* accelerators;

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (n = 0; n < npoints; ++n)
        for (j = 0; j < nnodes; ++j)
            Ar[n+npoints*j] = Ai[n+npoints*j] = Br[n+npoints*j] = Bi[n+npoints*j] = 0.0;

    /* Integration for each node as reference point */
    for (n = 0; n < npoints; ++n)
    {
        /* reference location */
        for (j = 0; j < 3; ++j)
            q[j] = points[n+j*npoints];

        /* Integration for each element */
        for (e = 0; e < nelements; ++e)
        {
            nvert = (int)elements[e];

            /* Collect element vertex nodes and coordinates */
            for (s = 0; s < nvert; ++s)
            {
                elem[s] = (int)elements[e+(s+1)*nelements];
                for (j = 0; j < 3; ++j)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
            gs = gauss_division(q, accelerators[e].center, dist);
            switch(nvert)
            {
            case 3:
                int_tri_lin(&g3[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
                break;
            case 4:
                int_quad_lin(&g4[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
                break;
            }

            for (s = 0; s < nvert; ++s)
            {
                Ar[n+npoints*elem[s]] += ar[s];
                Ai[n+npoints*elem[s]] += ai[s];
                Br[n+npoints*elem[s]] += br[s];
                Bi[n+npoints*elem[s]] += bi[s];
            }
        }
    }

    free(accelerators);
}

/* ------------------------------------------------------------------------ */
/* Compute far field system matrices of a bem model */
void matrix_field_lin2D(int nnodes,
                      const double *nodes,
                      int nelements,
                      const double *elements,
                      int npoints,
                      const double *points,
                      const gauss2D_t *g,
                      const double *dist,
                      double k,
                      double *Ar,
                      double *Ai,
                      double *Br,
                      double *Bi)
{
#define NDIM 2
#define NVERT 2
    double q[NDIM];
    int j, e, s, n, gs;
    int elem[NVERT], nvert;
    double nod[NVERT*NDIM];
    double ai[NVERT], bi[NVERT], ar[NVERT], br[NVERT];
    accelerator2D_t* accelerators;

    /* Compute element centres */
    accelerators = (accelerator2D_t *)calloc(nelements, sizeof(accelerator2D_t));
    init_accelerators2D(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (n = 0; n < npoints; ++n)
        for (j = 0; j < nnodes; ++j)
            Ar[n+npoints*j] = Ai[n+npoints*j] = Br[n+npoints*j] = Bi[n+npoints*j] = 0.0;

    /* Integration for each node as reference point */
    for (n = 0; n < npoints; ++n)
    {
        /* reference location */
        for (j = 0; j < NDIM; ++j)
            q[j] = points[n+j*npoints];

        /* Integration for each element */
        for (e = 0; e < nelements; ++e)
        {
            /* Collect element vertex nodes and coordinates */
            for (s = 0; s < NVERT; ++s)
            {
                elem[s] = (int)elements[e+s*nelements];
                for (j = 0; j < NDIM; ++j)
                    nod[s+NVERT*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
            gs = gauss_division2D(q, accelerators[e].center, dist);
            int_line_lin(&g[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);

            for (s = 0; s < NVERT; ++s)
            {
                Ar[n+npoints*elem[s]] += ar[s];
                Ai[n+npoints*elem[s]] += ai[s];
                Br[n+npoints*elem[s]] += br[s];
                Bi[n+npoints*elem[s]] += bi[s];
            }
        }
    }

    free(accelerators);
}

#include "bemHG_con.h"

#include "integral_direct.h"
#include "vector.h"
#include "mesh.h"

#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_const(int nnodes,
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
    accelerator_t *accelerators;
    const double *q;
    int j, e, s, n, gs;
    int elem[4], nvert;
    double nod[12];
    double ai, bi, ar, br;

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
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
                for (j = 0; j < 3; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            if (e == n)
            {
                switch(nvert)
                {
                case 4:
                    int_quad_const_sing(&g4[0], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                    break;
                case 3:
                    int_tri_const_sing(&g4[0], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                    break;
                }
            }
            else
            {
                gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                case 4:
                    int_quad_const(g4[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                    break;
                case 3:
                    int_tri_const(g3[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                    break;
                }
            }
            Ar[n+nelements*e] = ar;
            Ai[n+nelements*e] = ai;
            Br[n+nelements*e] = br;
            Bi[n+nelements*e] = bi;
        }
    }
    free(accelerators);
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
                        double k,
                        double *Ar,
                        double *Ai,
                        double *Br,
                        double *Bi)
{
    accelerator_t *accelerators;
    double q[3];
    double nod[12];
    int j, e, s, n, gs;
    int elem[4], nvert;
    double ai, bi, ar, br;

    /* Allocate space for element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each node as reference point */
    for (n = 0; n < npoints; n++)
    {
        /* reference location */
        for (j = 0; j < 3; j++)
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
                for (j = 0; j < 3; j++)
                    nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
            }
            /* perform the integration */
            gs = gauss_division(q, accelerators[e].center, dist);
            switch(nvert)
            {
            case 4:
                int_quad_const(g4[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                break;
            case 3:
                int_tri_const(g3[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
                break;
            }
            Ar[n+npoints*e] = ar;
            Ai[n+npoints*e] = ai;
            Br[n+npoints*e] = br;
            Bi[n+npoints*e] = bi;
        }
    }
    free(accelerators);
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
                              double k,
                              double *Ar,
                              double *Ai,
                              double *Br,
                              double *Bi)
{
    accelerator_t *accelerators;
    const double *q;
    int j, e, s, n, gs, p;
    int elem[4], nvert;
    double nod[12];
    double ai, bi, ar, br;

    /* Allocate space for element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
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
            for (j = 0; j < 3; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }

        if (e == n) /* singular case */
		{
			switch(nvert)
			{
            case 4:
				int_quad_const_sing(&g4[0], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
				break;
			case 3:
				int_tri_const_sing(&g4[0], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
				break;
			}
		}
        else /* regular case */
		{
			gs = gauss_division(q, accelerators[e].center, dist);
			switch(nvert)
			{
			case 4:
				int_quad_const(g4[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
				break;
			case 3:
				int_tri_const(g3[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
				break;
			}
		}

		Ar[p] = ar;
		Ai[p] = ai;
		Br[p] = br;
		Bi[p] = bi;
	}
    free(accelerators);
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
                               double k,
                               double *Ar,
                               double *Ai,
                               double *Br,
                               double *Bi)
{
    accelerator_t *accelerators;
    const double *q;
    double nod[12];
    int j, e, s, n, gs, p;
    int elem[4], nvert;
    double ai, bi, ar, br;

    /* Allocate space for element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
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
            for (j = 0; j < 3; j++)
                nod[s+nvert*j] = nodes[elem[s]+j*nnodes];
        }

        gs = gauss_division(q, accelerators[e].center, dist);
        switch(nvert)
        {
        case 4:
            int_quad_const(g4[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
            break;
        case 3:
            int_tri_const(g3[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi);
            break;
        }
        Ar[p] = ar;
        Ai[p] = ai;
        Br[p] = br;
        Bi[p] = bi;
    }

    free(accelerators);
}

#include "bemHG_li.h"

#include "integral_direct_bm.h"
#include "mesh.h"

#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_lin_bm(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
					 double alphar,
					 double alphai,
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
			/* Singular integrals */
            if (sing)
            {
                switch(nvert)
                {
				/*
                case 3:
                    int_tri_lin_sing(&g4[0], nod, &accelerators[e], corner, k, ar, ai, br, bi);
                    break;
				*/
                case 4:
					 /* Singular QUAD */
                    int_quad_lin_sing_bm(&g4[0], nod, &accelerators[e], corner, 
					 q, k, alphar, alphai, ar, ai, br, bi); */
                    break;
				default:
					/* Nincs más elem! */
					break;
                }
            }
            /* Non-singular integrals */
            else
            {
                gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
				/* 
				case 3:
					int_tri_lin(&g3[gs], nod, &accelerators[e], q, k, ar, ai, br, bi);
					break;
				*/
                case 4:
                    int_quad_lin_bm(&g4[gs], nod, &accelerators[e], q, q, 
					k, alphar, alphai, ar, ai, br, bi);
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
#include "bemHG_con_bm.h"

#include "integral_direct_bm.h"
#include "vector.h"
#include "mesh.h"

#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system matrices of a bem model */
void matrix_surf_const_bm(int nnodes,
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
                    /* int_quad_const_sing(&g4[0], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi); */
                    break;
                case 3:
                    int_tri_const_sing_bm(&g4[0], nod, &accelerators[e], q, q, k, 0, 1.0/k, &ar, &ai, &br, &bi);
                    break;
                }
            }
            else
            {
                gs = gauss_division(q, accelerators[e].center, dist);
                switch(nvert)
                {
                case 4:
                    /*int_quad_const(&g4[gs], nod, &accelerators[e], q, k, &ar, &ai, &br, &bi); */
                    break;
                case 3:
                    int_tri_const_bm(&g3[gs], nod, &accelerators[e], q, q, k, 0, 1.0/k,  &ar, &ai, &br, &bi);
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
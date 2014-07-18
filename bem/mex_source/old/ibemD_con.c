#include "ibemD_con.h"

#include "integral_indirect.h"
#include "vector.h"
#include "mesh.h"
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system B matrices of an indirect bem model */
void iBEM_D_surf_const(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       const gauss_t *g3,
                       const gauss_t *g4,
                       const double *dist,
                       double k, double *Br,
                       double *Bi)
{
    accelerator_t *accelerators;
    int j, ey, s, ex, gs;
    int elemx[4], elemy[4], nvertx, nverty;
    double nodx[12], nody[12];
    double bi, br;

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Integration for each element (x) as reference */
    for (ex = 0; ex < nelements; ex++)
    {
        nvertx = (int)elements[ex]; /* number of vertices of x elem */
        /* Collect element vertex nodes and coordinates */
        for (s = 0; s < nvertx; s++)
        {
            elemx[s] = (int)elements[ex+(s+1)*nelements];
            for (j = 0; j < 3; j++)
                nodx[s+nvertx*j] = nodes[elemx[s]+j*nnodes];
        }
        /* singular integration */
        switch(nvertx)
        {
        case 3:
            dintD_tri_const_sing(&g3[0], nodx, &accelerators[ex], &g4[0], k, &br, &bi);
            break;
        case 4:
            dintD_quad_const_sing(&g4[0], nodx, &accelerators[ex], k, &br, &bi);
            break;
        }
        Br[ex+nelements*ex] = br;
        Bi[ex+nelements*ex] = bi;

        /* Integration for the half of elements (y) */
        for (ey = ex+1; ey < nelements; ey++)
        {
            nverty = (int)elements[ey]; /* number of vertices of y elem */
            /* Collect element vertex nodes and coordinates */
            for (s = 0; s < nverty; s++)
            {
                elemy[s] = (int)elements[ey+(s+1)*nelements];
                for (j = 0; j < 3; j++)
                    nody[s+nverty*j] = nodes[elemy[s]+j*nnodes];
            }
            gs = gauss_division(accelerators[ex].center, accelerators[ey].center, dist);
            switch(nvertx)
            {
            case 4:
                switch (nverty)
                {
                case 4:
                    dintD_quadquad_const(&g4[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, &br, &bi);
                    break;
                case 3:
                    dintD_triquad_const(&g3[gs], nody, &accelerators[ey], &g4[gs], nodx, &accelerators[ex], k, &br, &bi);
                    break;
                }
                break;
            case 3:
                switch (nverty)
                {
                case 4:
                    dintD_triquad_const(&g3[gs], nodx, &accelerators[ex], &g4[gs], nody, &accelerators[ey], k, &br, &bi);
                    break;
                case 3:
                    dintD_tritri_const(&g3[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, &br, &bi);
                    break;
                }
            }
            Br[ex+nelements*ey] = br;
            Bi[ex+nelements*ey] = bi;
            Br[ey+nelements*ex] = br;
            Bi[ey+nelements*ex] = bi;
        }
    }
    free(accelerators);
}

#include "ibemB_li.h"

#include "integral_indirect.h"
#include "vector.h"
#include "mesh.h"
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system B matrices of an indirect bem model               */
void iBEM_B_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
                     double *Br,
                     double *Bi)
{
    accelerator_t *accelerators;
    int j, ey, s, ex, gs, ix, iy;
    int elemx[4], elemy[4], nvertx, nverty;
    double nodx[12], nody[12];
    double bi[16], br[16];

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (ix = 0; ix < nnodes; ix++)
        for (iy = 0; iy < nnodes; iy++)
            Br[ix+nnodes*iy] = Bi[ix+nnodes*iy] = 0.0;

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
            dint_tri_lin_sing(&g3[0], nodx, &accelerators[ex], &g4[0], k, br, bi);
            break;
        case 4:
            dint_quad_lin_sing(&g4[0], nodx, &accelerators[ex], k, &br, &bi);
            break;
        }
        for (ix = 0; ix < nvertx; ix++)
        {
            for (iy = 0; iy < nvertx; iy++)
            {
                Br[elemx[ix]+nnodes*elemx[iy]] += br[ix+nvertx*iy];
                Bi[elemx[ix]+nnodes*elemx[iy]] += bi[ix+nvertx*iy];
            }
        }

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
                    dint_quadquad_lin(&g4[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, br, bi);
                    break;
                case 3:
                    dint_triquad_lin(&g3[gs], nody, &accelerators[ey], &g4[gs], nodx, &accelerators[ex], k, br, bi);
                    break;
                }
                break;
            case 3:
                switch (nverty)
                {
                case 4:
                    dint_triquad_lin(&g3[gs], nodx, &accelerators[ex], &g4[gs], nody, &accelerators[ey], k, br, bi);
                    break;
                case 3:
                    dint_tritri_lin(&g3[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, br, bi);
                    break;
                }
            }

            for (ix = 0; ix < nvertx; ix++)
            {
                for (iy = 0; iy < nverty; iy++)
                {
                    Br[elemx[ix]+nnodes*elemy[iy]] += br[ix+nvertx*iy];
                    Br[elemy[iy]+nnodes*elemx[ix]] += br[ix+nvertx*iy];
                    Bi[elemx[ix]+nnodes*elemy[iy]] += bi[ix+nvertx*iy];
                    Bi[elemy[iy]+nnodes*elemx[ix]] += bi[ix+nvertx*iy];
                }
            }
        }
    }
    free(accelerators);
}

#include "ibemD_li.h"

#include "integral_indirect.h"
#include "vector.h"
#include "mesh.h"
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------ */
/* Compute surface system B matrices of an indirect bem model               */
void iBEM_D_surf_lin(int nnodes,
                     const double *nodes,
                     int nelements,
                     const double *elements,
                     const gauss_t *g3,
                     const gauss_t *g4,
                     const double *dist,
                     double k,
                     double *Dr,
                     double *Di)
{
    accelerator_t *accelerators;
    int j, ey, s, ex, gs, ix, iy;
    int elemx[4], elemy[4], nvertx, nverty;
    double nodx[12], nody[12];
    double dr[16], di[16];

    /* Compute element centres */
    accelerators = (accelerator_t *)calloc(nelements, sizeof(accelerator_t));
    init_accelerators(nnodes, nodes, nelements, elements, accelerators);

    /* Clear output matrices */
    for (ix = 0; ix < nnodes; ix++)
        for (iy = 0; iy < nnodes; iy++)
            Dr[ix+nnodes*iy] = Di[ix+nnodes*iy] = 0.0;

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
            dintD_tri_lin_sing(&g3[0], nodx, &accelerators[ex], &g4[0], k, dr, di);
            break;
        case 4:
            dintD_quad_lin_sing(&g4[0], nodx, &accelerators[ex], k, dr, di);
            break;
        }
        for (ix = 0; ix < nvertx; ix++)
        {
            for (iy = 0; iy < nvertx; iy++)
            {
                Dr[elemx[ix]+nnodes*elemx[iy]] += dr[ix+nvertx*iy];
                Di[elemx[ix]+nnodes*elemx[iy]] += di[ix+nvertx*iy];
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
                    dintD_quadquad_lin(&g4[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, dr, di);
                    break;
                case 3:
                    dintD_triquad_lin(&g3[gs], nody, &accelerators[ey], &g4[gs], nodx, &accelerators[ex], k, dr, di);
                    break;
                }
                break;
            case 3:
                switch (nverty)
                {
                case 4:
                    dintD_triquad_lin(&g3[gs], nodx, &accelerators[ex], &g4[gs], nody, &accelerators[ey], k, dr, di);
                    break;
                case 3:
                    dintD_tritri_lin_fast(&g3[gs], nodx, &accelerators[ex], nody, &accelerators[ey], k, dr, di);
                    break;
                }
            }

            for (ix = 0; ix < nvertx; ix++)
            {
                for (iy = 0; iy < nverty; iy++)
                {
                    Dr[elemx[ix]+nnodes*elemy[iy]] += dr[ix+nvertx*iy];
                    Dr[elemy[iy]+nnodes*elemx[ix]] += dr[ix+nvertx*iy];
                    Di[elemx[ix]+nnodes*elemy[iy]] += di[ix+nvertx*iy];
                    Di[elemy[iy]+nnodes*elemx[ix]] += di[ix+nvertx*iy];
                }
            }
        }
    }
    free(accelerators);
}

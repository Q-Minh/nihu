/*bemHG_const  Generate acoustic BEM system matrices
 *   B = ibemB_const(nodes, elements, g3, g4, dist, k, points)
 *   Computes the acoustic BEM system matrix B for the case of constant
 %   shape functions.
 * Parameters:
 *   nodes    : N x 3 matrix of xyz coordinates of model vertices
 *   elements : matrix of element node indices. Each row describes one
 *              element in the form [3 n1 n2 n3 0] for TRIA and
 *              [4 n1 n2 n3 n4] for QUAD elements.
 *   g3, g4   : Gaussian quadrature structure for TRIA and QUAD elements
 *              obtained from GAUSSQUAD2
 *   dist     : Distance limits governing integration density parameter
 *   k        : acoustic wave number
 *   points   : optional, if not given, then the acoustuc surface matrices
 *              are computed. If defined, then the field point matrices
 *              are returned. Points is a M x 3 matrix containing the xyz
 *              coordinates of the field points
 *
 * Peter Fiala
 * 2011
 */

/* $Make: mex -O -output ibemD_const ibemD_const.mex.c ibemD_con.c mesh.c integral_indirect.c quadrature.c element.c vector.c green.c $ */

#include "mex.h"
#include "ibemD_con.h"
#include "types.h"

/* ------------------------------------------------------------------------ */
/* The entry point of the mex function. */
/* This function reads Matlab parameters, allocates space for the output */
/* and calls C routines that perform the computations */
void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    /* input parameters */
    double *nodes, *elements, *points, *dist, k;
    int nnodes, nelements, npoints;
    gauss_t *g3, *g4;
    /* output parameters */
    double *Br, *Bi;
    /* local variables */
    int j;

    /* Allocate memory for Gaussian integration structures */
    g3 = (gauss_t *)calloc(4, sizeof(gauss_t));
    g4 = (gauss_t *)calloc(4, sizeof(gauss_t));

    /* transfer input parameters */
    nodes = mxGetPr(prhs[0]);
    nnodes = mxGetM(prhs[0]);
    elements = mxGetPr(prhs[1]);
    nelements = mxGetM(prhs[1]);
    for (j = 0; j < 4; j++)
    {
        g3[j].num = mxGetM(mxGetField(prhs[2], j, "N"));
        g3[j].N = mxGetPr(mxGetField(prhs[2], j, "N"));
        g3[j].Nxi = mxGetPr(mxGetField(prhs[2], j, "Nxi"));
        g3[j].Neta = mxGetPr(mxGetField(prhs[2], j, "Neta"));
        g3[j].w = mxGetPr(mxGetField(prhs[2], j, "w"));
        g3[j].xi = mxGetPr(mxGetField(prhs[2], j, "xi"));
    }
    for (j = 0; j < 4; j++)
    {
        g4[j].num = mxGetM(mxGetField(prhs[3], j, "N"));
        g4[j].N = mxGetPr(mxGetField(prhs[3], j, "N"));
        g4[j].Nxi = mxGetPr(mxGetField(prhs[3], j, "Nxi"));
        g4[j].Neta = mxGetPr(mxGetField(prhs[3], j, "Neta"));
        g4[j].w = mxGetPr(mxGetField(prhs[3], j, "w"));
        g4[j].xi = mxGetPr(mxGetField(prhs[3], j, "xi"));
    }
    dist = mxGetPr(prhs[4]);
    k = mxGetScalar(prhs[5]);

    if (nrhs == 6)
    {
        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(nelements, nelements, mxCOMPLEX);
        Br = mxGetPr(plhs[0]);
        Bi = mxGetPi(plhs[0]);

        /* call C subroutine */
        iBEM_D_surf_const(nnodes, nodes, nelements, elements,
                          g3, g4, dist, k, Br, Bi);
    }
    else if (nrhs == 7)
    {
        points = mxGetPr(prhs[6]);
        npoints = mxGetM(prhs[6]);

        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(npoints, nelements, mxCOMPLEX);
        Br = mxGetPr(plhs[0]);
        Bi = mxGetPi(plhs[0]);

        /* call C subroutine */
        /* matrix_field_const(nnodes, nodes, nelements, elements,
                npoints, points, g3, g4, dist, k, Br, Bi); */
    }

    free(g3);
    free(g4);
}
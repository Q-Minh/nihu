/* $Make: mex -O -output bemHG_const_sp bemHG_const_sp.mex.c bemHG_con.c mesh.c integral_direct.c quadrature.c element.c vector.c green.c$ */

#include "mex.h"
#include "bemHG_con.h"
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
    /* Allocate memory for Gaussian integration structures */
    gauss_t *g3 = new gauss_t[4];
    gauss_t *g4 = new gauss_t[4];

    /* transfer input parameters */
    double *nodes = mxGetPr(prhs[0]);
    int nnodes = mxGetM(prhs[0]);
    double *elements = mxGetPr(prhs[1]);
    int nelements = mxGetM(prhs[1]);
    for (int j = 0; j < 4; j++)
    {
        g3[j].num = mxGetM(mxGetField(prhs[2], j, "N"));
        g3[j].N = mxGetPr(mxGetField(prhs[2], j, "N"));
        g3[j].w = mxGetPr(mxGetField(prhs[2], j, "w"));

        g4[j].num = mxGetM(mxGetField(prhs[3], j, "N"));
        g4[j].N = mxGetPr(mxGetField(prhs[3], j, "N"));
        g4[j].Nxi = mxGetPr(mxGetField(prhs[3], j, "Nxi"));
        g4[j].Neta = mxGetPr(mxGetField(prhs[3], j, "Neta"));
        g4[j].w = mxGetPr(mxGetField(prhs[3], j, "w"));
    }
    double *dist = mxGetPr(prhs[4]);
    complex_scalar k(*mxGetPr(prhs[5]), *mxGetPi(prhs[5]));

    if (nrhs == 7)
    {
        double *pairs = mxGetPr(prhs[6]);
        int npairs = mxGetM(prhs[6]);

        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
        double *Ar = mxGetPr(plhs[0]);
        double *Ai = mxGetPi(plhs[0]);
        double *Br = mxGetPr(plhs[1]);
        double *Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        if (k.imag()==0)
            matrix_surf_const_sparse(nnodes, nodes, nelements, elements,
                                     npairs, pairs, g3, g4, dist, k.real(), Ar, Ai, Br, Bi);
        else
            matrix_surf_const_sparse(nnodes, nodes, nelements, elements,
                                     npairs, pairs, g3, g4, dist, k, Ar, Ai, Br, Bi);
    }
    else if (nrhs == 8)
    {
        double *points = mxGetPr(prhs[6]);
        int npoints = mxGetM(prhs[6]);

        double *pairs = mxGetPr(prhs[7]);
        int npairs = mxGetM(prhs[7]);

        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(npairs, 1, mxCOMPLEX);
        double *Ar = mxGetPr(plhs[0]);
        double *Ai = mxGetPi(plhs[0]);
        double *Br = mxGetPr(plhs[1]);
        double *Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        if (k.imag()==0)
            matrix_field_const_sparse(nnodes, nodes, nelements, elements,
                                      npoints, points, npairs, pairs, g3, g4, dist, k.real(), Ar, Ai, Br, Bi);
        else
            matrix_field_const_sparse(nnodes, nodes, nelements, elements,
                                      npoints, points, npairs, pairs, g3, g4, dist, k, Ar, Ai, Br, Bi);
    }

    delete [] g3;
    delete [] g4;

    return;
}

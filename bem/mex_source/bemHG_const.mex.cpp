/* $Make: mex -O -output bemHG_const bemHG_const.mex.cpp mesh.cpp element.cpp vector.cpp$ */

#include "mex.h"
#include "bemHG_con.hpp"
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
        g3[j].Nxi = mxGetPr(mxGetField(prhs[2], j, "Nxi"));
        g3[j].Neta = mxGetPr(mxGetField(prhs[2], j, "Neta"));
        g3[j].w = mxGetPr(mxGetField(prhs[2], j, "w"));

        g4[j].num = mxGetM(mxGetField(prhs[3], j, "N"));
        g4[j].N = mxGetPr(mxGetField(prhs[3], j, "N"));
        g4[j].Nxi = mxGetPr(mxGetField(prhs[3], j, "Nxi"));
        g4[j].Neta = mxGetPr(mxGetField(prhs[3], j, "Neta"));
        g4[j].w = mxGetPr(mxGetField(prhs[3], j, "w"));
    }
    double *dist = mxGetPr(prhs[4]);
    complex_scalar k(*mxGetPr(prhs[5]), *mxGetPi(prhs[5]));

    if (nrhs == 6)
    {
        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(nelements, nelements, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(nelements, nelements, mxCOMPLEX);
        double *Ar = mxGetPr(plhs[0]);
        double *Ai = mxGetPi(plhs[0]);
        double *Br = mxGetPr(plhs[1]);
        double *Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        if (k.imag() == 0.0)
            matrix_surf_const(nnodes, nodes, nelements, elements,
                              g3, g4, dist, k.real(), Ar, Ai, Br, Bi);
        else
            matrix_surf_const(nnodes, nodes, nelements, elements,
                              g3, g4, dist, k, Ar, Ai, Br, Bi);
    }
    else if (nrhs == 7)
    {
        double *points = mxGetPr(prhs[6]);
        int npoints = mxGetM(prhs[6]);

        /* Allocate output parameters */
        plhs[0] = mxCreateDoubleMatrix(npoints, nelements, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(npoints, nelements, mxCOMPLEX);
        double *Ar = mxGetPr(plhs[0]);
        double *Ai = mxGetPi(plhs[0]);
        double *Br = mxGetPr(plhs[1]);
        double *Bi = mxGetPi(plhs[1]);

        /* call C subroutine */
        if (k.imag() == 0.0)
            matrix_field_const(nnodes, nodes, nelements, elements,
                               npoints, points, g3, g4, dist, k.real(), Ar, Ai, Br, Bi);
        else
            matrix_field_const(nnodes, nodes, nelements, elements,
                               npoints, points, g3, g4, dist, k, Ar, Ai, Br, Bi);
    }

    delete [] g3;
    delete [] g4;
    return;
}

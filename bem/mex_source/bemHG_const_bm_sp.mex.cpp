/* $Make: mex -O -output bemHG_const_bm_sp bemHG_const_bm_sp.mex.cpp mesh.cpp element.cpp vector.cpp$ */

#include "mex.h"
#include "bemHG_con_bm.hpp"
#include "types.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum {NGAUSS = 4};

    /* Allocate memory for Gaussian integration structures */
    gauss_t *g3 = new gauss_t[NGAUSS];
    gauss_t *g4 = new gauss_t[NGAUSS];

    /* transfer input parameters */
    double *nodes = mxGetPr(prhs[0]);
    int nnodes = mxGetM(prhs[0]);
    double *elements = mxGetPr(prhs[1]);
    int nelements = mxGetM(prhs[1]);
    for (int j = 0; j < NGAUSS; j++)
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
    complex_scalar alpha(*mxGetPr(prhs[6]), *mxGetPi(prhs[6]));

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
    if (k.imag() == 0.0)
        matrix_surf_const_bm_sparse(nnodes, nodes, nelements, elements,
                                    npairs, pairs, g3, g4, dist, k.real(), alpha, Ar, Ai, Br, Bi);
    else
        matrix_surf_const_bm_sparse(nnodes, nodes, nelements, elements,
                                    npairs, pairs, g3, g4, dist, k, alpha, Ar, Ai, Br, Bi);

    delete [] g3;
    delete [] g4;

    return;
}


% $Make: mex -O -output Boonen Boonen13.mex.cpp$ */

#include <mex.h>

#include "bem.hpp"

typedef tmp::vector<quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef mesh_t::x_t x_t;
typedef bem<elem_type_vector, isoparametric_field> bem_t;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *nodes = mxGetPr(prhs[0]);
unsigned nnodes = mxGetM(prhs[0]);
double *elements = mxGetPr(prhs[1]);
unsigned nelements = mxGetM(prhs[1]);

%
CMesh<double, 3> Mesh;
Mesh.buildFromMex(nnodes, nodes, nelements, elements);

double k = mxGetScalar(prhs[2]);
double precision = mxGetScalar(prhs[4]);

BEM_Type Bem;
Bem.buildField(CField0<double, Result_Type, Descriptor_Type >::ISOPARAM, Mesh);
Bem.setWaveNumber(k);
Bem.setPrecision(precision);

double *points = mxGetPr(prhs[3]);
unsigned npoints = mxGetM(prhs[3]);

plhs[0] = mxCreateDoubleMatrix(nnodes, nnodes, mxCOMPLEX);
double *Gr = mxGetPr(plhs[0]);
double *Gi = mxGetPi(plhs[0]);

for (unsigned iNode = 0; iNode < nnodes; ++iNode)
{
Matrix<Result_Type, Dynamic, 1> G = Bem.surfaceVector(iNode, Vector3d(points[iNode], points[iNode+nnodes], points[iNode+2*nnodes]));
for (unsigned i = 0; i < nnodes; ++i)
{
Gr[iNode+i*nnodes] = G[i].real();
Gi[iNode+i*nnodes] = G[i].imag();
}
}

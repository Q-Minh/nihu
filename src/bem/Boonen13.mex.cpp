#include <mex.h>
#include "bem.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef mesh_t::x_t x_t;
typedef bem<elem_type_vector, isoparametric_field> bem_t;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *nodes = mxGetPr(prhs[0]);
    unsigned nnodes = mxGetM(prhs[0]);
    double *elements = mxGetPr(prhs[1]);
    unsigned nelements = mxGetM(prhs[1]);

	/* build mesh */
	mesh_t mesh;
	mesh.build_from_mex<4+1>(nodes, nnodes, elements, nelements);
	bem_t b(mesh);

    double k = mxGetScalar(prhs[2]);

	double *points = mxGetPr(prhs[3]);
    unsigned npoints = mxGetM(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(npoints, nnodes, mxCOMPLEX);
    double *Gr = mxGetPr(plhs[0]);
    double *Gi = mxGetPi(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(npoints, nnodes, mxCOMPLEX);
    double *Hr = mxGetPr(plhs[1]);
    double *Hi = mxGetPi(plhs[1]);

	for (unsigned iPoint = 0; iPoint < npoints; ++iPoint)
	{
		x_t x0;
		x0 << points[iPoint], points[iPoint+npoints], points[iPoint+2*npoints];
		auto result = b.eval(x0, k);

		for (unsigned i = 0; i < nnodes; ++i)
		{
			Gr[iPoint+i*npoints] = result(i,0).real();
			Gi[iPoint+i*npoints] = result(i,0).imag();
			Hr[iPoint+i*npoints] = result(i,1).real();
			Hi[iPoint+i*npoints] = result(i,1).imag();
		}
	}
}


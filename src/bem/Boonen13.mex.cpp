#include <mex.h>
#include "bem.hpp"
#include "mex_matrix.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef green_HG_kernel kernel_t;
typedef bem<elem_type_vector, isoparametric_field, kernel_t> bem_t;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *nodes = mxGetPr(prhs[0]);
    size_t nnodes = mxGetM(prhs[0]);
    double *elements = mxGetPr(prhs[1]);
    size_t nelements = mxGetM(prhs[1]);

	mesh_t mesh;
	mesh.build_from_mex<4+1>(nodes, nnodes, elements, nelements);

	bem_t b(mesh);

    double k = mxGetScalar(prhs[2]);

	double *points = mxGetPr(prhs[3]);
    size_t npoints = mxGetM(prhs[3]);
	
	mex_complex_matrix<kernel_t::num_elements> output(npoints, nnodes);
    for (size_t i = 0; i < kernel_t::num_elements; ++i)
        plhs[i] = output.get_matrix(i);

	unsigned prevproc = 0;
	for (unsigned iPoint = 0; iPoint < npoints; ++iPoint)
	{
		typedef mesh_t::x_t x_t;
		x_t x0;

		x0 << points[iPoint], points[iPoint+npoints], points[iPoint+2*npoints];
		auto result = b.eval(x0, k);

		for (unsigned i = 0; i < nnodes; ++i)
			output(iPoint, i) += result.row(i);

		unsigned proc = 80*iPoint/npoints;
		if (proc > prevproc)
		{
			mexPrintf("-", proc);
			mexEvalString("drawnow;");
			prevproc = proc;
		}
	}
}


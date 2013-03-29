#include <mex.h>
#include "weighted_residual.hpp"
#include "mex_matrix.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef function_space<mesh_t, isoparametric_field> test_space_t;
typedef function_space<mesh_t, isoparametric_field> trial_space_t;
typedef green_G_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> weighted_residual_t;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mex_real_matrix nodes(prhs[0]);
    mex_real_matrix elements(prhs[1]);

	mesh_t mesh;
	mesh.build_from_mex<4+1>(nodes, elements);

    double k = mxGetScalar(prhs[2]);
	kernel_t::set_wave_number(k);

	mex_complex_matrix output(nodes.rows(), nodes.rows());
    plhs[0] = output.get_matrix();

	test_space_t test_space(mesh);
	trial_space_t trial_space(mesh);
	weighted_residual<kernel_t, test_space_t, trial_space_t> wr(test_space, trial_space);
	wr.eval(output);
}


#include "weighted_residual.hpp"
#include "../util/mex_matrix.h"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;
typedef function_space<mesh_t, constant_field, dirac_field> test_space_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;
typedef helmholtz_HG_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> wr_t;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // read mesh from Matlab and build
    mex::real_matrix nodes(prhs[0]);
    mex::real_matrix elements(prhs[1]);
	mesh_t mesh;
	mesh.build_from_mex<4+1>(nodes, elements);

    // initialise kernel
    double k = mxGetScalar(prhs[2]);
    kernel_t::set_wave_number(k);
    
    // initialise function spaces
    test_space_t test(mesh);
    trial_space_t trial(mesh);
    
    // allocate space for output matrix
    mex::complex_matrix G(test.get_num_dofs(), trial.get_num_dofs(), plhs[0]);
    mex::complex_matrix H(test.get_num_dofs(), trial.get_num_dofs(), plhs[1]);
    couple<mex::complex_matrix> result(G, H);

    // evaluate weighted residual into result
    wr_t wr(test, trial);
    wr.eval(result);
}


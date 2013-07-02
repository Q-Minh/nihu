#include "weighted_residual.hpp"
#include "../util/mex_matrix.h"
#include "../library/helmholtz_kernel.hpp"

// #define GALERKIN_ISO

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef mesh<elem_type_vector> mesh_t;

typedef helmholtz_GH_kernel kernel_t;

#ifdef COLLOC_CONSTANT
typedef function_space_view<mesh_t, constant_field> test_space_t;
typedef test_space_t trial_space_t;
typedef weighted_residual<true, kernel_t, test_space_t, trial_space_t> wr_t;
#endif

#ifdef GALERKIN_CONSTANT
typedef function_space_view<mesh_t, constant_field> test_space_t;
typedef test_space_t trial_space_t;
typedef weighted_residual<false, kernel_t, test_space_t, trial_space_t> wr_t;
#endif

#ifdef GALERKIN_ISO
typedef function_space_view<mesh_t, isoparametric_field> test_space_t;
typedef test_space_t trial_space_t;
typedef weighted_residual<false, kernel_t, test_space_t, trial_space_t> wr_t;
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // read mesh from Matlab and build
    mex::real_matrix nodes(prhs[0]);
    mex::real_matrix elements(prhs[1]);
	mesh_t msh(nodes, elements);

    // initialise kernel
    kernel_t kernel;
    double k = mxGetScalar(prhs[2]);
    kernel.set_wave_number(k);
    
    // initialise function spaces
    test_space_t test(msh);
    trial_space_t trial(msh);
    
    // allocate space for output matrix
    mex::complex_matrix G(test.get_num_dofs(), trial.get_num_dofs(), plhs[0]);
    mex::complex_matrix H(test.get_num_dofs(), trial.get_num_dofs(), plhs[1]);
    couple<mex::complex_matrix> result(G, H);
    
    mex::real_matrix n_kernel_eval(1, 1, plhs[2]);

    // evaluate weighted residual into result
    wr_t::eval(result, kernel, test, trial);
    
    // export number of kernel evaluations
    n_kernel_eval(0,0) = kernel.get_num_evaluations();
}


#include "weighted_residual.hpp"
#include "../util/mex_matrix.h"

// #define GALERKIN_ISO

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef Mesh<elem_type_vector> mesh_t;

#ifdef COLLOC_CONSTANT
typedef function_space<mesh_t, constant_field, dirac_field> test_space_t;
typedef function_space<mesh_t, constant_field, function_field> trial_space_t;
#endif

#ifdef GALERKIN_CONSTANT
typedef function_space<mesh_t, constant_field, function_field> test_space_t;
typedef test_space_t trial_space_t;
#endif

#ifdef GALERKIN_ISO
typedef function_space<mesh_t, isoparametric_field, function_field> test_space_t;
typedef test_space_t trial_space_t;
#endif

typedef unit_kernel kernel_t;
typedef weighted_residual<kernel_t, test_space_t, trial_space_t> wr_t;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // read mesh from Matlab and build
    mex::real_matrix nodes(prhs[0]);
    mex::real_matrix elements(prhs[1]);
	mesh_t mesh(nodes, elements);

    // initialise function spaces
    test_space_t test(mesh);
    trial_space_t trial(mesh);
    
    // allocate space for output matrix
    mex::real_matrix G(test.get_num_dofs(), trial.get_num_dofs(), plhs[0]);
    
    mex::real_matrix n_kernel_eval(1, 1, plhs[1]);

    // evaluate weighted residual into result
    wr_t wr(test, trial);
    wr.eval(G);
    
    // export number of kernel evaluations
    n_kernel_eval(0,0) = kernel_t::get_num_evaluations();
}


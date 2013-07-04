#include "integral_operator.hpp"
#include "../util/mex_matrix.h"
#include "../library/helmholtz_kernel.hpp"

typedef tmp::vector<tria_1_elem, quad_1_elem> elem_type_vector;
typedef mesh<elem_type_vector> mesh_t;

#ifdef COLLOCATIONAL
    typedef formalism::collocational FORMALISM;
#else
    typedef formalism::general FORMALISM;
#endif

#ifdef CONSTANT
    typedef field_option::constant FIELD_OPTION;
#else
    typedef field_option::isoparametric FIELD_OPTION;
#endif
    
typedef function_space_view<mesh_t, FIELD_OPTION> func_sp_t;

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    // read mesh from Matlab and build
    mex::real_matrix nodes(prhs[0]);
    mex::real_matrix elements(prhs[1]);
	mesh_t msh(nodes, elements);

    // initialise function spaces
    auto func_sp = create_function_space_view(msh, field_option::constant());
    
    // allocate space for output matrix
    int n = func_sp.get_num_dofs();
    mex::complex_matrix G(n, n, plhs[0]);
    mex::complex_matrix H(n, n, plhs[1]);
    auto result = create_couple(G, H);
    
    // initialise kernel and boundary operator
    double k = mxGetScalar(prhs[2]);
    helmholtz_GH_kernel kernel;
    kernel.set_wave_number(k);
    auto b_op = create_integral_operator(kernel);
    
    // evaluate weighted residual into result
    result += b_op(func_sp, func_sp, FORMALISM());
    
    // export number of kernel evaluations
    mex::real_matrix n_kernel_eval(1, 1, plhs[2]);
    n_kernel_eval(0,0) = b_op.get_kernel().get_num_evaluations();
}

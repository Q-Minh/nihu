#include "util/mex_matrix.hpp"
#include "core/mesh.hpp"
#include "library/lib_element.hpp"
#include "library/laplace_kernel.hpp"
#include "core/weighted_residual.hpp"

typedef NiHu::mex::real_matrix<double> dMatrix;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
// D = domain_integral_test(nodes, delements, belements)
{
	dMatrix nodes(prhs[0]), dElements(prhs[1]), bElements(prhs[2]);
	auto domain = NiHu::create_mesh(nodes, dElements, NiHu::tria_1_volume_tag());
	auto boundary = NiHu::create_mesh(nodes, bElements, NiHu::line_1_tag());

	auto const &dspace = NiHu::constant_view(domain);
	auto const &bspace = NiHu::constant_view(boundary);

	auto K = NiHu::create_integral_operator(NiHu::laplace_2d_SLP_kernel());

	dMatrix D(boundary.get_num_elements(), domain.get_num_elements(), plhs[0]);

	D << NiHu::dirac(bspace) * K[dspace];
}


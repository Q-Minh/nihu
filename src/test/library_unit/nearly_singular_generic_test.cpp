#include "library/lib_element.hpp"
#include "library/nearly_singular_collocational.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"

int main(void)
{
	typedef NiHu::tria_2_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 2.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 2.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
	typedef NiHu::tria_0_shape_set nset_t;
	typedef NiHu::field<elem_t, nset_t, NiHu::_1d> trial_field_t;
	trial_field_t trial_field(elem, trial_field_t::dofs_t::Zero());

	typedef NiHu::laplace_3d_SLP_kernel kernel_t;
	kernel_t kernel;

	typedef NiHu::nearly_singular_collocational<trial_field_t, kernel_t, 15U, 15U> nsg_t;
	nsg_t nsg(trial_field, kernel);

	Eigen::Matrix<double, 1, 1> integral;
	integral.setZero();

	elem_t::x_t x;
	x << .2, .3, .5;

	kernel_t::test_input_t tsi(x);

	nsg.integrate(integral, tsi);

	std::cout << integral << std::endl;

	return 0;
}

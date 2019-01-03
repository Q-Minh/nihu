#include "library/lib_element.hpp"
#include "library/nearly_singular_generic.hpp"
#include "library/laplace_kernel.hpp"

int main(void)
{
	typedef NiHu::tria_2_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 2.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 2.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
	typedef NiHu::tria_1_shape_set nset_t;
	typedef NiHu::field<elem_t, nset_t, NiHu::_1d> field_t;

	typedef NiHu::laplace_3d_SLP_kernel kernel_t;
	kernel_t kernel;

	typedef NiHu::nearly_singular_generic<field_t, kernel_t, 5, 5> nsg_t;
	nsg_t nsg(elem, kernel);

	Eigen::Matrix<double, 1, 3> integral;
	integral.setZero();

	elem_t::x_t x;
	x << .2, .3, .5;

	kernel_t::test_input_t tsi(x);

	nsg.integrate(integral, tsi);

	return 0;
}

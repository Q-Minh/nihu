#include "library/lib_element.hpp"
#include "library/nearly_singular_collocational.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"

typedef NiHu::helmholtz_3d_SLP_kernel<double> kernel_t;
double wave_number = 1;
kernel_t kernel(wave_number);
typedef kernel_t::result_t Scalar_t;

void test_quadratic()
{
	typedef NiHu::tria_2_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 2.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 2.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
	typedef NiHu::tria_2_shape_set nset_t;
	typedef NiHu::field<elem_t, nset_t, NiHu::_1d> trial_field_t;
	trial_field_t trial_field(elem, trial_field_t::dofs_t::Zero());

	typedef NiHu::nearly_singular_collocational<trial_field_t, kernel_t, 15U, 15U> nsg_t;
	nsg_t nsg(trial_field, kernel);

	Eigen::Matrix<Scalar_t, 1, nset_t::num_nodes> integral;
	integral.setZero();

	elem_t::x_t x;
	x << .2, .3, .5;

	kernel_t::test_input_t tsi(x);

	nsg.integrate(integral, tsi);

	std::cout << integral << std::endl;
	std::cout << integral.sum() << std::endl;
}


void test_linear()
{
	typedef NiHu::tria_1_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 2.0, 0.0,
		0.0, 0.0, 2.0,
		0.0, 0.0, 0.0;
	elem_t elem(coords);
	typedef NiHu::tria_1_shape_set nset_t;
	typedef NiHu::field<elem_t, nset_t, NiHu::_1d> trial_field_t;
	trial_field_t trial_field(elem, trial_field_t::dofs_t::Zero());

	typedef NiHu::nearly_singular_collocational<trial_field_t, kernel_t, 15U, 15U> nsg_t;
	nsg_t nsg(trial_field, kernel);

	Eigen::Matrix<Scalar_t, 1, nset_t::num_nodes> integral;
	integral.setZero();

	elem_t::x_t x;
	x << .2, .3, .5;

	kernel_t::test_input_t tsi(x);

	nsg.integrate(integral, tsi);

	std::cout << integral << std::endl;
	std::cout << integral.sum() << std::endl;

	std::cout << NiHu::helmholtz_3d_SLP_collocation_constant_plane_nearly_singular::eval(
		elem, x, wave_number
	) << std::endl;
}


int main(void)
{
	std::cout << "quadratic" << std::endl;
	test_quadratic();
	std::cout << std::endl;
	std::cout << "linear" << std::endl;
	test_linear();

	return 0;
}

#include "library/lib_element.hpp"
#include "library/nearly_singular_collocational.hpp"
#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/lib_shape.hpp"
#include "core/double_integral.hpp"

typedef NiHu::helmholtz_3d_DLP_kernel<double> kernel_t;
double wave_number = 1;
kernel_t kernel(wave_number);
typedef kernel_t::result_t Scalar_t;
typedef NiHu::quad_1_gauss_shape_set nset_t;


void test_fields()
{
	typedef NiHu::quad_1_elem elem_t;
	elem_t::coords_t trial_coords;
	trial_coords <<
		0.0, 2.0, 2.0, 0.0,
		0.0, 0.0, 2.0, 2.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t trial_elem(trial_coords);
	typedef NiHu::field<elem_t, NiHu::quad_1_gauss_shape_set> trial_field_t;
	trial_field_t trial_field(trial_elem, trial_field_t::dofs_t::Constant(0));

	elem_t::coords_t test_coords;
	test_coords <<
		0.0, 0.0, 0.0, 0.0,
		0.0, 2.0, 2.0, 0.0,
		0.0, 0.0, 2.0, 2.0;
	elem_t test_elem(test_coords);
	typedef NiHu::field<elem_t, NiHu::quad_1_gauss_shape_set> test_field_orig_t;
	test_field_orig_t test_field_orig(test_elem, test_field_orig_t::dofs_t::Constant(10));
	typedef NiHu::dirac_field<test_field_orig_t> test_field_t;
	test_field_t test_field = NiHu::dirac(test_field_orig);

	typedef Eigen::Matrix<Scalar_t, 4, 4> result_t;

	result_t result;
	result.setZero();
	NiHu::nearly_singular_integral<kernel_t, test_field_t, trial_field_t>::eval(
		result, kernel, test_field, trial_field
	);
	std::cout << result << std::endl;

	std::cout << std::endl;


	result_t res2;
	res2.setZero();

	for (unsigned i = 0; i < 4; ++i)
	{
		typedef NiHu::nearly_singular_collocational<trial_field_t, kernel_t, 15U, 15U> nsg_t;
		nsg_t nsg(trial_field, kernel);

		elem_t::x_t x = test_elem.get_x(test_field_t::nset_t::corner_at(i));
		kernel_t::test_input_t tsi(x);
		nsg.integrate(res2.row(i), tsi);
	}

	std::cout << res2 << std::endl;

	std::cout << std::endl;

	result_t res3;
	res3.setZero();
	for (unsigned i = 0; i < 4; ++i)
	{
		elem_t::x_t x = test_elem.get_x(test_field_t::nset_t::corner_at(i));
		res3(i,0) = NiHu::helmholtz_3d_DLP_collocation_constant_plane_nearly_singular::eval(
			trial_elem, x, wave_number
		);
	}

	std::cout << res3 << std::endl;
	std::cout << res2.rowwise().sum() << std::endl;

	std::cout << std::endl;

	result_t result4;
	result4 = NiHu::double_integral<kernel_t, test_field_t, trial_field_t>::eval(
		kernel, test_field, trial_field, std::true_type()
	);
	std::cout << result4 << std::endl;
}



#if 0
void test_quadratic()
{
	typedef NiHu::quad_2_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 1.0,
		0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
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
	typedef NiHu::quad_1_elem elem_t;
	elem_t::coords_t coords;
	coords <<
		0.0, 2.0, 2.0, 0.0,
		0.0, 0.0, 2.0, 2.0,
		0.0, 0.0, 0.0, 0.0;
	elem_t elem(coords);
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

#endif


int main(void)
{
#if 0
	std::cout << "quadratic" << std::endl;
	test_quadratic();
	std::cout << std::endl;
	std::cout << "linear" << std::endl;
	test_linear();
#endif

	test_fields();

	return 0;
}

#include <core/weighted_residual.hpp>
#include <library/laplace_kernel.hpp>
#include <library/laplace_singular_integrals.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;

typedef laplace_2d_SLP_kernel kernel_t;
typedef line_1_elem elem_t;
typedef elem_t::coords_t coords_t;

auto kernel = kernel_t();

void constant_collocation_tester(void)
{
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef dirac_field<trial_field_t> test_field_t;
	typedef double_integral<kernel_t, test_field_t, trial_field_t> integral_t;

	coords_t n1;
	n1 <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(n1);

	trial_field_t const &trial_field = create_field_view(elem, field_option::constant());
	test_field_t const &test_field = dirac(trial_field);

	std::cout << integral_t::eval(kernel, test_field, trial_field, std::true_type()) << std::endl;
}

void galerkin_face_tester(void)
{
	typedef field_view<elem_t, field_option::constant> trial_field_t;
	typedef trial_field_t test_field_t;

	coords_t n1;
	n1 <<
		0.0, 1.0,
		0.0, 0.0;
	elem_t elem(n1);

	trial_field_t const &trial_field = create_field_view(elem, field_option::constant());
	test_field_t const &test_field = trial_field;

	typedef double_integral<kernel_t, test_field_t, trial_field_t> integral_t;
	std::cout << integral_t::eval(kernel, test_field, trial_field, std::true_type()) << std::endl;
}

int main(void)
{
	// constant_collocation_tester();

	galerkin_face_tester();

	return 0;
}

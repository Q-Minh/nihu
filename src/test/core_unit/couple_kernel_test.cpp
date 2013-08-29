#include <iostream>
#include "bem/kernel.hpp"
#include "library/laplace_kernel.hpp"

int main(void)
{
	// create a set of kernels
	auto K = create_couple_kernel(
		laplace_3d_SLP_kernel(),
		laplace_3d_DLP_kernel(),
		laplace_3d_SLP_kernel()
	);

	// create an element
	typename quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;
	quad_1_elem elem(coords);

	// create test and trial kernel inputs
	auto xi1 = quad_domain::get_center();
	auto xi2 = quad_domain::get_corner(0);
	kernel_traits<decltype(K)>::test_input_t tst_i(elem, xi1);
	kernel_traits<decltype(K)>::trial_input_t trl_i(elem, xi2);

	// evaluate kernel
	auto res = K(tst_i, trl_i);

	// weight kernel by a matrix
	double a = 3.0;
	auto wres = a * res * 2;

	std::cout << wres.get<0>() << std::endl;
	std::cout << wres.get<1>() << std::endl;
	std::cout << wres.get<2>() << std::endl;

	return 0;
}


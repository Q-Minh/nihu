#include "library/laplace_kernel.hpp"
#include <iostream>
#include <chrono>

unsigned const N = 1e6;

template <class kernel>
void tester(kernel_base<kernel> const &k)
{
	quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;

	quad_1_elem elem(coords);
	auto xi1 = quad_1_elem::domain_t::get_corner(0);
	auto xi2 = quad_1_elem::domain_t::get_corner(1);

	auto start = std::chrono::steady_clock::now();
	for (unsigned n = 0; n < N; ++n)
	{
		typename kernel_base<kernel>::test_input_t x(elem, xi1);
		typename kernel_base<kernel>::trial_input_t y(elem, xi2);
		k(x,y);
	}
	auto stop = std::chrono::steady_clock::now();
	auto diff = stop - start;
	std::cout << std::chrono::duration<double, std::micro>(diff).count() << " usec" << std::endl ;
}

int main(void)
{
	std::cout << N << " evaluations: " << std::endl;

	std::cout << "laplace G kernel with wall output:      | ";
	tester(laplace_3d_SLP_kernel());

	std::cout << "laplace H kernel with wall output:      | ";
	tester(laplace_3d_DLP_kernel());

	std::cout << "couple G H kernel with wall output:     | ";
	tester(create_couple_kernel(laplace_3d_DLP_kernel(), laplace_3d_SLP_kernel()));

	return 0;
}


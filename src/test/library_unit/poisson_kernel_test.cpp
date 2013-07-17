#include "library/poisson_kernel.hpp"
#include <iostream>
#include <chrono>

unsigned const N = 1e4;

template <class kernel>
void tester(kernel const &k)
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
		typename kernel::test_input_t x(elem, xi1);
		typename kernel::trial_input_t y(elem, xi2);
		k(x,y);
	}
	auto stop = std::chrono::steady_clock::now();
	auto diff = stop - start;
	std::cout << std::chrono::duration <double, std::micro> (diff).count() << " usec" << std::endl ;
}

int main(void)
{
	std::cout << N << " evaluations: " << std::endl;
	std::cout << "poisson G kernel with wall output:      | ";
	tester(poisson_G_kernel());
	std::cout << "poisson G kernel with immediate output: | ";
	tester(poisson_G_kernel_immediate());
	std::cout << "poisson H kernel with wall output:      | ";
	tester(poisson_H_kernel());
	return 0;
}


#include "library/poisson_kernel.hpp"
#include <iostream>
#include <chrono>

template <class kernel>
void tester(kernel const &k)
{
	typedef typename kernel::test_input_t test_input_t;
	typedef typename kernel::trial_input_t trial_input_t;

	quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;

	quad_1_elem elem(coords);
	auto xi1 = quad_1_elem::domain_t::get_center();
	auto xi2 = quad_1_elem::domain_t::get_corner(0);

	unsigned N = 10000;

	auto start = std::chrono::steady_clock::now();
	for (unsigned n = 0; n < N; ++n)
	{
		test_input_t x(elem, xi1);
		trial_input_t y(elem, xi2);
		k(x,y);
//		output_t out(x, y);
	}
	auto diff = std::chrono::steady_clock::now() - start;
	std::cout << N << " evaluations: " << std::chrono::duration <double, std::micro> (diff).count() << " usec" << std::endl ;
}

int main(void)
{
	std::cout << "poisson G kernel with wall output: ";
	tester(poisson_G_kernel());
	std::cout << "poisson G kernel with immediate output: ";
	tester(poisson_G_kernel_immediate());
	std::cout << "poisson H kernel with wall output: ";
	tester(poisson_H_kernel());
	return 0;
}


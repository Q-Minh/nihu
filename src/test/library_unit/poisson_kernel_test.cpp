#include "library/poisson_kernel.hpp"
#include <iostream>
#include <chrono>

template <class kernel>
void tester(void)
{
	typedef typename kernel::test_input_t test_input_t;
	typedef typename kernel::trial_input_t trial_input_t;
	typedef typename kernel::output_t output_t;

	quad_1_elem::coords_t coords;
	coords <<
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 0;

	quad_1_elem elem(coords);
	auto xi1 = quad_1_elem::domain_t::get_center();
	auto xi2 = quad_1_elem::domain_t::get_corner(0);

	unsigned N = 1000;

	auto start = std::chrono::steady_clock::now();
	for (unsigned n = 0; n < N; ++n)
	{
		test_input_t x(elem, xi1);
		trial_input_t y(elem, xi2);
		output_t out(x, y);
	}
	auto diff = std::chrono::steady_clock::now() - start;
	std::cout << std::chrono::duration <double, std::micro> (diff).count() << std::endl;
}

int main(void)
{
	tester<poisson_G_kernel>();
	tester<poisson_H_kernel>();
	return 0;
}


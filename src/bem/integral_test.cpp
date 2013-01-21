#include <iostream>
#include <algorithm>
#include <vector>

#include "integral.hpp"
#include "kernel.hpp"

static unsigned const NELEM = 100000;

template <class elem_t, class field_gen_t, class kernel_t, unsigned N>
struct tester
{
	typedef typename elem_t::coords_t coords_t;
	typedef field<elem_t, field_gen_t> field_t;
	typedef integral<field_t, kernel_t, N> integral_t;

	static void eval(void)
	{
		gauss_quad<typename elem_t::domain_t, N>::init();
		kernel_t::set_wave_number(1.0);

		std::vector<elem_t *> elements;
		for (unsigned i = 0; i < NELEM; ++i)
			elements.push_back(new elem_t(elem_t::coords_t::Random()));

		std::for_each(elements.begin(), elements.end(), 
			[] (elem_t const *e) {
				integral_t::eval(field_t(*e));
			}
		);

		std::for_each(elements.begin(), elements.end(), 
			[] (elem_t const *e) { delete e; }
		);
	}
};

int main(void)
{
	tester<quad_1_elem, constant_field, green_kernel, 5>::eval();
	return 0;
}


#include <iostream>
#include <algorithm>
#include <iterator>

#include "../bem/element.hpp"
#include "../bem/kernel.hpp"
#include "../bem/elem_accelerator.hpp"
#include "../bem/nset_accelerator.hpp"

int main(void)
{
	typedef quad_1_elem elem_t;
	typedef green_G_kernel kernel_t;

	typedef gauss_quadrature<elem_t::domain_t> quadrature_t;

	typedef elem_t::lset_t nset_t;
	typedef nset_accelerator<nset_t> accel_t;
	typedef nset_pool<nset_t> nset_pool_t;

	typedef kernel_t::input_t kernel_input_t;

	std::vector<quadrature_t> q;
	q.push_back(quadrature_t(1));
	q.push_back(quadrature_t(2));
	q.push_back(quadrature_t(3));

	nset_pool_t np(q.begin(), q.end());

	elem_t e(elem_t::coords_t::Random());
	typedef elem_accelerator<kernel_input_t> elem_accel_t;

	elem_pool<elem_t, kernel_input_t> ep(e, q.begin(), q.end());

	std::for_each(np.begin(), np.end(),
		[] (accel_t const &a) {
		std::copy(a.begin(), a.end(), std::ostream_iterator<nset_t::L_t>(std::cout, "\n\n"));
	});

	std::for_each(ep.begin(), ep.end(),
		[] (elem_accel_t const &a) {
		std::for_each(a.begin(), a.end(),
			[] (kernel_input_t const &k) {
			std::cout <<k.get_x() <<std::endl;
		});
	});

	return 0;
}


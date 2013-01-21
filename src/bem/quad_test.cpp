#include "quadrature.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"
#include <iostream>

template <class Q>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			Q::init();
			typename Q::xivec_t xi = Q::get_xi();
			typename Q::weightvec_t w = Q::get_weight();
			std::cout << "size: " << Q::size << std::endl;
			std::cout << "base points:" << std::endl << xi << std::endl;
			std::cout << "weights:" << std::endl << w << std::endl;
			std::cout << "sum of weights: " << w.sum() << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<
		gauss_quad<tria_domain, 1>,
		gauss_quad<tria_domain, 2>,
		gauss_quad<tria_domain, 3>,
		gauss_quad<tria_domain, 5>
	> quadSequence;

	tmp::call_each<quadSequence, tester<tmp::_1> >();

	return 0;
}


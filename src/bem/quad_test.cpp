#include "quadrature.hpp"
#include "../tmp/sequence.hpp"
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
	typedef tiny<
		gauss_quad<line_domain, 10>,
		gauss_quad<quad_domain, 4>
	> quadSequence;

	tmp::call_each<quadSequence, tester<_1> >();

	return 0;
}


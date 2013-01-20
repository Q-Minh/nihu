#include <iostream>
#include "element.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/control.hpp"

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			ElemType e(0, ElemType::nodes_t::Zero(), ElemType::coords_t::Random());
			typename ElemType::xi_t xi = ElemType::xi_t::Zero();
			typename ElemType::x_t x = e.get_x(xi);
			typename ElemType::x_t n = e.get_normal(xi);
			
			std::cout << "xi: " << xi << std::endl;
			std::cout << "x: " << x << std::endl;
			std::cout << "n: " << n << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<
		line_1_elem,
		tria_1_elem,
		parallelogram_elem,
		quad_1_elem
	> elemVector;
	tmp::call_each<elemVector, tester<tmp::_1> >();
	return 0;
}


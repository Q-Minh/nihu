#include <iostream>
#include "bem/element.hpp"
#include "tmp/sequence.hpp"
#include "tmp/control.hpp"

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl << "Elem type id: " << ElemType::id << std::endl;

			typename ElemType::xi_t xi = ElemType::domain_t::get_center();
			std::cout << "xi: " << xi.transpose() << std::endl;

			ElemType e(ElemType::coords_t::Random());
			std::cout << "coords: " << e.get_coords() << std::endl;

			typename ElemType::x_t c = e.get_center();
			std::cout << "center: " << c << std::endl;

			typename ElemType::x_t x = e.get_x(xi);
			std::cout << "x(xi): " << x << std::endl;

			typename ElemType::x_t n = e.get_normal(xi);
			std::cout << "n(xi): " << n << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<tria_1_elem, quad_1_elem, tria_2_elem, quad_2_elem>,
		tester<tmp::_1>
	>();

	return 0;
}


#include "bem/element.hpp"
#include "tmp/sequence.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

#include <iostream>

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl;
			std::cout << "Elem type ID: " << ElemType::id << std::endl;
			std::cout << "===================" << std::endl;

			typename ElemType::xi_t xi = ElemType::domain_t::get_center();
			std::cout << "xi center: " << xi.transpose() << std::endl;

			ElemType e(ElemType::coords_t::Random());
			
			typename ElemType::x_t c = e.get_center();
			std::cout << "x  center: " << c.transpose() << std::endl;
			std::cout << "nodal coords: " << std::endl << e.get_coords() << std::endl;
			
			typename ElemType::x_t x = e.get_x(xi);
			std::cout << "x(xi): " << x.transpose() << std::endl;

			typename ElemType::x_t n = e.get_normal(xi);
			std::cout << "n(xi): " << n.transpose() << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			tria_1_elem, quad_1_elem,
			tria_2_elem, quad_2_elem
		>,
		tester<tmp::_1>
	>();

	return 0;
}


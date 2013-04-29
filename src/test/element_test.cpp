#include <iostream>
#include "../bem/element.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/control.hpp"

template <class ElemType>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			typename ElemType::xi_t xi = ElemType::domain_t::get_center();
			std::cout << "xi: " << xi << std::endl;

			ElemType e(ElemType::coords_t::Random());
			std::cout << "coords: " << e.get_coords() << std::endl;

			typename ElemType::dx_t dx = e.get_dx(xi);
			std::cout << "dx: " << dx << std::endl;

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
	tester<quad_2_elem>::type t;
	t();

	return 0;
}

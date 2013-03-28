#include <iostream>

#include "elem_accelerator.hpp"

#include "../tmp/sequence.hpp"
#include "../tmp/control.hpp"

template <class ElemType>
struct tester
{
	struct type
	{
		typedef ElemType elem_t;
		typedef location<typename elem_t::x_t> input_t;
		typedef elem_accelerator<input_t> elem_accelerator_t;
		typedef gauss_quadrature<typename elem_t::domain_t> quadrature_t;

		void operator() (void)
		{
			ElemType e(ElemType::coords_t::Random());
			quadrature_t q(2);

			elem_accelerator_t ea(e, q);

			std::for_each(
				ea.begin(),
				ea.end(),
				[] (input_t const &a) { std::cout << a.get_x() << std::endl; }
			);
		}
	};
};

int main(void)
{
	typedef tmp::vector<
		tria_1_elem,
		quad_1_elem
	> elemVector;

	tmp::call_each<elemVector, tester<tmp::_1> >();

	return 0;
}

